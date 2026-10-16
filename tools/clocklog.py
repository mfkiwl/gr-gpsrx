#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The radio's sample clock against GPS, once per second, to a CSV: a software GPSDO's logbook.

  python tools/clocklog.py --csv /somewhere/private/clock.csv --seconds 3600 \\
      --fix-file /somewhere/private/fix.json

Radio -> Receiver (C++ channels), as apps/gpsrx_live.py, plus a sink on the 'fix' port that writes
one row per fix (the receiver solves once per stream second, so rows are 1 s apart in GPS time) and,
before the first fix or during an outage of more than 2 s, a 'no fix' row per wall second.

Columns: utc (wall clock, ISO), gps_tow (s, valid fixes), n_sats, fix_valid, clock_drift_ppb (the
receiver's 30 s straight-line fit, so neighbouring rows share most of their data), residual_rms_m,
clock_offset_s (sample/fs - GPS time at the fix: the clock's PHASE, what an Allan deviation should
be computed from), resync (1 when the receiver saw lost samples and restarted its clock history).
No position goes in the CSV; the position goes only to --fix-file.

The bias-T is written on and READ BACK before the flowgraph starts; if the readback is not 'true'
the run stops there. (The readback is the driver's setting, not a current measurement.)"""
import argparse
import csv
import datetime
import os
import sys
import threading
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "examples"))
import pmt  # noqa: E402
from gnuradio import gr, soapy  # noqa: E402
try:
    from gnuradio import gpsrx
except ImportError:
    import _devpath  # noqa: E402,F401
    from gnuradio import gpsrx

COLUMNS = ["utc", "gps_tow", "n_sats", "fix_valid", "clock_drift_ppb", "residual_rms_m", "clock_offset_s", "resync"]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="milliseconds")


class clock_sink(gr.basic_block):
    """Writes the clock fields of every fix to a CSV, and nothing about where the antenna is."""

    def __init__(self, path, obs_path=""):
        gr.basic_block.__init__(self, name="clock_sink", in_sig=None, out_sig=None)
        self.message_port_register_in(pmt.intern("fix"))
        self.set_msg_handler(pmt.intern("fix"), self.on_fix)
        self.message_port_register_in(pmt.intern("obs"))
        self.set_msg_handler(pmt.intern("obs"), self.on_obs)
        # raw per-channel observables (no fix needed): for the one-satellite method, timed by the wall clock
        self._obs = None
        if obs_path:
            self._obs = open(obs_path, "w", newline="")
            self._obs_w = csv.writer(self._obs)
            self._obs_w.writerow(["utc", "slot", "prn", "sys", "carrier_hz", "cn0_db", "lock"])
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        self._fh = open(path, "w", newline="")
        self._w = csv.DictWriter(self._fh, fieldnames=COLUMNS)
        self._w.writeheader()
        self._lock = threading.Lock()
        self.last_fix_wall = time.time()
        self.n_rows = self.n_valid = 0

    def write(self, row):
        with self._lock:
            self._w.writerow({k: row.get(k, "") for k in COLUMNS})
            self._fh.flush()
            self.n_rows += 1
            self.n_valid += int(row.get("fix_valid", 0))

    def on_fix(self, msg):
        f = pmt.to_python(msg)
        if not isinstance(f, dict):
            return
        self.last_fix_wall = time.time()
        valid = bool(f.get("ok") and f.get("valid"))
        drift = f.get("clock_drift_ppb")
        self.write(dict(
            utc=utc_now(),
            gps_tow=f"{f['gps_tow']:.9f}" if valid and f.get("gps_tow") is not None else "",
            n_sats=int(f.get("n", 0)),
            fix_valid=int(valid),
            clock_drift_ppb=f"{drift:.3f}" if valid and drift is not None else "",
            residual_rms_m=f"{f['rms_m']:.3f}" if f.get("rms_m") is not None else "",
            clock_offset_s=f"{f['clock_offset_s']:.12f}" if valid and f.get("clock_offset_s") is not None else "",
            resync=int(bool(f.get("resync")))))

    def on_obs(self, msg):
        if self._obs is None:
            return
        d = pmt.to_python(msg)
        if not isinstance(d, dict) or "carrier_hz" not in d:
            return
        with self._lock:
            self._obs_w.writerow([utc_now(), int(d.get("slot", -1)), int(d.get("prn", 0)), str(d.get("sys", "GPS")),
                                  f"{float(d['carrier_hz']):.3f}", f"{float(d.get('cn0_db', 0)):.1f}", f"{float(d.get('lock', 0)):.2f}"])
            self._obs.flush()

    def no_fix_row(self, n_tracking):
        self.write(dict(utc=utc_now(), n_sats=n_tracking, fix_valid=0, resync=0))

    def close(self):
        with self._lock:
            self._fh.close()
            if self._obs is not None:
                self._obs.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="the log; keep it outside the repo")
    ap.add_argument("--seconds", type=float, default=60.0)
    ap.add_argument("--driver", default="sdrplay")
    ap.add_argument("--dev-args", default="")
    ap.add_argument("--antenna", default="Antenna B", help="the RSPdx's bias-T powers Antenna B only")
    ap.add_argument("--rate", type=float, default=2.048e6)
    ap.add_argument("--gain", type=float, default=40.0)
    ap.add_argument("--agc", type=int, default=1)
    ap.add_argument("--no-bias-t", action="store_true", help="passive antenna: skip the bias-T (and its check)")
    ap.add_argument("--channels", type=int, default=8)
    ap.add_argument("--fix-file", default=os.path.join(HERE, "..", "lab_local", "clocklog_fix.json"))
    ap.add_argument("--eph-file", default=os.path.join(HERE, "..", "lab_local", "eph_live.json"))
    ap.add_argument("--iono-file", default=os.path.join(HERE, "..", "lab_local", "iono_terms.json"))
    ap.add_argument("--obs-file", default="", help="also log raw per-channel Doppler (utc, prn, carrier_hz, C/N0) here - works with NO fix; a private file")
    ap.add_argument("--rinex-file", default="", help="also write RINEX 3 observations (per-satellite Doppler) here; a private file")
    ap.add_argument("--every", type=float, default=10.0, help="console progress line every N s")
    a = ap.parse_args()

    tb = gr.top_block()
    src = soapy.source(f"driver={a.driver}" + (f",{a.dev_args}" if a.dev_args else ""), "fc32", 1, "", "", [""], [""])
    src.set_sample_rate(0, a.rate)
    src.set_bandwidth(0, 2.5e6)
    if a.antenna:
        src.set_antenna(0, a.antenna)
    src.set_frequency(0, 1575.42e6)
    src.set_gain_mode(0, bool(a.agc))
    src.set_gain(0, a.gain)
    src.set_min_output_buffer(int(16 * a.rate * 4e-3))
    if not a.no_bias_t:
        if a.driver == "sdrplay" and a.antenna != "Antenna B":
            print(f"bias-T requested on {a.antenna or 'the default port'}: the RSPdx powers Antenna B only - stopping", flush=True)
            return 2
        src.write_setting("biasT_ctrl", "true")
        readback = str(src.read_setting("biasT_ctrl")).strip().lower()
        print(f"bias-T readback: biasT_ctrl = {readback!r}", flush=True)
        if readback != "true":
            print("bias-T did NOT read back on - stopping before the run", flush=True)
            return 2

    rx = gpsrx.receiver(a.rate, n_channels=a.channels, fix_file=a.fix_file, eph_file=a.eph_file,
                        iono_file=a.iono_file, hold=False, engine="cpp", rinex_file=a.rinex_file)
    st = gpsrx.status_sink(every_s=3600.0)          # for the tracking count only; its table is not printed
    sink = clock_sink(a.csv, obs_path=a.obs_file)
    tb.connect(src, rx)
    for port in ("sky", "status", "obs", "nav", "fix"):
        tb.msg_connect(rx, port, st, port)
    tb.msg_connect(rx, "fix", sink, "fix")
    tb.msg_connect(rx, "obs", sink, "obs")

    t0 = time.time()
    next_report, next_tick = t0 + a.every, t0 + 1.0
    print(f"clocklog: {a.seconds:.0f} s -> {a.csv}; position -> {a.fix_file} (not printed)", flush=True)
    tb.start()
    try:
        while time.time() - t0 < a.seconds:
            time.sleep(0.05)
            now = time.time()
            # fixes arrive with +-100 ms of scheduler jitter, so one wall second can bring none and the
            # next two: only a 2 s silence is an outage, then a 'no fix' row for each second of it
            if now >= next_tick:
                next_tick += 1.0
                if now - sink.last_fix_wall > 2.0:
                    sink.no_fix_row(sum(1 for c in st.chan.values() if c.get("what") == "tracking"))
            if now >= next_report:
                next_report += a.every
                print(f"  {now - t0:5.0f} s: {sink.n_rows} rows, {sink.n_valid} valid fixes", flush=True)
    except KeyboardInterrupt:
        pass
    tb.stop()
    tb.wait()
    if not a.no_bias_t:
        try:
            src.write_setting("biasT_ctrl", "false")
        except Exception:
            pass
    sink.close()
    print(f"clocklog: {sink.n_rows} rows, {sink.n_valid} valid fixes in {time.time() - t0:.0f} s", flush=True)
    return 0 if sink.n_valid else 1


if __name__ == "__main__":
    sys.exit(main())
