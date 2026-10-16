#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The crystal's frequency error from carrier Doppler alone, one satellite at a time: the
known-position method, checked against the PVT solution's drift on the same data.

  python tools/doppler_check.py run.obs --eph eph.json --fix-file fix.json [--csv run.csv]

For a receiver that is not moving, every satellite's measured Doppler is the geometric Doppler
(from its orbit and the receiver's position) plus one common offset: the crystal's error, which
shifts the LO by eps * f_L1. So for each satellite,

    eps = -(D_measured - D_predicted) / f_L1          (same sign convention as tools/wwv.py)

and every satellite should give the SAME eps. Agreement between satellites is the test of the
method; agreement with the PVT drift (--csv, the receiver's 30 s fit) is the test of its accuracy.
Unlike the PVT drift this needs no fix: one tracked satellite, an ephemeris for it and a position
good to a few km (a km of position error is ~0.1 Hz of Doppler, under 0.1 ppb).

Inputs are the RINEX observations tools/clocklog.py writes with --rinex-file (D1C column), the
receiver's ephemeris store (--eph-file of the same run) and the receiver position from the fix
file. The position is read and never printed."""
import argparse
import csv
import datetime
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))
from gpsrx import pvt  # noqa: E402

F_L1 = 1575.42e6
LAMBDA = pvt.C / F_L1


def read_obs_tow(path):
    """[(tow, [(sys, prn, D1C, S1C)])]: the epoch's GPS time of week from its date line."""
    epochs = []
    with open(path) as f:
        for line in f:
            if line.rstrip().endswith("END OF HEADER"):
                break
        cur = None
        for line in f:
            if line.startswith(">"):
                p = line.split()
                d = datetime.date(int(p[1]), int(p[2]), int(p[3]))
                dow = (d.weekday() + 1) % 7                       # GPS week starts on Sunday
                tow = dow * 86400 + int(p[4]) * 3600 + int(p[5]) * 60 + float(p[6])
                cur = (tow, [])
                epochs.append(cur)
            elif cur is not None and line[:1] == "G":
                vals = [float(line[3 + 16 * i: 3 + 16 * i + 14]) for i in range(4)]
                cur[1].append((line[0], int(line[1:3]), vals[2], vals[3]))
    return epochs


def read_obs_csv(path, utc_to_gps_s):
    """clocklog.py --obs-file rows -> [(tow, [(sys, prn, carrier_hz, cn0)])], one epoch per wall second.
    The wall clock (NTP) stands in for GPS time: a second of error is under 0.5 ppb of Doppler."""
    by_sec = {}
    for r in csv.DictReader(open(path)):
        u = datetime.datetime.fromisoformat(r["utc"])
        dow = (u.weekday() + 1) % 7
        sow = dow * 86400 + u.hour * 3600 + u.minute * 60 + u.second
        if r["sys"] != "GPS" or float(r["lock"]) < 0.5:
            continue
        by_sec.setdefault(sow, {})[int(r["prn"])] = (float(r["carrier_hz"]), float(r["cn0_db"]))
    return [(sow + utc_to_gps_s, [("G", prn, d, c) for prn, (d, c) in sorted(v.items())]) for sow, v in sorted(by_sec.items())]


def predicted_doppler(eph, t_rx, rx_ecef):
    """Geometric Doppler (Hz, positive = closing) of a satellite seen from a fixed receiver at
    receive time t_rx: range rate from the orbit at transmit time, by central difference."""
    r = np.linalg.norm(pvt.sat_ecef(eph, t_rx) - rx_ecef)
    t_tx = t_rx - r / pvt.C
    h = 0.5
    p0, p1 = pvt.sat_ecef(eph, t_tx - h), pvt.sat_ecef(eph, t_tx + h)
    rate = (np.linalg.norm(p1 - rx_ecef) - np.linalg.norm(p0 - rx_ecef)) / (2 * h)
    return -rate / LAMBDA - eph.get("af1", 0.0) * F_L1      # the SV clock's own rate, tiny


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("obs", help="RINEX observations from clocklog.py --rinex-file, or the raw CSV from --obs-file")
    ap.add_argument("--utc-to-gps", type=float, default=16.9,
                    help="for the raw CSV: GPS time minus the row's wall time (18 s of leap seconds less ~1.1 s of receiver latency, measured)")
    ap.add_argument("--eph", required=True, help="the run's ephemeris store (--eph-file)")
    ap.add_argument("--fix-file", required=True, help="fix json: the receiver position (not printed)")
    ap.add_argument("--csv", default="", help="the run's clocklog CSV: compare with the PVT drift")
    ap.add_argument("--min-cn0", type=float, default=30.0)
    a = ap.parse_args()

    ephs = {int(v["prn"]): v for k, v in json.load(open(a.eph)).items() if k.startswith("G")}
    fx = json.load(open(a.fix_file))
    fx = fx.get("fix", fx)
    rx = np.array(fx.get("ecef_mean") or fx["ecef"], float)
    epochs = read_obs_csv(a.obs, a.utc_to_gps) if a.obs.endswith(".csv") else read_obs_tow(a.obs)
    if not epochs:
        sys.exit("no epochs")

    per = {}                                    # prn -> [(tow, eps_ppb)]
    for tow, obs in epochs:
        for sysc, prn, d_meas, cn0 in obs:
            if prn not in ephs or cn0 < a.min_cn0:
                continue
            d_pred = predicted_doppler(ephs[prn], tow, rx)
            per.setdefault(prn, []).append((tow, -(d_meas - d_pred) / F_L1 * 1e9))

    print(f"{len(epochs)} epochs, {len(per)} satellites with an ephemeris and C/N0 >= {a.min_cn0:.0f} dB-Hz")
    print("| PRN | epochs | crystal error from this satellite ALONE (ppb) | sd over epochs |\n|---:|---:|---:|---:|")
    means = []
    for prn in sorted(per):
        e = np.array([v for _, v in per[prn]])
        means.append(e.mean())
        print(f"| G{prn:02d} | {len(e)} | {e.mean():+.1f} | {e.std(ddof=1) if len(e) > 1 else 0:.1f} |")
    means = np.array(means)
    print(f"\nall satellites: mean {means.mean():+.1f} ppb, satellite-to-satellite sd {means.std(ddof=1) if len(means) > 1 else 0:.1f} ppb")
    # per epoch, the mean over satellites: the one-number-per-second estimate
    by_epoch = {}
    for prn, vals in per.items():
        for tow, v in vals:
            by_epoch.setdefault(tow, []).append(v)
    ep = np.array([np.mean(v) for _, v in sorted(by_epoch.items())])
    print(f"per-epoch mean over satellites: {ep.mean():+.1f} ppb, epoch-to-epoch sd {ep.std(ddof=1):.1f} ppb (n={len(ep)})")

    if a.csv:
        rows = [r for r in csv.DictReader(open(a.csv)) if r["fix_valid"] == "1" and r["clock_drift_ppb"]]
        tows = {round(float(r["gps_tow"])) for r in rows}
        inside = [float(r["clock_drift_ppb"]) for r in rows if round(float(r["gps_tow"])) in {round(t) for t in by_epoch}]
        if inside:
            print(f"PVT drift over the same epochs: {np.mean(inside):+.1f} ppb (sd {np.std(inside, ddof=1):.1f}); "
                  f"Doppler-only minus PVT: {ep.mean() - np.mean(inside):+.1f} ppb")
        else:
            print(f"CSV has {len(tows)} valid epochs but none overlap the observations")
    return 0


if __name__ == "__main__":
    sys.exit(main())
