#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The radio's frequency error from a standard-frequency station's carrier (WWV 5/10/15 MHz, CHU
7.85 MHz): an independent check on the GPS number, on the LO path at HF.

  python tools/wwv.py --freq 10e6 --seconds 120 --antenna "Antenna A"
  python tools/wwv.py --freq 10e6 --seconds 120 --antenna "Antenna A" --corr-ppm -0.75   # with the correction

Tunes IF_HZ below the carrier (away from the zero-IF DC spike), mixes the carrier to 0 Hz on the
NOMINAL sample grid, averages down to 1 kS/s in memory (no IQ is written anywhere), then one long
Hann-windowed FFT finds the carrier within +-SEARCH_HZ. With the crystal off by eps, the LO and the
sample clock are both off by eps, and the carrier lands at -eps * f_carrier whatever the IF:

    eps_ppb = -offset_hz / f_carrier * 1e9

Skywave adds ionospheric Doppler (~0.1 Hz at 10 MHz, i.e. ~10 ppb, varying over minutes), so the
10 s slices are reported too: their spread is the measurement's honest uncertainty.

The bias-T is written OFF and read back first: this is for HF antennas that must not see DC."""
import argparse
import sys
import time

import numpy as np
import SoapySDR
from SoapySDR import SOAPY_SDR_CF32, SOAPY_SDR_RX

IF_HZ = 25e3
OUT_RATE = 1000.0
SEARCH_HZ = 40.0
SLICE_S = 10.0


def capture(sdr, fs, seconds, shift_hz=0.0):
    """1 kS/s complex baseband of the carrier region, mixed on the nominal sample grid; shift_hz is
    an extra digital shift (the software correction)."""
    dec = int(round(fs / OUT_RATE))
    st = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CF32)
    buf = np.empty(65536, np.complex64)
    need = int(seconds * fs)
    n, carry, out, overflows = 0, np.empty(0, np.complex128), [], 0
    w = -2j * np.pi * IF_HZ / fs
    ws = -2j * np.pi * shift_hz / fs
    sdr.activateStream(st)
    t0 = time.time()
    try:
        while n < need:
            r = sdr.readStream(st, [buf], len(buf), timeoutUs=1000000)
            if r.ret < 0:
                overflows += int(r.ret == SoapySDR.SOAPY_SDR_OVERFLOW)
                if r.ret != SoapySDR.SOAPY_SDR_OVERFLOW and time.time() - t0 > 5:
                    raise RuntimeError(f"readStream error {r.ret}")
                continue
            k = np.arange(n, n + r.ret)
            x = buf[:r.ret].astype(np.complex128) * np.exp(w * (k % int(fs)))   # IF_HZ * fs-periodic: exact
            if shift_hz:
                x *= np.exp(ws * k)
            n += r.ret
            x = np.concatenate([carry, x])
            m = len(x) // dec * dec
            out.append(x[:m].reshape(-1, dec).mean(axis=1))
            carry = x[m:]
    finally:
        sdr.deactivateStream(st)
        sdr.closeStream(st)
    return np.concatenate(out), overflows, time.time() - t0


def peak(z, rate):
    """Carrier frequency (Hz) within +-SEARCH_HZ and its SNR (dB over the median bin)."""
    N = len(z)
    nfft = 1 << int(np.ceil(np.log2(N * 4)))
    Z = np.abs(np.fft.fftshift(np.fft.fft(z * np.hanning(N), nfft))) ** 2
    f = np.fft.fftshift(np.fft.fftfreq(nfft, 1.0 / rate))
    win = np.abs(f) <= SEARCH_HZ
    i = np.flatnonzero(win)[np.argmax(Z[win])]
    a, b, c = np.log(Z[i - 1:i + 2] + 1e-300)
    d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
    snr = 10 * np.log10(Z[i] / np.median(Z[win]))
    return float(f[i] + d * (f[1] - f[0])), float(snr)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--freq", type=float, default=10e6)
    ap.add_argument("--seconds", type=float, default=120.0)
    ap.add_argument("--antenna", default="Antenna A")
    ap.add_argument("--rate", type=float, default=2e6)
    ap.add_argument("--gain", type=float, default=-1, help="dB; -1 = AGC")
    ap.add_argument("--corr-ppm", type=float, default=0.0, help="driver frequency correction (ppm) to apply")
    ap.add_argument("--soft-ppb", type=float, default=0.0,
                    help="software correction: the clock error (ppb) to remove, as a digital shift of soft_ppb*1e-9*freq Hz")
    ap.add_argument("--driver", default="sdrplay")
    a = ap.parse_args()

    sdr = SoapySDR.Device(f"driver={a.driver}")
    sdr.writeSetting("biasT_ctrl", "false")
    bt = str(sdr.readSetting("biasT_ctrl")).strip().lower()
    print(f"bias-T readback: biasT_ctrl = {bt!r}", flush=True)
    if bt != "false":
        print("bias-T did NOT read back off - stopping before anything streams", flush=True)
        return 2
    sdr.setAntenna(SOAPY_SDR_RX, 0, a.antenna)
    sdr.setSampleRate(SOAPY_SDR_RX, 0, a.rate)
    sdr.setBandwidth(SOAPY_SDR_RX, 0, 200e3)
    if a.corr_ppm:
        sdr.setFrequencyCorrection(SOAPY_SDR_RX, 0, a.corr_ppm)
    sdr.setFrequency(SOAPY_SDR_RX, 0, a.freq - IF_HZ)
    if a.gain < 0:
        sdr.setGainMode(SOAPY_SDR_RX, 0, True)
    else:
        sdr.setGainMode(SOAPY_SDR_RX, 0, False)
        sdr.setGain(SOAPY_SDR_RX, 0, a.gain)
    print(f"antenna {sdr.getAntenna(SOAPY_SDR_RX, 0)}; tuned {sdr.getFrequency(SOAPY_SDR_RX, 0) / 1e6:.6f} MHz "
          f"(carrier {a.freq / 1e6:.4f} MHz at +{IF_HZ / 1e3:.0f} kHz); rate {sdr.getSampleRate(SOAPY_SDR_RX, 0) / 1e6:.3f} MS/s; "
          f"correction {sdr.getFrequencyCorrection(SOAPY_SDR_RX, 0):+.4f} ppm; {a.seconds:.0f} s", flush=True)

    expect_hz = -a.soft_ppb * 1e-9 * a.freq          # where that clock error puts the carrier; shifted out below
    if a.soft_ppb:
        print(f"software correction: {a.soft_ppb:+.1f} ppb -> digital shift of {-expect_hz:+.4f} Hz", flush=True)
    z, ovf, wall = capture(sdr, a.rate, a.seconds, shift_hz=expect_hz)
    if a.corr_ppm:
        sdr.setFrequencyCorrection(SOAPY_SDR_RX, 0, 0.0)
    f_all, snr_all = peak(z, OUT_RATE)
    L = int(SLICE_S * OUT_RATE)
    sl = [peak(z[k:k + L], OUT_RATE) for k in range(0, len(z) - L + 1, L)]
    fs_ = np.array([s[0] for s in sl])
    # a slice in a fade can lock onto a sideband or a spur: keep slices with SNR >= 15 dB that sit
    # within 0.5 Hz of the median (0.5 Hz = 50 ppb at 10 MHz, far beyond any ionospheric Doppler)
    good = np.array([s[1] for s in sl]) >= 15.0
    if good.any():
        good &= np.abs(fs_ - np.median(fs_[good])) <= 0.5
    print(f"captured {len(z) / OUT_RATE:.1f} s in {wall:.1f} s wall; overflows {ovf}")
    print(f"whole record: carrier offset {f_all:+.4f} Hz, SNR {snr_all:.1f} dB  ->  {-f_all / a.freq * 1e9:+.1f} ppb")
    print(f"{SLICE_S:.0f} s slices (Hz; ? = rejected): " + " ".join(f"{f:+.3f}{'' if g else '?'}" for f, g in zip(fs_, good)))
    if good.sum() >= 3:
        m, s = fs_[good].mean(), fs_[good].std(ddof=1)
        print(f"kept slices: {good.sum()}/{len(sl)}, mean {m:+.4f} Hz, sd {s:.4f} Hz  ->  "
              f"{-m / a.freq * 1e9:+.1f} +- {s / np.sqrt(good.sum()) / a.freq * 1e9:.1f} ppb "
              f"(slice-to-slice sd {s / a.freq * 1e9:.1f} ppb)")
    if snr_all < 15.0:
        print("WEAK: whole-record SNR under 15 dB - the peak may not be the carrier")
    return 0


if __name__ == "__main__":
    sys.exit(main())
