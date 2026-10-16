#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-or-later
"""The sample clock's frequency error from a tools/clocklog.py CSV: the number, its stability, and
whether it has settled.

  python tools/allan.py clock.csv --out /somewhere/report_dir             # baseline: stats + ADEV
  python tools/allan.py cold.csv --out dir --warmup                       # + time to settle within 10 ppb
  python tools/allan.py clock.csv --out dir --temps room.csv              # + ppb per degree C

Everything is computed from the clock PHASE (clock_offset_s: sample/fs - GPS time, one per GPS
second), whose slope is the fractional frequency error. The receiver's clock_drift_ppb column is a
30 s straight-line fit of the same phase, so neighbouring values are not independent: it is
summarised (mean, sd, min, max), but the Allan deviation and the trend come from the phase.

Lost samples restart the receiver's clock history (resync = 1) and step the phase, so the valid
fixes are cut into segments wherever GPS time does not advance by one second or a resync shows,
and the Allan deviation pools its second differences over the segments: a step is never counted
as clock noise. The output names no position; the CSV holds none."""
import argparse
import csv
import datetime
import json
import os
import sys

import numpy as np

TAUS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]
MIN_INDEP = 5                 # a tau needs ~5 independent intervals (terms / tau) to name the floor
BLOCK_S = 300                 # frequency blocks for the trend
DRIFT_FIT_S = 30              # the receiver's drift fit window: earlier values are partial fits
SETTLE_PPB = 10.0             # warm-up: within this of the final value, and staying there


def load(path):
    rows = list(csv.DictReader(open(path)))
    if not rows:
        sys.exit(f"{path}: empty")
    t_wall = np.array([datetime.datetime.fromisoformat(r["utc"]).timestamp() for r in rows])
    valid = np.array([r["fix_valid"] == "1" for r in rows])
    return rows, t_wall, valid


def segments(rows, valid):
    """[(tow, phase, drift, wall index)] on an exact 1 s GPS-time grid, one per unbroken run.

    The receiver solves on a code epoch, so when the satellite set changes the solve instant moves
    by tens of ms and a step reads 0.85 or 1.18 s: not lost data. Each grid second takes the
    nearest fix within 0.5 s, carried to the grid instant along the run's straight-line rate (the
    carry is <= 0.5 s x a few ppb of rate error: < 1 ns, under the 1 s noise). Seconds with no fix
    are NaN holes. A run breaks only at a resync (samples lost: the phase stepped) or a gap > 5 s."""
    runs, cur = [], []
    for i, r in enumerate(rows):
        if r.get("resync") == "1":
            if cur:
                runs.append(cur)
            cur = []
            continue
        if not valid[i] or not r["clock_offset_s"]:
            continue
        tow, x = float(r["gps_tow"]), float(r["clock_offset_s"])
        if cur:
            step = tow - cur[-1][0]
            if step < -302400:                       # a week roll
                step += 604800
                tow += 604800
            if step <= 0 or step > 5.0:
                runs.append(cur)
                cur = []
        d = float(r["clock_drift_ppb"]) if r["clock_drift_ppb"] else np.nan
        cur.append((tow, x, d, i))
    if cur:
        runs.append(cur)
    out = []
    for run in runs:
        a = np.array(run)
        tow, x = a[:, 0], a[:, 1]
        if len(tow) < 3:
            continue
        rate = np.polyfit(tow - tow[0], x, 1)[0]
        grid = tow[0] + np.arange(int(np.floor(tow[-1] - tow[0])) + 1)
        j = np.clip(np.searchsorted(tow, grid), 1, len(tow) - 1)
        j = np.where(np.abs(tow[j - 1] - grid) <= np.abs(tow[j] - grid), j - 1, j)
        near = np.abs(tow[j] - grid) <= 0.5
        xg = np.where(near, x[j] + rate * (grid - tow[j]), np.nan)
        dg = np.where(near, a[j, 2], np.nan)
        out.append((grid, xg, dg, a[j, 3].astype(int)))
    return out


def oadev(segs, taus):
    """Overlapping Allan deviation from phase (s), tau0 = 1 s, pooled over segments.
    sigma^2(m) = sum (x[i+2m] - 2 x[i+m] + x[i])^2 / (2 m^2 tau0^2 N_terms)."""
    res = []
    for m in taus:
        ss, n = 0.0, 0
        for tow, x, _, _ in segs:
            if len(x) > 2 * m:
                d2 = x[2 * m:] - 2 * x[m:-m] + x[:-2 * m]
                d2 = d2[~np.isnan(d2)]                   # a hole in any of the three: that term is skipped
                ss += float(np.sum(d2 ** 2))
                n += len(d2)
        if n == 0:
            res.append(dict(tau=m, adev=None, n=0))
            continue
        adev = np.sqrt(ss / (2.0 * m * m * n))
        # approximate 1-sigma: the overlapping estimate has about n/m independent terms
        err = adev / np.sqrt(max(n / m, 1.0))
        res.append(dict(tau=m, adev=float(adev), err=float(err), n=int(n)))
    return res


def phase_freq(tow, x):
    """Fractional frequency (ppb) and its formal 1-sigma from a straight line through the phase."""
    ok = ~np.isnan(x)
    tow, x = tow[ok], x[ok]
    p, cov = np.polyfit(tow - tow[0], x, 1, cov=True)
    return p[0] * 1e9, np.sqrt(cov[0, 0]) * 1e9


def blocks(segs, t0_wall, t_wall):
    """Frequency (ppb) of each BLOCK_S block of contiguous phase, by GPS-time start."""
    out = []
    for tow, x, _, idx in segs:
        for k in range(0, len(x) - BLOCK_S + 1, BLOCK_S):
            sl = slice(k, k + BLOCK_S)
            f, _ = phase_freq(tow[sl], x[sl])
            out.append(dict(t_min=(t_wall[idx[k]] - t0_wall) / 60.0, f_ppb=float(f)))
    return out


def settled_drift(segs):
    """The receiver's drift values once its 30 s fit is full (after each start / resync)."""
    vals, tt = [], []
    for tow, _, d, idx in segs:
        keep = (tow - tow[0]) >= DRIFT_FIT_S
        keep &= ~np.isnan(d)
        vals.append(d[keep])
        tt.append(idx[keep])
    return np.concatenate(vals) if vals else np.array([]), np.concatenate(tt) if tt else np.array([], int)


def trend(blk, adev_at):
    """Slope of the last 30 min of block frequencies: ppb per hour, its 1-sigma, and a verdict."""
    tail = [b for b in blk if b["t_min"] >= blk[-1]["t_min"] - 30.0 + 1e-9] if blk else []
    if len(tail) < 3:
        return dict(slope_ppb_per_h=None, verdict="too short to tell (needs 3 blocks of 5 min in the last 30 min)")
    t = np.array([b["t_min"] for b in tail]) / 60.0
    f = np.array([b["f_ppb"] for b in tail])
    p, cov = np.polyfit(t, f, 1, cov=True)
    slope, se = float(p[0]), float(np.sqrt(cov[0, 0]))
    moved = abs(slope) * (t[-1] - t[0])
    noise = (adev_at or 0.0) * 1e9
    still = abs(slope) > 2 * se and moved > noise
    return dict(slope_ppb_per_h=slope, slope_se=se, moved_ppb=float(moved), noise_ppb=noise,
                verdict=("STILL TRENDING" if still else "settled") +
                f": {slope:+.1f} +- {se:.1f} ppb/h over the last {len(tail)} blocks "
                f"(moved {moved:.1f} ppb; ADEV(500 s, or 200 s) = {noise:.1f} ppb)")


def warmup(drift, idx, t_wall, t0_wall):
    """Time from the start of the log (or --t0, the plug-in) until the drift enters, and stays within, SETTLE_PPB of its
    final value (the mean of the last 10 minutes); and how far it moved on the way."""
    if len(drift) < 600:
        return dict(verdict="too short: needs at least 10 minutes of settled drift")
    t_min = (t_wall[idx] - t0_wall) / 60.0
    final = float(np.mean(drift[t_min >= t_min[-1] - 10.0]))
    outside = np.where(np.abs(drift - final) > SETTLE_PPB)[0]
    if len(outside) == 0:
        settle = float(t_min[0])
    elif outside[-1] + 1 < len(drift):
        settle = float(t_min[outside[-1] + 1])
    else:
        settle = None
    return dict(final_ppb=final, first_ppb=float(drift[0]), first_fix_min=float(t_min[0]),
                settle_min=settle, moved_ppb=float(final - drift[0]),
                span_ppb=float(np.max(drift) - np.min(drift)),
                verdict=(f"settled within {SETTLE_PPB:.0f} ppb of {final:.1f} ppb at {settle:.1f} min after the start"
                         if settle is not None else "never settled within the run"))


def temp_fit(drift, idx, t_wall, temps_path):
    """Drift against room temperature: temps CSV of utc,temp_c, interpolated to each fix."""
    tr = list(csv.DictReader(open(temps_path)))
    tt = np.array([datetime.datetime.fromisoformat(r["utc"]).timestamp() for r in tr])
    tc = np.array([float(r["temp_c"]) for r in tr])
    o = np.argsort(tt)
    tt, tc = tt[o], tc[o]
    tw = t_wall[idx]
    inside = (tw >= tt[0]) & (tw <= tt[-1])
    if inside.sum() < 60 or np.ptp(tc) < 0.5:
        return dict(verdict="not enough temperature range or overlap to fit (needs >= 0.5 C and 60 fixes)")
    T = np.interp(tw[inside], tt, tc)
    p, cov = np.polyfit(T, drift[inside], 1, cov=True)
    return dict(ppb_per_c=float(p[0]), se=float(np.sqrt(cov[0, 0])), temp_range_c=float(np.ptp(tc)),
                n=int(inside.sum()), verdict=f"{p[0]:+.1f} ppb/C (formal +-{np.sqrt(cov[0, 0]):.1f}; "
                f"neighbouring drift values are correlated, so the real uncertainty is larger)")


def plots(out, adev, drift, idx, blk, t_wall, t0_wall, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ink, ink2, grid, surf, blue = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb", "#2a78d6"
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": ink2, "axes.labelcolor": ink, "xtick.color": ink2,
                         "ytick.color": ink2, "axes.spines.top": False, "axes.spines.right": False})

    pts = [a for a in adev if a["adev"]]
    tau = np.array([a["tau"] for a in pts])
    s = np.array([a["adev"] for a in pts])
    e = np.array([a["err"] for a in pts])
    fig, ax = plt.subplots(figsize=(6.4, 4.4), facecolor=surf)
    ax.set_facecolor(surf)
    ax.plot(tau, tau[0] * s[0] / tau, ls="--", lw=1, color=ink2)
    ax.text(tau[-1], tau[0] * s[0] / tau[-1], "  slope -1\n  (white phase noise)", color=ink2, fontsize=8, va="center")
    ax.errorbar(tau, s, yerr=e, color=blue, lw=2, marker="none", capsize=0, elinewidth=1)
    solid = np.array([a["n"] / a["tau"] >= MIN_INDEP for a in pts])
    ax.plot(tau[solid], s[solid], ls="none", marker="o", ms=7, color=blue, markeredgecolor=surf, markeredgewidth=1.5)
    ax.plot(tau[~solid], s[~solid], ls="none", marker="o", ms=7, color=surf, markeredgecolor=blue, markeredgewidth=1.5)
    if (~solid).any():
        ax.text(0.02, 0.03, f"hollow: fewer than {MIN_INDEP} independent intervals (too few to trust)",
                transform=ax.transAxes, color=ink2, fontsize=8)
    k = int(np.argmin(np.where(solid, s, np.inf)))
    ax.annotate(f"floor {s[k] * 1e9:.2g} ppb at {tau[k]} s", (tau[k], s[k]), textcoords="offset points",
                xytext=(0, -18), ha="center", color=ink, fontsize=9)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="major", color=grid, lw=0.8)
    ax.set_xlabel("averaging time tau (s)")
    ax.set_ylabel("overlapping Allan deviation")
    ax.set_title(f"{label}: sample clock against GPS", loc="left", color=ink, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "adev.png"), dpi=150, facecolor=surf)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 3.6), facecolor=surf)
    ax.set_facecolor(surf)
    ax.plot((t_wall[idx] - t0_wall) / 60.0, drift, color=blue, lw=1.2)
    if blk:
        ax.plot([b["t_min"] + BLOCK_S / 120.0 for b in blk], [b["f_ppb"] for b in blk], ls="none", marker="o", ms=6,
                color=ink, markeredgecolor=surf, markeredgewidth=1.5)
        ax.text(0.99, 0.02, "line: receiver drift (30 s fit)   dots: 5-min frequency from phase", transform=ax.transAxes,
                ha="right", va="bottom", color=ink2, fontsize=8)
    ax.grid(True, color=grid, lw=0.8)
    ax.set_xlabel("minutes since start")
    ax.set_ylabel("clock error (ppb)")
    ax.set_title(f"{label}: frequency error over the run", loc="left", color=ink, fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "drift.png"), dpi=150, facecolor=surf)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", required=True, help="directory for summary.json, adev.png, drift.png")
    ap.add_argument("--label", default="RSPdx")
    ap.add_argument("--warmup", action="store_true", help="report time to settle within 10 ppb of the final value")
    ap.add_argument("--t0", default="", help="ISO time to count minutes from (the plug-in), instead of the first row")
    ap.add_argument("--temps", default="", help="CSV of utc,temp_c (room temperature) to fit drift against")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    rows, t_wall, valid = load(a.csv)
    t0_wall = datetime.datetime.fromisoformat(a.t0).timestamp() if a.t0 else t_wall[0]
    segs = segments(rows, valid)
    if not segs:
        sys.exit("no valid fixes")
    longest = max(segs, key=lambda s: len(s[0]))
    f_all, f_se = phase_freq(longest[0], longest[1])
    adev = oadev(segs, TAUS)
    blk = blocks(segs, t0_wall, t_wall)
    drift, idx = settled_drift(segs)
    a300 = next((x["adev"] for x in adev if x["tau"] == 200 and x["adev"]), None)   # nearest listed tau to a block
    a300 = next((x["adev"] for x in adev if x["tau"] == 500 and x["adev"]), a300)
    best = min((x for x in adev if x["adev"] and x["n"] / x["tau"] >= MIN_INDEP), key=lambda x: x["adev"])
    bf = np.array([b["f_ppb"] for b in blk])

    S = dict(
        csv=os.path.basename(a.csv),
        log_start_utc=rows[0]["utc"], log_end_utc=rows[-1]["utc"], rows=len(rows),
        valid_fixes=int(valid.sum()), grid_holes=int(sum(np.isnan(sg[1]).sum() for sg in segs)), fix_availability=float(valid.mean()),
        segments=len(segs), resyncs=sum(r.get("resync") == "1" for r in rows),
        n_sats_median=float(np.median([int(r["n_sats"]) for r in rows if r["fix_valid"] == "1"])),
        n_sats_min=int(min(int(r["n_sats"]) for r in rows if r["fix_valid"] == "1")),
        n_sats_max=int(max(int(r["n_sats"]) for r in rows if r["fix_valid"] == "1")),
        residual_rms_m_median=float(np.median([float(r["residual_rms_m"]) for r in rows if r["fix_valid"] == "1"])),
        freq_from_phase_ppb=float(f_all), freq_from_phase_formal_se_ppb=float(f_se),
        freq_from_phase_segment_s=int(len(longest[0])),
        drift_mean_ppb=float(np.mean(drift)), drift_sd_ppb=float(np.std(drift, ddof=1)),
        drift_min_ppb=float(np.min(drift)), drift_max_ppb=float(np.max(drift)), drift_n=int(len(drift)),
        block_freqs=blk, block_sd_ppb=float(np.std(bf, ddof=1)) if len(bf) > 1 else None,
        adev=adev, adev_floor=dict(tau=best["tau"], adev_ppb=best["adev"] * 1e9),
        trend=trend(blk, a300))
    if a.warmup:
        S["warmup"] = warmup(drift, idx, t_wall, t0_wall)
    if a.temps:
        S["temperature"] = temp_fit(drift, idx, t_wall, a.temps)
    json.dump(S, open(os.path.join(a.out, "summary.json"), "w"), indent=1)
    plots(a.out, adev, drift, idx, blk, t_wall, t0_wall, a.label)

    print(f"{S['rows']} rows, {S['valid_fixes']} valid fixes ({100 * S['fix_availability']:.1f}%), "
          f"{S['segments']} segment(s), {S['grid_holes']} empty grid second(s), {S['resyncs']} resync(s); satellites {S['n_sats_min']}-{S['n_sats_max']} "
          f"(median {S['n_sats_median']:.0f}); residual rms median {S['residual_rms_m_median']:.2f} m")
    print(f"frequency from phase (longest segment, holes skipped, {S['freq_from_phase_segment_s']} s): {f_all:.2f} ppb "
          f"(formal +-{f_se:.3f})")
    print(f"receiver drift column (after its 30 s fit fills): mean {S['drift_mean_ppb']:.2f}, sd {S['drift_sd_ppb']:.2f}, "
          f"min {S['drift_min_ppb']:.2f}, max {S['drift_max_ppb']:.2f} ppb  (n={S['drift_n']})")
    if blk:
        print("5-min blocks (ppb): " + ", ".join(f"{b['f_ppb']:.1f}" for b in blk))
    print(f"trend: {S['trend']['verdict']}")
    print("\n| tau (s) | ADEV | ADEV (ppb) | ~1-sigma (ppb) | terms |\n|---:|---:|---:|---:|---:|")
    for x in adev:
        if x["adev"]:
            print(f"| {x['tau']} | {x['adev']:.2e} | {x['adev'] * 1e9:.3f} | {x['err'] * 1e9:.3f} | {x['n']} |")
        else:
            print(f"| {x['tau']} | - | - | - | 0 |")
    print(f"\nADEV floor: {S['adev_floor']['adev_ppb']:.3f} ppb at tau = {S['adev_floor']['tau']} s "
          f"(among taus with >= {MIN_INDEP} independent intervals)")
    for k in ("warmup", "temperature"):
        if k in S:
            print(f"{k}: {S[k]['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
