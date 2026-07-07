#!/usr/bin/env python3
"""
convergence_analysis.py  --  per-leg ABFE sampling-convergence diagnostic
=========================================================================

Answers two questions for an Amber TI/MBAR leg, using the REAL prod.out files:

  Q1. Has the leg's dG PLATEAUED by the run length?  (cumulative-time scan)
  Q2. Do the FORWARD and REVERSE halves agree?       (time-reversal check)

If both hold, the current per-window simulation time is sufficient; if dG is
still drifting at the final time, or forward/reverse disagree beyond the
combined error bar, the leg is under-sampled and needs MORE ns (not fewer).

This is read-only. It never writes into the leg dirs except the output CSV/PNG
you point --out at (default: <leg_dir>/convergence.{csv,png,log}).

Usage
-----
    python convergence_analysis.py <leg_dir> [--temp 298.0] \
        [--fractions 0.2 0.4 0.6 0.8 1.0] [--estimator auto|ti|bar|mbar] \
        [--out <prefix>]

<leg_dir> contains lambda_<x.xxx>/prod.out for every window.

Method
------
* Parse each window once with alchemlyb.parsing.amber.extract_u_nk / extract_dHdl
  (the sanitized production path). dHdl drives TI; u_nk drives BAR/MBAR.
* CUMULATIVE scan: for each fraction f in --fractions, keep the FIRST f of every
  window's samples and recompute dG. A converged leg flattens as f -> 1.0.
* FORWARD/REVERSE: dG from the first 50% vs the last 50% of every window.
  |dG_fwd - dG_rev| should be within sqrt(err_fwd^2 + err_rev^2) (we flag > 1 kT
  ~ 0.6 kcal/mol OR > combined 1-sigma, whichever is larger).
* Estimator selection: 'auto' uses MBAR if it solves, else BAR, else TI. TI is
  always reported because it survives the high-lambda charging singularity that
  makes BAR/MBAR fail on the (pre-v2.5.70) decharge leg.

Everything is logged to BOTH stdout and <out>.log so it can be uploaded.
"""
# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.
import os, re, sys, glob, math, argparse, datetime

KB = 0.0019872041  # kcal/mol/K


class Tee:
    def __init__(self, path):
        self.f = open(path, "w"); self.s = sys.stdout
    def write(self, x): self.s.write(x); self.f.write(x)
    def flush(self): self.s.flush(); self.f.flush()


def _ti_dG(dvdl_by_lam, lams):
    """Trapezoid TI over mean dV/dl, extrapolating absent 0/1 endpoints."""
    import numpy as np
    xs = sorted(dvdl_by_lam)
    ys = [dvdl_by_lam[x] for x in xs]
    if xs and xs[0] > 0.0:
        y0 = ys[0] + (ys[0]-ys[1])*(0.0-xs[0])/(xs[0]-xs[1]) if len(xs) > 1 else ys[0]
        xs = [0.0]+xs; ys = [y0]+ys
    if xs and xs[-1] < 1.0:
        y1 = ys[-1] + (ys[-1]-ys[-2])*(1.0-xs[-1])/(xs[-1]-xs[-2]) if len(xs) > 1 else ys[-1]
        xs = xs+[1.0]; ys = ys+[y1]
    _trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy>=2 renamed trapz
    return float(_trapz(ys, xs))


def _bar_dG(u_nk, kT):
    """Adjacent-pair BAR over an alchemlyb u_nk frame. Robust to empty/degenerate."""
    import numpy as np
    from pymbar.other_estimators import bar as bar_fn
    cols = list(u_nk.columns)
    src = u_nk.index.get_level_values(-1).to_numpy()
    uniq = sorted(set(src.tolist()))
    def near(v): return min(uniq, key=lambda s: abs(float(s)-float(v)))
    dG = 0.0; var = 0.0; npairs = 0
    for i in range(len(cols)-1):
        ci, cj = cols[i], cols[i+1]
        li = ci[-1] if isinstance(ci, tuple) else ci
        lj = cj[-1] if isinstance(cj, tuple) else cj
        mi = src == near(li); mj = src == near(lj)
        wF = (u_nk[cj].to_numpy()[mi] - u_nk[ci].to_numpy()[mi])
        wR = (u_nk[ci].to_numpy()[mj] - u_nk[cj].to_numpy()[mj])
        wF = wF[np.isfinite(wF)]; wR = wR[np.isfinite(wR)]
        if len(wF) == 0 or len(wR) == 0:
            return None, None  # cannot do BAR on this slice
        try:
            o = bar_fn(wF, wR, relative_tolerance=1e-4, maximum_iterations=10000)
            df, de = o['Delta_f'], o['dDelta_f']
        except Exception:
            try:
                o = bar_fn(wF, wR, relative_tolerance=1e-4, maximum_iterations=10000,
                           compute_uncertainty=False)
                df, de = o['Delta_f'], 0.0
            except Exception:
                return None, None
        dG += df*kT; var += ((de if np.isfinite(de) else 0.0)*kT)**2; npairs += 1
    return dG, math.sqrt(var)


def _mbar_dG(u_nk, kT):
    try:
        from alchemlyb.estimators import MBAR
        m = MBAR(); m.fit(u_nk)
        return float(m.delta_f_.iloc[0, -1])*kT, float(m.d_delta_f_.iloc[0, -1])*kT
    except Exception:
        return None, None


def main():
    import numpy as np
    ap = argparse.ArgumentParser()
    ap.add_argument("legdir")
    ap.add_argument("--temp", type=float, default=298.0)
    ap.add_argument("--fractions", type=float, nargs="+",
                    default=[0.2, 0.4, 0.6, 0.8, 1.0])
    ap.add_argument("--estimator", choices=["auto", "ti", "bar", "mbar"],
                    default="auto")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    pref = args.out or os.path.join(args.legdir, "convergence")
    sys.stdout = Tee(pref + ".log")
    kT = KB * args.temp

    print("="*72)
    print("CONVERGENCE ANALYSIS")
    print("="*72)
    print("when    :", datetime.datetime.now().isoformat())
    print("leg_dir :", os.path.abspath(args.legdir))
    print("temp    :", args.temp, " kT =", round(kT, 5), "kcal/mol")
    print("frac    :", args.fractions, " estimator:", args.estimator)

    try:
        import pandas as pd
        from alchemlyb.parsing.amber import extract_u_nk, extract_dHdl
    except Exception as e:
        print("FATAL: need pandas + alchemlyb:", e); sys.exit(1)

    lam_dirs = sorted(glob.glob(os.path.join(args.legdir, "lambda_*")))
    lams = sorted(float(re.search(r'lambda_([\d.]+)', d).group(1)) for d in lam_dirs)
    print(f"\n{len(lams)} windows: " + ", ".join(f"{l:.3f}" for l in lams))

    # Parse once; store full per-window frames.
    u_full, h_full, found = {}, {}, []
    for l in lams:
        p = os.path.join(args.legdir, f"lambda_{l:.3f}", "prod.out")
        if not os.path.exists(p):
            print(f"  [skip] missing {p}"); continue
        try:
            _u = extract_u_nk(p, T=args.temp)
        except Exception as e:
            print(f"  [warn] u_nk parse failed lambda={l:.3f}: {e}"); _u = None
        # v2.5.73 BUGFIX: a None/empty u_nk (e.g. fully-decoupled endpoint
        # lambda=1.000, or a half-decoupled window alchemlyb returns no u_nk for)
        # must NOT be added to `found`, else len(u_full[l]) crashes (line 166).
        if _u is None or len(_u) == 0:
            print(f"  [skip] no u_nk at lambda={l:.3f} (endpoint); dHdl-only for TI")
            try:
                h = extract_dHdl(p, T=args.temp)
                if h is not None and len(h) > 0:
                    h_full[l] = h
            except Exception:
                pass
            continue
        u_full[l] = _u
        # dHdl may legitimately be absent at decharge endpoints (lambda 0/1)
        try:
            h = extract_dHdl(p, T=args.temp)
            if h is not None and len(h) > 0:
                h_full[l] = h
        except Exception as e:
            print(f"  [warn] dHdl absent lambda={l:.3f}: {e}")
        found.append(l)
    if not found:
        print("No parseable windows."); sys.exit(1)
    found_h = sorted(h_full)          # windows usable for TI
    # frame count comes from u_nk (always present); used for time slicing
    nframes = {l: len(u_full[l]) for l in found}
    print("u_nk windows :", [f"{l:.3f}" for l in found])
    print("dHdl windows :", [f"{l:.3f}" for l in found_h],
          "  (endpoints without dV/dl are TI-extrapolated)")
    print("frames/window:", {f"{l:.3f}": nframes[l] for l in found})

    # ---- cumulative-time scan ----
    print("\n" + "-"*72)
    print("Q1  CUMULATIVE-TIME SCAN  (dG using first f of every window)")
    print("-"*72)
    print(f"{'frac':>6} {'ns~':>6} {'TI':>10} {'BAR':>14} {'MBAR':>14}")
    rows = []
    # estimate ns from frames assuming the leg's nominal length is the last fraction
    for f in args.fractions:
        dvdl = {}
        u_parts = []
        for l in found:
            n = max(2, int(round(nframes[l]*f)))
            if l in h_full:
                nh = max(2, int(round(len(h_full[l])*f)))
                dvdl[l] = float(h_full[l].iloc[:nh].to_numpy().mean())
            u_parts.append(u_full[l].iloc[:n])
        u_cat = pd.concat(u_parts)
        ti = _ti_dG(dvdl, found_h)
        bar, bar_e = _bar_dG(u_cat, kT)
        mb, mb_e = _mbar_dG(u_cat, kT)
        rows.append((f, ti, bar, bar_e, mb, mb_e))
        bs = f"{bar:8.3f}+/-{bar_e:.2f}" if bar is not None else "      n/a     "
        ms = f"{mb:8.3f}+/-{mb_e:.2f}" if mb is not None else "      n/a     "
        print(f"{f:>6.2f} {f*5:>6.2f} {ti:>10.3f} {bs:>14} {ms:>14}")

    # ---- forward / reverse ----
    print("\n" + "-"*72)
    print("Q2  FORWARD vs REVERSE HALVES  (time-reversal consistency)")
    print("-"*72)
    def half(which):
        dvdl = {}; parts = []
        for l in found:
            n = nframes[l]; u = u_full[l]
            sl = slice(0, n//2) if which == "fwd" else slice(n//2, n)
            if l in h_full:
                h = h_full[l]; nh = len(h)
                hsl = slice(0, nh//2) if which == "fwd" else slice(nh//2, nh)
                dvdl[l] = float(h.iloc[hsl].to_numpy().mean())
            parts.append(u.iloc[sl])
        uc = pd.concat(parts)
        return _ti_dG(dvdl, found_h), _bar_dG(uc, kT), _mbar_dG(uc, kT)
    ti_f, (bar_f, bef), (mb_f, mef) = half("fwd")
    ti_r, (bar_r, ber), (mb_r, mer) = half("rev")
    print(f"  TI    fwd={ti_f:8.3f}  rev={ti_r:8.3f}  |diff|={abs(ti_f-ti_r):.3f}")
    if bar_f is not None and bar_r is not None:
        print(f"  BAR   fwd={bar_f:8.3f}  rev={bar_r:8.3f}  |diff|={abs(bar_f-bar_r):.3f}")
    if mb_f is not None and mb_r is not None:
        print(f"  MBAR  fwd={mb_f:8.3f}  rev={mb_r:8.3f}  |diff|={abs(mb_f-mb_r):.3f}")

    # ---- verdict ----
    # choose headline estimator
    last = rows[-1]
    cand = {"ti": last[1], "bar": last[2], "mbar": last[4]}
    if args.estimator == "auto":
        head = "mbar" if cand["mbar"] is not None else ("bar" if cand["bar"] is not None else "ti")
    else:
        head = args.estimator
    # plateau: |dG(1.0) - dG(0.8)| small relative to spread
    plateau_drift = None
    if len(rows) >= 2:
        col = {"ti": 1, "bar": 2, "mbar": 4}[head]
        a, b = rows[-2][col], rows[-1][col]
        if a is not None and b is not None:
            plateau_drift = abs(b - a)
    fr_diff = {"ti": abs(ti_f-ti_r),
               "bar": (abs(bar_f-bar_r) if bar_f is not None and bar_r is not None else None),
               "mbar": (abs(mb_f-mb_r) if mb_f is not None and mb_r is not None else None)}[head]

    print("\n" + "="*72); print("VERDICT"); print("="*72)
    print(f"headline estimator: {head.upper()}  (auto = MBAR>BAR>TI by solvability)")
    GATE = 0.6  # kcal/mol ~ 1 kT
    def tag(v): 
        return "n/a" if v is None else (f"{v:.3f}  ({'OK' if v <= GATE else 'HIGH'})")
    print(f"last-step plateau drift |dG(f=last)-dG(f=prev)| = {tag(plateau_drift)}")
    print(f"forward/reverse |dG_fwd - dG_rev|              = {tag(fr_diff)}")
    converged = (plateau_drift is not None and plateau_drift <= GATE and
                 fr_diff is not None and fr_diff <= GATE)
    if converged:
        print("\n=> CONVERGED at the current run length. Cutting time may be safe;\n"
              "   the current per-window ns is sufficient for this leg.")
    else:
        print("\n=> NOT demonstrably converged. The leg is still drifting and/or\n"
              "   forward!=reverse. Do NOT cut sampling time; consider MORE ns\n"
              "   (or denser lambda / replica exchange) for this leg.")
    if head == "ti" and cand["mbar"] is None:
        print("   NOTE: MBAR/BAR unsolvable here (overlap collapse) -- this leg has\n"
              "   a protocol problem; convergence of TI alone is not sufficient.")

    # ---- CSV + plot ----
    import csv
    with open(pref + ".csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["fraction", "approx_ns", "TI", "BAR", "BAR_err", "MBAR", "MBAR_err"])
        for (f, ti, bar, be, mb, me) in rows:
            w.writerow([f, f*5, f"{ti:.4f}",
                        "" if bar is None else f"{bar:.4f}",
                        "" if be is None else f"{be:.4f}",
                        "" if mb is None else f"{mb:.4f}",
                        "" if me is None else f"{me:.4f}"])
        w.writerow([])
        w.writerow(["half", "TI", "BAR", "MBAR"])
        w.writerow(["forward", f"{ti_f:.4f}",
                    "" if bar_f is None else f"{bar_f:.4f}",
                    "" if mb_f is None else f"{mb_f:.4f}"])
        w.writerow(["reverse", f"{ti_r:.4f}",
                    "" if bar_r is None else f"{bar_r:.4f}",
                    "" if mb_r is None else f"{mb_r:.4f}"])
    print(f"\n[csv] {pref}.csv")

    try:
        import matplotlib; matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fr = [r[0] for r in rows]
        plt.figure(figsize=(7, 4.5))
        plt.plot(fr, [r[1] for r in rows], "o-", label="TI")
        if all(r[2] is not None for r in rows):
            plt.plot(fr, [r[2] for r in rows], "s-", label="BAR")
        if all(r[4] is not None for r in rows):
            plt.plot(fr, [r[4] for r in rows], "^-", label="MBAR")
        plt.xlabel("fraction of trajectory used (1.0 = full run)")
        plt.ylabel("cumulative dG (kcal/mol)")
        plt.title(f"{os.path.basename(os.path.abspath(args.legdir))} convergence")
        plt.grid(alpha=0.3); plt.legend(); plt.tight_layout()
        plt.savefig(pref + ".png", dpi=140); plt.close()
        print(f"[png] {pref}.png")
    except Exception as e:
        print("[plot skipped]", e)

    print("\nLog:", pref + ".log")


if __name__ == "__main__":
    main()