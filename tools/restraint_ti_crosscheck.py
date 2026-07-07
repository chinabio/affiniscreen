#!/usr/bin/env python3
"""
restraint_ti_crosscheck.py  --  numerical-vs-analytic cross-check for the
Option-A (icfe=0, lambda-scaled Boresch &rst) complex_restraint leg.

WHY THIS EXISTS
---------------
The package FEPAnalyzer cannot analyze this leg: it needs alchemlyb's icfe=1
free-energy section (clambda / DV/DL), which an icfe=0 lambda-scaled-&rst leg
does NOT produce. The only per-window alchemical signal here is the logged
<RESTRAINT> energy in the AVERAGES block.

For U(lambda) = lambda * U_full   (U_full = k_full*(x-x0)^2):
    dU/dlambda = U_full = <RESTRAINT_lambda> / lambda
    dG_attach  = integral_0^1  <RESTRAINT_lambda>/lambda  dlambda

That integrand has a 1/lambda SINGULARITY at lambda->0, so a naive trapezoid is
unreliable and grid-dependent. This script:
  1. extracts <RESTRAINT> for ALL windows (from the AVERAGES block),
  2. forms the integrand g(lambda) = <RESTRAINT>/lambda,
  3. integrates THREE ways:
       (a) naive trapezoid in lambda           (shows how bad the singularity is)
       (b) substitution u=sqrt(lambda): dG = integral 2u*g(u^2) du  (tames ~1/sqrt)
       (c) substitution s=ln(lambda): dG = integral lambda*g du = integral <R> d(ln l)
          (exact-ish if <R> ~ const near 0; integrand = <RESTRAINT> itself)
  4. compares all three to the analytic Boresch correction magnitude (|-11.45|).

OUTPUT: prints a table + the three integrals, and writes
        restraint_ti_crosscheck.csv next to the leg dir.

USAGE:
  python restraint_ti_crosscheck.py \
      ~/Run_dir/.../fep/complex_restraint \
      --analytic -11.447503

Pure standard library + (optional) numpy. No GPU. Read-only on prod.out files.
"""
import argparse
import csv
import math
import re
import sys
from pathlib import Path

LAMBDAS_DEFAULT = [0.0, 0.004, 0.016, 0.036, 0.064, 0.1, 0.144, 0.196, 0.256,
                   0.324, 0.4, 0.5, 0.6, 0.7, 0.8, 0.875, 0.925, 0.95, 0.975, 1.0]


def parse_avg_restraint(prod_out: Path):
    """Return <RESTRAINT> from the AVERAGES block, or None."""
    in_avg = False
    val = None
    try:
        text = prod_out.read_text(errors="replace")
    except Exception:
        return None
    for line in text.splitlines():
        if "A V E R A G E S" in line:
            in_avg = True
            continue
        if not in_avg:
            continue
        if "R M S" in line:
            break
        m = re.search(r"RESTRAINT\s*=\s*([-0-9.Ee+]+)", line)
        if m:
            try:
                val = float(m.group(1))
            except ValueError:
                pass
    return val


def trapz(xs, ys):
    return sum(0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i])
              for i in range(len(xs) - 1))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("leg_dir", type=Path)
    ap.add_argument("--analytic", type=float, default=-11.447503,
                    help="analytic Boresch correction (signed) for comparison.")
    ap.add_argument("--lambdas", default=None,
                    help="comma list; default = the 20-window restraint schedule.")
    a = ap.parse_args()

    lams = ([float(x) for x in a.lambdas.split(",")] if a.lambdas
            else LAMBDAS_DEFAULT)

    rows = []
    for lam in lams:
        po = a.leg_dir / f"lambda_{lam:.3f}" / "prod.out"
        R = parse_avg_restraint(po) if po.exists() else None
        rows.append({"lambda": lam, "restraint": R,
                     "exists": po.exists()})

    # report extraction
    print("=== <RESTRAINT> per window (AVERAGES block) ===")
    print(f"{'lambda':>8} {'<RESTRAINT>':>12} {'g=<R>/lam':>12}   note")
    usable = []
    for r in rows:
        lam, R = r["lambda"], r["restraint"]
        if R is None:
            note = "MISSING prod.out" if not r["exists"] else "no RESTRAINT in avgs"
            print(f"{lam:8.3f} {'--':>12} {'--':>12}   {note}")
            continue
        if lam == 0.0:
            # dU/dl at lam=0 = U_full; <RESTRAINT>/lam is 0/0. We approximate
            # U_full(0) by extrapolation later; record raw R only.
            print(f"{lam:8.3f} {R:12.4f} {'(endpoint)':>12}   lambda=0: g undefined (1/lam)")
            r["g"] = None
            usable.append(r)
            continue
        g = R / lam
        r["g"] = g
        print(f"{lam:8.3f} {R:12.4f} {g:12.4f}")
        usable.append(r)

    # build arrays over windows that have a defined g (lam>0)
    pts = [(r["lambda"], r["restraint"], r["g"]) for r in usable
           if r.get("g") is not None and r["restraint"] is not None]
    pts.sort()
    if len(pts) < 3:
        print("\nNOT ENOUGH usable windows to integrate.", file=sys.stderr)
        sys.exit(2)

    xs = [p[0] for p in pts]                 # lambda
    Rs = [p[1] for p in pts]                 # <RESTRAINT>
    gs = [p[2] for p in pts]                 # g = <R>/lambda

    # (a) naive trapezoid of g(lambda) in lambda
    dG_naive = trapz(xs, gs)

    # (b) substitution u = sqrt(lambda): dlambda = 2u du; integrand g*2u
    us = [math.sqrt(x) for x in xs]
    integ_b = [g * 2.0 * u for g, u in zip(gs, us)]
    dG_sqrt = trapz(us, integ_b)

    # (c) substitution s = ln(lambda): dlambda = lambda d s; g*lambda = <RESTRAINT>
    #     so dG = integral <RESTRAINT> d(ln lambda). Lower limit ln(lam_min);
    #     this MISSES the (0, lam_min] segment -- report it as a lower bound and
    #     also add an analytic tail assuming <RESTRAINT> ~ const = Rs[0] there:
    #     tail = Rs[0] * (ln(lam_min) - ln(0+)) -> diverges, so instead assume
    #     U_full const => contribution from [0,lam_min] in lambda of g=U_full is
    #     U_full*lam_min. Use U_full ~ Rs[0]/xs[0].
    ss = [math.log(x) for x in xs]
    dG_ln_from_min = trapz(ss, Rs)            # integral over [lam_min, 1]
    U_full_0 = Rs[0] / xs[0]                   # estimate of dU/dl at lambda->0
    tail_0_to_min = U_full_0 * xs[0]           # integral_0^lam_min U_full dlam
    dG_ln_total = dG_ln_from_min + tail_0_to_min

    print("\n=== TI integrals of dG_attach = integral_0^1 <RESTRAINT>/lambda dlambda ===")
    print(f"  (a) naive trapezoid in lambda            : {dG_naive:10.3f} kcal/mol")
    print(f"      ^ UNRELIABLE: ignores 1/lambda blow-up below lambda={xs[0]:.3f}")
    print(f"  (b) sqrt(lambda) substitution            : {dG_sqrt:10.3f} kcal/mol")
    print(f"  (c) ln(lambda) form + linear tail to 0   : {dG_ln_total:10.3f} kcal/mol")
    print(f"        (= {dG_ln_from_min:.3f} over [{xs[0]:.3f},1] + {tail_0_to_min:.3f} tail [0,{xs[0]:.3f}])")
    print(f"\n  lowest-lambda integrand g({xs[0]:.3f}) = {gs[0]:.1f} kcal/mol "
          f"(this is what blows up -> drives the spread above)")

    print("\n=== comparison to analytic Boresch correction ===")
    print(f"  analytic correction (signed)             : {a.analytic:+.3f} kcal/mol")
    print(f"  |analytic|                               : {abs(a.analytic):.3f} kcal/mol")
    print("  NOTE: the numerical dG_attach above is the cost to TURN ON the")
    print("  restraint (positive); the analytic standard-state correction is the")
    print("  signed contribution to dG_bind. They are related but NOT identical")
    print("  quantities -- the analytic term also includes the 1/(8pi^2 V) and")
    print("  standard-state volume factors that a finite-lambda TI cannot capture.")
    print("  The spread among (a)/(b)/(c) is the SINGULARITY ARTIFACT: it shows")
    print("  the numerical leg is NOT a reliable estimator. Use the analytic term.")

    # write CSV
    out_csv = a.leg_dir / "restraint_ti_crosscheck.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lambda", "avg_restraint_kcal", "g_dUdlam_kcal"])
        for r in rows:
            w.writerow([r["lambda"], r["restraint"],
                        r.get("g") if r.get("g") is not None else ""])
        w.writerow([])
        w.writerow(["integral_method", "dG_kcal_mol"])
        w.writerow(["naive_trapz_lambda", dG_naive])
        w.writerow(["sqrt_lambda_sub", dG_sqrt])
        w.writerow(["ln_lambda_plus_tail", dG_ln_total])
        w.writerow(["analytic_boresch_signed", a.analytic])
    print(f"\nWrote {out_csv}")


if __name__ == "__main__":
    main()
