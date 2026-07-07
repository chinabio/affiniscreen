#!/usr/bin/env python
"""Generate a styled HTML report from a FINAL_RESULTS_MMPBSA.dat file.

Usage (CLI):
    python mmpbsa_report.py path/to/FINAL_RESULTS_MMPBSA.dat
    python mmpbsa_report.py path/to/workdir          # auto-finds mmgbsa/FINAL_RESULTS_MMPBSA.dat
    python mmpbsa_report.py path/to/batch_dir        # generates one report per lig_*

Usage (programmatic):
    from mmpbsa_report import generate_report
    html_path = generate_report("/path/to/FINAL_RESULTS_MMPBSA.dat")
    # writes <same dir>/FINAL_RESULTS.report.html and returns its Path
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import sys, re, html, argparse, base64, io, datetime
from pathlib import Path

# Module-level unicode constants. These are referenced from inside f-string
# {expressions}, where on Python <=3.11 a literal backslash escape like
# "\u2014" is a SyntaxError ("f-string expression part cannot include a
# backslash"). Hoisting them out of the f-string expression sidesteps that
# rule and keeps this file source-compatible with 3.8 - 3.13+.
EMDASH  = "\u2014"   # em-dash
DELTA   = "\u0394"   # Greek capital delta
PLUSMIN = "\u00b1"   # plus-minus sign

# ----------------------------- parsing -----------------------------

ENERGY_TERMS = ["BOND", "ANGLE", "DIHED", "VDWAALS", "EEL",
                "1-4 VDW", "1-4 EEL", "EGB", "ESURF"]
SCALAR_TERMS = ["G gas", "G solv", "TOTAL",
                "DELTA G gas", "DELTA G solv", "DELTA TOTAL"]

# Section headers we will look for. Order matters when splitting the file.
_SECTIONS = [
    ("Complex:",  "complex"),
    ("Receptor:", "receptor"),
    ("Ligand:",   "ligand"),
    ("Differences (Complex - Receptor - Ligand):", "delta"),
]

def _extract_section(text, header):
    """Return the body lines belonging to a section, or None.

    Body = all lines from the section header up to the next section header
    in `_SECTIONS` (or end of file). Robust to blank lines inside the body.
    """
    idx = text.find(header)
    if idx < 0:
        return None
    start = idx + len(header)
    # Find the next section header *after* this one.
    other_headers = [h for h, _ in _SECTIONS if h != header]
    cut = len(text)
    for h in other_headers:
        j = text.find(h, start)
        if 0 <= j < cut:
            cut = j
    return text[start:cut]

def _parse_terms_from_body(body):
    """Pull energy terms out of one section body."""
    terms = {}
    for line in body.splitlines():
        line = line.rstrip()
        if not line.strip():
            continue
        # Try the longest term names first so 'DELTA TOTAL' wins over 'TOTAL'.
        for t in sorted(ENERGY_TERMS + SCALAR_TERMS, key=len, reverse=True):
            if line.startswith(t + " ") or line.startswith(t + "\t"):
                nums = re.findall(r"-?\d+\.\d+", line[len(t):])
                if len(nums) >= 3:
                    terms[t] = {
                        "avg": float(nums[0]),
                        "std": float(nums[1]),
                        "sem": float(nums[2]),
                    }
                break
    return terms

def parse_mmpbsa_dat(path):
    """Extract complex/receptor/ligand/delta blocks and metadata."""
    path = Path(path)
    text = path.read_text(errors="replace")

    meta = {
        "path": str(path),
        "run_date": None,
        "complex_prmtop": None,
        "trajectory": None,
        "frames": None,
        "ligand_resname": None,
        "igb": None,
        "saltcon": None,
    }
    m = re.search(r"\|\s*Run on\s+(.+)", text)
    if m: meta["run_date"] = m.group(1).strip()
    m = re.search(r"\|Solvated complex topology file:\s*(\S+)", text)
    if m: meta["complex_prmtop"] = m.group(1).strip()
    m = re.search(r"\|Initial mdcrd.*?:\s*(\S+)", text)
    if m: meta["trajectory"] = m.group(1).strip()
    m = re.search(r"using\s+([\d.]+)\s+complex frames", text)
    if m: meta["frames"] = int(float(m.group(1)))
    m = re.search(r'Ligand residue name is\s+"([^"]+)"', text)
    if m: meta["ligand_resname"] = m.group(1)
    m = re.search(r"igb\s*=\s*(\d+)", text)
    if m: meta["igb"] = int(m.group(1))
    m = re.search(r"saltcon\s*=\s*([\d.]+)", text)
    if m: meta["saltcon"] = float(m.group(1))

    blocks = {}
    for header, key in _SECTIONS:
        body = _extract_section(text, header)
        if body is not None:
            blocks[key] = _parse_terms_from_body(body)

    return {"meta": meta, "blocks": blocks}


# ----------------------------- interpretation -----------------------------

INTERPRET = {
    "VDWAALS": {
        "name": "Van der Waals",
        "what": "Shape / packing complementarity",
        "good": "Dominant driver. Tight steric fit in the pocket; many close contacts. "
                "Suggests the binding mode is geometrically well-defined.",
        "neutral": "Modest VDW contribution; some shape fit but not the main driver.",
        "weak": "Weak VDW - poor pocket fit; consider re-examining the docking pose.",
    },
    "EEL": {
        "name": "Coulomb electrostatic",
        "what": "Direct charge-charge / dipole interactions in the complex",
        "good": "Strong favorable electrostatics - likely salt bridges or directed H-bonds.",
        "neutral": "Modest favorable electrostatics - likely some H-bonds or charge-pair contacts. Not the main driver.",
        "weak": "Near-zero or unfavorable Coulomb - little polar interaction in the binding mode.",
    },
    "EGB": {
        "name": "GB polar solvation",
        "what": "Cost of removing waters from the pocket and partially desolvating the ligand",
        "good": "",
        "neutral": "Standard desolvation penalty for burying polar groups on binding.",
        "weak": "",
    },
    "ESURF": {
        "name": "Non-polar solvation (cavity + dispersion)",
        "what": "Hydrophobic contribution from burying surface area",
        "good": "Favorable hydrophobic burial - adds to the binding affinity beyond VDW.",
        "neutral": "Small favorable hydrophobic contribution.",
        "weak": "Weak / unfavorable - little hydrophobic surface buried.",
    },
}

def classify_term(term, val):
    if term == "VDWAALS":
        if val < -30: return "good"
        if val < -10: return "neutral"
        return "weak"
    if term == "EEL":
        if val < -30: return "good"
        if val < -5:  return "neutral"
        return "weak"
    if term == "ESURF":
        if val < -5: return "good"
        if val < -2: return "neutral"
        return "weak"
    return "neutral"

def overall_verdict(dG, sem):
    if dG < -40: return ("Strong predicted binder",   "#2b8a3e")
    if dG < -25: return ("Moderate predicted binder", "#1971c2")
    if dG < -10: return ("Weak predicted binder",     "#f08c00")
    return       ("Unfavorable / non-binder",         "#c92a2a")


# ----------------------------- plots -----------------------------

def make_waterfall_png_b64(delta):
    """Return a base64-encoded PNG of the energy decomposition bar chart."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    terms = [t for t in ("VDWAALS", "EEL", "EGB", "ESURF") if t in delta]
    vals = [delta[t]["avg"] for t in terms]
    errs = [delta[t]["sem"] for t in terms]
    total = delta.get("DELTA TOTAL", {}).get("avg")
    total_err = delta.get("DELTA TOTAL", {}).get("sem")

    if total is not None:
        terms.append("DELTA\nTOTAL"); vals.append(total); errs.append(total_err or 0.0)

    colors = []
    for t, v in zip(terms, vals):
        if t.startswith("DELTA"):
            colors.append("#1971c2")
        elif v < 0:
            colors.append("#2b8a3e")
        else:
            colors.append("#c92a2a")

    fig, ax = plt.subplots(figsize=(7.5, 3.5))
    y = list(range(len(terms)))[::-1]
    bars = ax.barh(y, vals, xerr=errs, color=colors, ecolor="#444",
                   capsize=4, edgecolor="black", linewidth=0.6)
    ax.set_yticks(y); ax.set_yticklabels(terms, fontsize=10)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("Contribution to \u0394G_bind  (kcal/mol)", fontsize=10)
    ax.set_title("MM/GBSA energy decomposition", fontsize=11, weight="bold")
    ax.grid(axis="x", alpha=0.25, linestyle="--")
    for bar, val in zip(bars, vals):
        x_text = val + (1.2 if val >= 0 else -1.2)
        ha = "left" if val >= 0 else "right"
        ax.text(x_text, bar.get_y() + bar.get_height()/2,
                f"{val:+.1f}", va="center", ha=ha, fontsize=9, color="#222")
    ax.set_xlim(min(vals)*1.25 - 5, max(vals)*1.25 + 5)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


# ----------------------------- HTML -----------------------------

CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       max-width: 980px; margin: 24px auto; padding: 0 20px; color: #222; line-height: 1.5; }
h1 { border-bottom: 3px solid #1971c2; padding-bottom: 6px; }
h2 { color: #1971c2; margin-top: 28px; }
.hero { padding: 18px 22px; border-radius: 8px; margin: 16px 0; color: white;
        font-size: 1.15em; }
.hero .big { font-size: 1.9em; font-weight: 700; display: block; margin-top: 4px; }
table { border-collapse: collapse; width: 100%; margin: 10px 0 22px; font-size: 0.94em; }
th, td { padding: 7px 10px; border-bottom: 1px solid #e0e0e0; text-align: left;
         vertical-align: top; }
th { background: #f5f7fa; font-weight: 600; }
tr:nth-child(even) td { background: #fafbfc; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
.good    { color: #2b8a3e; font-weight: 600; }
.bad     { color: #c92a2a; font-weight: 600; }
.neutral { color: #495057; font-weight: 600; }
.meta { font-size: 0.85em; color: #666; }
.caveat { background: #fff3bf; border-left: 4px solid #f59f00;
          padding: 10px 14px; border-radius: 4px; margin: 14px 0; font-size: 0.93em; }
.footer { font-size: 0.8em; color: #888; border-top: 1px solid #e0e0e0;
          padding-top: 10px; margin-top: 30px; }
img.chart { max-width: 100%; border: 1px solid #e0e0e0; border-radius: 6px;
            padding: 8px; background: white; }
"""

def fmt_kcal(v, plus=True):
    if v is None: return EMDASH
    return f"{v:+.2f}" if plus else f"{v:.2f}"

def build_html(parsed, src_path):
    meta = parsed["meta"]
    delta = parsed["blocks"].get("delta", {})

    dT = delta.get("DELTA TOTAL")
    if not dT:
        return (f"<html><body><h1>Cannot parse</h1>"
                f"<p>No DELTA TOTAL block found in {html.escape(str(src_path))}.</p>"
                f"</body></html>")

    dG  = dT["avg"]; std = dT["std"]; sem = dT["sem"]
    dGgas  = delta.get("DELTA G gas",  {}).get("avg")
    dGsolv = delta.get("DELTA G solv", {}).get("avg")
    verdict, vcolor = overall_verdict(dG, sem)

    term_rows_html = []
    for term in ("VDWAALS", "EEL", "EGB", "ESURF"):
        d = delta.get(term)
        if not d: continue
        v = d["avg"]
        info = INTERPRET[term]
        cls = "good" if v < 0 else "bad"
        bucket = "neutral" if term == "EGB" else classify_term(term, v)
        read = info.get(bucket) or info.get("neutral") or ""
        term_rows_html.append(
            f"<tr>"
            f"<td><b>{term}</b><br><span class='meta'>{html.escape(info['name'])}</span></td>"
            f"<td class='num {cls}'>{fmt_kcal(v)}</td>"
            f"<td>{html.escape(info['what'])}</td>"
            f"<td>{html.escape(read)}</td>"
            f"</tr>")

    # Pre-compute conditional classes outside f-string (Python <=3.11 limitation).
    gas_cls  = "good" if (dGgas  or 0) < 0 else "bad"
    solv_cls = "good" if (dGsolv or 0) < 0 else "bad"
    bind_cls = "good" if dG < 0 else "bad"
    summary_rows_html = [
        f"<tr><td><b>{DELTA}G_gas</b><br><span class='meta'>Gas-phase interaction</span></td>"
        f"<td class='num {gas_cls}'>{fmt_kcal(dGgas)}</td>"
        f"<td colspan='2'>Sum of VDW + EEL + bonded terms (in single-trajectory MM/GBSA, "
        f"bonded terms cancel exactly).</td></tr>",
        f"<tr><td><b>{DELTA}G_solv</b><br><span class='meta'>Net solvation cost</span></td>"
        f"<td class='num {solv_cls}'>{fmt_kcal(dGsolv)}</td>"
        f"<td colspan='2'>EGB + ESURF: net cost of taking the ligand and the binding pocket "
        f"out of bulk water.</td></tr>",
        f"<tr><td><b>{DELTA}G_bind (TOTAL)</b><br><span class='meta'>Predicted binding free energy</span></td>"
        f"<td class='num {bind_cls}'>{fmt_kcal(dG)}</td>"
        f"<td colspan='2'><b>{DELTA}G_gas + {DELTA}G_solv.</b> Use this number for relative ranking "
        f"across your screen.</td></tr>",
    ]

    try:
        chart_b64 = make_waterfall_png_b64(delta)
        chart_html = f'<img class="chart" src="data:image/png;base64,{chart_b64}" alt="Energy decomposition chart">'
    except Exception as e:
        chart_html = f"<p class='meta'><i>(Chart unavailable: {html.escape(str(e))})</i></p>"

    ligand_name = Path(src_path).parent.parent.name or "(unknown)"
    title = f"MM/GBSA Report - {ligand_name}"

    band = ('strong micro-to-nanomolar' if dG < -30
            else 'micro-to-millimolar'  if dG < -15
            else 'very weak / no')

    # Hoist out of f-string expressions (no backslashes inside {...} on 3.11).
    frames_cell         = meta.get("frames", EMDASH)
    igb_cell            = meta.get("igb",    EMDASH)
    salt_cell           = meta.get("saltcon", EMDASH)
    term_rows_joined    = "".join(term_rows_html)
    summary_rows_joined = "".join(summary_rows_html)
    gen_time            = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    return f"""<!doctype html>
<html><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>{CSS}</style></head>
<body>

<h1>MM/GBSA Binding-Energy Report</h1>
<p class="meta">
  Source file: <code>{html.escape(str(src_path))}</code><br>
  Ligand workdir: <code>{html.escape(ligand_name)}</code><br>
  Generated: {gen_time}
</p>

<div class="hero" style="background:{vcolor};">
  {html.escape(verdict)}
  <span class="big">{DELTA}G<sub>bind</sub> = {dG:+.2f} {PLUSMIN} {sem:.2f} kcal/mol</span>
</div>

<h2>Headline numbers</h2>
<table>
  <tr><th>Quantity</th><th>Value</th><th>Interpretation</th></tr>
  <tr><td><b>{DELTA}G_bind</b> (DELTA TOTAL)</td>
      <td class="num">{dG:+.2f} kcal/mol</td>
      <td>Predicted binding free energy</td></tr>
  <tr><td>Std. dev. across frames</td>
      <td class="num">{PLUSMIN}{std:.2f} kcal/mol</td>
      <td>Frame-to-frame fluctuation in the trajectory</td></tr>
  <tr><td>Std. err. of mean</td>
      <td class="num">{PLUSMIN}{sem:.2f} kcal/mol</td>
      <td>Statistical precision of the average (lower = tighter ranking)</td></tr>
  <tr><td>Frames analyzed</td>
      <td class="num">{frames_cell}</td>
      <td>Production-trajectory snapshots used</td></tr>
  <tr><td>Method</td>
      <td class="num">igb={igb_cell}, salt={salt_cell} M</td>
      <td>Single-trajectory MM/GBSA</td></tr>
</table>

<h2>Term-by-term decomposition</h2>
{chart_html}

<table>
  <tr><th>Term</th><th>{DELTA} (kcal/mol)</th><th>What it means</th><th>Read</th></tr>
  {term_rows_joined}
</table>

<h2>Aggregated free-energy components</h2>
<table>
  <tr><th>Component</th><th>Value</th><th colspan="2">Notes</th></tr>
  {summary_rows_joined}
</table>

<div class="caveat">
  <b>Caveat on absolute magnitude.</b> MM/GBSA {DELTA}G_bind values from <code>igb=8</code>
  are typically 2-5x more negative than experimental {DELTA}G_bind, because the method
  ignores configurational entropy (no normal-mode step), uses implicit solvent,
  and inherits force-field / GB-radii systematics. Use this number for
  <b>relative ranking within a screen</b>, not as a literal absolute Kd.
  A {DELTA}G_bind of {dG:+.1f} kcal/mol corresponds roughly to a
  <b>{band}</b> binder experimentally, but calibrate against known actives if available.
</div>

<h2>How to use this number</h2>
<table>
  <tr><th>Comparison</th><th>What it means</th></tr>
  <tr><td>vs. other ligands in this screen</td>
      <td>Differences {PLUSMIN}1 kcal/mol are statistically meaningful (SEM={sem:.2f}).
          Use to rank-order; not for absolute Kd.</td></tr>
  <tr><td>vs. a known active</td>
      <td>If you have an active with measured Kd, compare its MM/GBSA score to
          this one - the offset is your calibration.</td></tr>
  <tr><td>vs. experimental {DELTA}G_bind</td>
      <td><b>Not directly comparable.</b> Magnitude inflated; sign is informative.</td></tr>
</table>

<div class="footer">
  Generated by <code>mmpbsa_report.py</code> from the AffiniScreen.
  {DELTA}G_bind = {DELTA}G_gas + {DELTA}G_solv (single-trajectory MM/GBSA, no entropy correction).
</div>

</body></html>
"""


# ----------------------------- importable API -----------------------------

def generate_report(dat_path, out_html=None):
    """Generate one HTML report from one FINAL_RESULTS_MMPBSA.dat file."""
    dat_path = Path(dat_path)
    if not dat_path.exists():
        raise FileNotFoundError(dat_path)
    parsed = parse_mmpbsa_dat(dat_path)
    delta = parsed["blocks"].get("delta", {})
    if "DELTA TOTAL" not in delta:
        raise RuntimeError(f"No DELTA TOTAL block in {dat_path} - MMPBSA may have failed")
    html_out = build_html(parsed, dat_path)
    if out_html is None:
        out_html = dat_path.with_name("FINAL_RESULTS.report.html")
    out_html = Path(out_html)
    out_html.write_text(html_out)
    return out_html


# ----------------------------- CLI driver -----------------------------

def resolve_targets(arg):
    arg = Path(arg).expanduser().resolve()
    if arg.is_file() and arg.suffix.lower() == ".dat":
        # final49: any .dat file is accepted directly (covers _MMGBSA naming).
        yield arg, arg.with_name("FINAL_RESULTS.report.html"); return
    if arg.is_dir():
        # final49: accept either FINAL_RESULTS_MMPBSA.dat (pipeline) or
        # FINAL_RESULTS_MMGBSA.dat (older/manual runs). Prefer MMPBSA if both.
        def _find_dat(d):
            for nm in ("FINAL_RESULTS_MMPBSA.dat", "FINAL_RESULTS_MMGBSA.dat"):
                p = d / nm
                if p.exists():
                    return p
            return None
        single = _find_dat(arg / "mmgbsa")
        if single is None:
            single = _find_dat(arg)        # also tolerate dat sitting in workdir root
        if single is not None:
            yield single, single.with_name("FINAL_RESULTS.report.html"); return
        batch_hits = sorted(arg.glob("lig_*/mmgbsa/FINAL_RESULTS_MMPBSA.dat")) \
            + sorted(arg.glob("lig_*/mmgbsa/FINAL_RESULTS_MMGBSA.dat"))
        if batch_hits:
            for hit in batch_hits:
                yield hit, hit.with_name("FINAL_RESULTS.report.html")
            return
    sys.exit(f"ERROR: cannot find FINAL_RESULTS_MM[PG]BSA.dat at or under {arg}")

def main():
    ap = argparse.ArgumentParser(
        description="Generate a styled HTML report from FINAL_RESULTS_MMPBSA.dat")
    ap.add_argument("path", type=Path,
                    help="A .dat file, a single-ligand workdir, or a batch dir.")
    ap.add_argument("--print-summary", action="store_true")
    a = ap.parse_args()

    n = 0
    for dat, out in resolve_targets(a.path):
        try:
            html_path = generate_report(dat, out)
            n += 1
            print(f"[OK] {html_path}")
            if a.print_summary:
                parsed = parse_mmpbsa_dat(dat)
                dT = parsed["blocks"].get("delta", {}).get("DELTA TOTAL", {})
                dG, sem = dT.get("avg"), dT.get("sem")
                if dG is not None:
                    print(f"     {DELTA}G_bind = {dG:+.2f} +/- {sem:.2f} kcal/mol")
        except Exception as e:
            print(f"[FAIL] {dat}  ({e})", file=sys.stderr)
    print(f"\nWrote {n} report(s).")

if __name__ == "__main__":
    main()
