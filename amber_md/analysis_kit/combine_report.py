#!/usr/bin/env python3
"""combine_report.py  (v3 - 2026-05-26)

Per-ligand combined HTML.

v3 changes:
  * NEW: parse_mmpbsa_dat() reads dG / std / SEM directly from
    mmgbsa/FINAL_RESULTS_MMPBSA.dat (the raw MMPBSA.py output), which is
    written by MMPBSA.py the instant the calculation finishes.
  * dG source order: raw .dat first, HTML report second, None last.
    This decouples the combined report from the timing/success of the
    HTML report generator -- previously the kit would render with
    "dG = n/a" if it ran in the same minute MM/GBSA finished, because
    the HTML report was rendered by a separate Python step that hadn't
    fired yet (and on v2.5.1 and earlier, was crashing on PYTHONPATH).
  * One-line stderr log: `[combine_report] dG source = dat|html|none`
    so this is debuggable next time without grepping HTML.
  * extract_mmgbsa() unchanged -- still used as the fallback and still
    the source for per-term VDWAALS/EEL/EGB/ESURF decomposition.

Inputs auto-discovered under <ligand_workdir>:
  mmgbsa/FINAL_RESULTS_MMPBSA.dat       <-- NEW preferred source for dG
  mmgbsa/FINAL_RESULTS.report.html      <-- fallback for dG; only source
                                            for per-term decomposition
  complex_solv_view.html
  analysis/rmsd_backbone.dat   rmsd_ligand.dat   rmsf_byres.dat
  analysis/hbond_avg.dat       hbond_count.dat
  analysis/contacts_byres.dat  contacts_series.dat (optional)

Outputs:
  analysis/rmsd.png  rmsf.png  hbond_chart.png  contact_heatmap.png
  analysis/COMBINED_REPORT.html
  analysis/summary.json   <-- consumed by build_screen_summary.py
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import base64, io, json, re, sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_dat(p: Path) -> np.ndarray:
    return np.loadtxt(p, comments="#") if p.exists() and p.stat().st_size else np.empty((0, 0))


def fig_to_uri(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def fig_uri_from_path(p: Path) -> str:
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()


def plot_rmsd(bb, lig, png):
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.plot(bb[:, 0], bb[:, 1], lw=1.4, label="Protein backbone")
    if lig.size:
        ax.plot(lig[:, 0], lig[:, 1], lw=1.4, color="#c92a2a", label="Ligand heavy atoms")
    ax.set_xlabel("Frame"); ax.set_ylabel("RMSD (A)")
    ax.set_title("RMSD vs. first frame"); ax.grid(alpha=.3); ax.legend(frameon=False)
    fig.savefig(png, dpi=130, bbox_inches="tight")
    return fig_to_uri(fig)


def plot_rmsf(rmsf, png):
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.bar(rmsf[:, 0], rmsf[:, 1], width=1.0, color="#1971c2")
    ax.set_xlabel("Residue"); ax.set_ylabel("RMSF (A)")
    ax.set_title("Per-residue Calpha fluctuation"); ax.grid(alpha=.3, axis="y")
    fig.savefig(png, dpi=130, bbox_inches="tight")
    return fig_to_uri(fig)


def plot_hbond_bars(hb_rows, png):
    if not hb_rows:
        return None
    rows = sorted(hb_rows, key=lambda r: -r[1])[:10]
    labels = [r[0] for r in rows]; fracs = [r[1] for r in rows]
    fig, ax = plt.subplots(figsize=(8, max(2.4, 0.35 * len(rows))))
    ax.barh(range(len(rows)), fracs, color="#2b8a3e")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Occupancy (fraction of frames)")
    ax.set_xlim(0, 1); ax.set_title("Top-10 ligand H-bonds")
    ax.grid(alpha=.3, axis="x")
    fig.savefig(png, dpi=130, bbox_inches="tight")
    return fig_to_uri(fig)


def plot_contact_heatmap(series_path, png):
    """contacts_series.dat: col0=frame, col1..N = per-contact 0/1.
    Aggregate to per-residue fraction by parsing header labels like :12@OD1_:LIG@N1."""
    if not series_path.exists() or series_path.stat().st_size == 0:
        return None, []
    with open(series_path) as fh:
        header = fh.readline().lstrip("#").split()
    data = np.loadtxt(series_path, comments="#")
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 2:
        return None, []
    frames = data[:, 0]; cols = data[:, 1:]
    res_for_col = []
    for name in header[1:]:
        m = re.search(r":(\d+)", name)
        res_for_col.append(int(m.group(1)) if m else -1)
    residues = sorted({r for r in res_for_col if r > 0})
    if not residues:
        return None, []
    occ = np.zeros((len(residues), data.shape[0]))
    for i, r in enumerate(residues):
        idx = [j for j, rr in enumerate(res_for_col) if rr == r]
        if idx:
            occ[i] = (cols[:, idx] > 0).any(axis=1).astype(float)
    keep = occ.mean(axis=1) > 0.05
    residues = [r for r, k in zip(residues, keep) if k]
    occ = occ[keep]
    if not residues:
        return None, []
    fig, ax = plt.subplots(figsize=(9, max(2.4, 0.22 * len(residues))))
    im = ax.imshow(occ, aspect="auto", cmap="viridis",
                   extent=[frames.min(), frames.max(), len(residues), 0],
                   vmin=0, vmax=1, interpolation="nearest")
    ax.set_yticks(np.arange(len(residues)) + 0.5)
    ax.set_yticklabels([f"R{r}" for r in residues], fontsize=8)
    ax.set_xlabel("Frame"); ax.set_title("Ligand-residue contact occupancy (5 A)")
    fig.colorbar(im, ax=ax, label="contact (0/1)")
    fig.savefig(png, dpi=130, bbox_inches="tight"); plt.close(fig)
    summary = sorted(zip(residues, occ.mean(axis=1).tolist()), key=lambda x: -x[1])
    return fig_uri_from_path(png), summary


def parse_hbond_avg(path: Path):
    """cpptraj hbond avgout: cols = Acceptor DonorH Donor Frames Frac AvgDist AvgAng"""
    if not path.exists():
        return []
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        acc, donH, don = parts[0], parts[1], parts[2]
        try:
            frac = float(parts[4])
            dist = float(parts[5])
            ang  = float(parts[6]) if len(parts) > 6 else float("nan")
        except ValueError:
            continue
        prot = acc if ":LIG" not in acc else don
        lig  = don if ":LIG" in don else acc
        rows.append((f"{prot}  ...  {lig}", frac, dist, ang))
    return rows


# ----------------------------------------------------------------------
# v3: MM/GBSA dG -- raw .dat preferred, HTML report fallback.
# ----------------------------------------------------------------------

def parse_mmpbsa_dat(path: Path) -> dict:
    """Read dG / std / SEM straight from FINAL_RESULTS_MMPBSA.dat.

    Format from MMPBSA.py looks like:
        DELTA TOTAL           -29.1266       2.5049       0.2505
    Same regex shape as amber_md.batch_aggregate.parse_mmpbsa_dat, so
    anything that ranks across ligands and the combined report agree
    by construction.

    Returns dict with keys 'dG', 'std', 'sem'. Empty dict on miss --
    callers should treat that as "fall through to the HTML report".
    Never raises.
    """
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        text = path.read_text(errors="replace")
    except Exception:
        return {}
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", text)
    if m:
        return {"dG":  float(m.group(1)),
                "std": float(m.group(2)),
                "sem": float(m.group(3))}
    # Looser fallback for unusual MMPBSA.py output formats.
    m = re.search(r"Delta\s+Total\s+Energy\s*[:=]?\s*(-?\d+\.\d+)",
                  text, re.IGNORECASE)
    if m:
        return {"dG":  float(m.group(1)),
                "std": None,
                "sem": None}
    return {}


def extract_mmgbsa(html: str) -> dict:
    """Parse dG / SEM / per-term decomposition out of the v2.4.8+ HTML
    report. Kept for two reasons:
      1. Fallback source for dG when the raw .dat isn't there (legacy runs).
      2. Sole source for per-term decomposition (VDWAALS, EEL, EGB, ESURF).
    """
    out = {}
    m = re.search(r"&Delta;G<sub>bind</sub>\s*=\s*(-?\d+\.\d+)\s*(?:&plusmn;|\xb1)\s*(\d+\.\d+)", html)
    if not m:
        m = re.search(r"\u0394G<sub>bind</sub>\s*=\s*(-?\d+\.\d+)\s*\u00b1\s*(\d+\.\d+)", html)
    if m:
        out["dG"], out["sem"] = float(m.group(1)), float(m.group(2))
    for term in ("VDWAALS", "EEL", "EGB", "ESURF"):
        m = re.search(rf"<b>{term}</b>.*?(-?\+?\d+\.\d+)", html, re.S)
        if m:
            out[term] = float(m.group(1).lstrip("+"))
    return out


def _resolve_mmgbsa(ligand_dir: Path) -> tuple[dict, str]:
    """Return (energy_dict, source_tag).
    source_tag is 'dat' | 'html' | 'merged' | 'none' and is logged to stderr.

    Order of preference:
      1. raw .dat for dG/std/SEM (always available right after MMPBSA.py).
      2. HTML report for whatever .dat couldn't supply, especially the
         per-term decomposition (VDWAALS/EEL/EGB/ESURF).
    If both are present we MERGE: dG/std/SEM from .dat, terms from HTML.
    """
    dat_path  = ligand_dir / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
    html_path = ligand_dir / "mmgbsa" / "FINAL_RESULTS.report.html"

    dat_dict = parse_mmpbsa_dat(dat_path)
    html_dict = {}
    if html_path.exists():
        try:
            html_dict = extract_mmgbsa(html_path.read_text(errors="replace"))
        except Exception:
            html_dict = {}

    if dat_dict and html_dict:
        merged = dict(html_dict)          # start with terms etc.
        merged.update(dat_dict)           # raw .dat wins for dG / std / sem
        return merged, "merged"
    if dat_dict:
        return dat_dict, "dat"
    if html_dict:
        return html_dict, "html"
    return {}, "none"


def build(ligand_dir: Path) -> Path:
    A = ligand_dir / "analysis"; A.mkdir(exist_ok=True)

    bb   = load_dat(A / "rmsd_backbone.dat")
    lig  = load_dat(A / "rmsd_ligand.dat")
    rmsf = load_dat(A / "rmsf_byres.dat")

    rmsd_uri = plot_rmsd(bb, lig, A / "rmsd.png") if bb.size else None
    rmsf_uri = plot_rmsf(rmsf,    A / "rmsf.png") if rmsf.size else None

    hb_rows  = parse_hbond_avg(A / "hbond_avg.dat")
    hb_uri   = plot_hbond_bars([(l, f) for l, f, *_ in hb_rows], A / "hbond_chart.png")

    cmap_uri, contact_summary = plot_contact_heatmap(A / "contacts_series.dat",
                                                     A / "contact_heatmap.png")

    # v3: raw .dat first, HTML fallback, merged when both present.
    mm, mm_source = _resolve_mmgbsa(ligand_dir)
    sys.stderr.write(f"[combine_report] dG source = {mm_source}\n")
    mm_path = ligand_dir / "mmgbsa" / "FINAL_RESULTS.report.html"

    bb_mean  = float(np.mean(bb[:, 1])) if bb.size else None
    bb_max   = float(np.max(bb[:, 1]))  if bb.size else None
    lig_mean = float(np.mean(lig[:, 1])) if lig.size else None
    rmsf_top = sorted(((int(r), float(v)) for r, v in rmsf), key=lambda x: -x[1])[:5] if rmsf.size else []

    summary = {
        "ligand":          ligand_dir.name,
        "dG":              mm.get("dG"),
        "std":             mm.get("std"),      # v3: now populated from raw .dat
        "sem":             mm.get("sem"),
        "dG_source":       mm_source,          # v3: 'dat' | 'html' | 'merged' | 'none'
        "VDWAALS":         mm.get("VDWAALS"),
        "EEL":             mm.get("EEL"),
        "EGB":             mm.get("EGB"),
        "ESURF":           mm.get("ESURF"),
        "rmsd_bb_mean":    bb_mean,
        "rmsd_bb_max":     bb_max,
        "rmsd_lig_mean":   lig_mean,
        "n_hbonds_strong": sum(1 for _, f, *_ in hb_rows if f >= 0.5),
        "top_contacts":    contact_summary[:5],
        "report_path":     f"{ligand_dir.name}/analysis/COMBINED_REPORT.html",
    }
    (A / "summary.json").write_text(json.dumps(summary, indent=2))

    # dG hero string: prefer SEM if present, fall back to std, then bare value.
    if summary["dG"] is None:
        dG_html = "n/a"
    elif summary.get("sem") is not None:
        dG_html = f"{summary['dG']:+.2f} +/- {summary['sem']:.2f} kcal/mol"
    elif summary.get("std") is not None:
        dG_html = f"{summary['dG']:+.2f} +/- {summary['std']:.2f} kcal/mol (std)"
    else:
        dG_html = f"{summary['dG']:+.2f} kcal/mol"

    term_rows = "".join(
        f"<tr><td>{t}</td><td class='num'>{mm[t]:+.2f}</td></tr>"
        for t in ("VDWAALS", "EEL", "EGB", "ESURF") if t in mm)
    rmsf_rows = "".join(f"<tr><td>{r}</td><td class='num'>{v:.2f}</td></tr>" for r, v in rmsf_top)
    hb_table_rows = "".join(
        f"<tr><td>{l}</td><td class='num'>{f:.2f}</td>"
        f"<td class='num'>{d:.2f}</td><td class='num'>{a:.1f}</td></tr>"
        for l, f, d, a in sorted(hb_rows, key=lambda r: -r[1])[:15]
    ) or "<tr><td colspan=4><i>No persistent H-bonds detected.</i></td></tr>"
    contact_rows = "".join(
        f"<tr><td>Res {r}</td><td class='num'>{occ:.2f}</td></tr>"
        for r, occ in contact_summary[:10]
    ) or "<tr><td colspan=2><i>No persistent contacts &ge; 5%.</i></td></tr>"

    ngl_link = ('<a href="../complex_solv_view.html">3-D viewer (NGL)</a>'
                if (ligand_dir / "complex_solv_view.html").exists() else "")
    mm_link  = ('<a href="../mmgbsa/FINAL_RESULTS.report.html">Full MM/GBSA report</a>'
                if mm_path.exists() else "")
    hb_img      = f'<img class="chart" src="{hb_uri}">'    if hb_uri    else "<p><i>No H-bond chart.</i></p>"
    contact_img = f'<img class="chart" src="{cmap_uri}">'  if cmap_uri  else "<p><i>No contact map.</i></p>"
    rmsd_img    = f'<img class="chart" src="{rmsd_uri}">'  if rmsd_uri  else ""
    rmsf_img    = f'<img class="chart" src="{rmsf_uri}">'  if rmsf_uri  else ""

    # v3: badge the dG hero with the source so users can see at a glance
    # whether the report fell back to a stale HTML or had no data at all.
    hero_bg = "#2b8a3e" if summary["dG"] is not None else "#868e96"
    src_badge = f" <small style='opacity:.75;font-size:.55em'>[source: {mm_source}]</small>"

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Combined report - {ligand_dir.name}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1050px;margin:24px auto;padding:0 20px;color:#222;line-height:1.5}}
 h1{{border-bottom:3px solid #1971c2;padding-bottom:6px}}
 h2{{color:#1971c2;margin-top:28px}}
 table{{border-collapse:collapse;width:100%;margin:10px 0 22px;font-size:.94em}}
 th,td{{padding:6px 10px;border-bottom:1px solid #e0e0e0;text-align:left;vertical-align:top}}
 th{{background:#f5f7fa}} td.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
 .grid{{display:grid;grid-template-columns:1fr 1fr;gap:20px}}
 .card{{border:1px solid #e0e0e0;border-radius:8px;padding:14px 18px;background:#fafbfc}}
 .hero{{padding:18px 22px;border-radius:8px;color:white;background:{hero_bg};font-size:1.1em}}
 .hero .big{{font-size:1.8em;font-weight:700;display:block;margin-top:4px}}
 img.chart{{max-width:100%;border:1px solid #e0e0e0;border-radius:6px;padding:6px;background:white}}
 code{{background:#f1f3f5;padding:1px 5px;border-radius:3px}}
 .btn{{display:inline-block;margin:4px 8px 4px 0;padding:6px 14px;border-radius:5px;background:#1971c2;color:white;text-decoration:none;font-size:.9em}}
 .btn:hover{{background:#1864ab}}
</style></head><body>

<h1>Combined MD report - <code>{ligand_dir.name}</code></h1>
<div class="hero">Predicted binding free energy{src_badge}
  <span class="big">&Delta;G<sub>bind</sub> = {dG_html}</span></div>

<h2>Quick links</h2>
<p>{ngl_link} &nbsp;&middot;&nbsp; {mm_link} &nbsp;&middot;&nbsp;
   <a href="../../screen_summary.html">screen summary</a><br>
  <a class="btn" href="load_pymol.pml">PyMOL session</a>
  <a class="btn" href="load_vmd.tcl">VMD session</a></p>

<div class="grid">
  <div class="card"><h2 style="margin-top:0">Trajectory stability</h2>
    <table><tr><th>Metric</th><th>Value</th></tr>
      <tr><td>Mean backbone RMSD</td><td class='num'>{bb_mean:.2f} A</td></tr>
      <tr><td>Max  backbone RMSD</td><td class='num'>{bb_max:.2f} A</td></tr>
      {'<tr><td>Mean ligand RMSD</td><td class="num">%.2f A</td></tr>'%lig_mean if lig_mean is not None else ''}
      <tr><td>Persistent H-bonds (occ &ge; 0.5)</td><td class='num'>{summary['n_hbonds_strong']}</td></tr>
    </table></div>
  <div class="card"><h2 style="margin-top:0">Top-5 most flexible residues</h2>
    <table><tr><th>Residue #</th><th>RMSF (A)</th></tr>{rmsf_rows}</table></div>
</div>

<h2>RMSD over time</h2>{rmsd_img}
<h2>Per-residue RMSF</h2>{rmsf_img}

<h2>Ligand H-bonds (top 15 by occupancy)</h2>
<table><tr><th>Pair (protein ... ligand)</th><th>Occupancy</th>
  <th>Avg dist (A)</th><th>Avg angle (deg)</th></tr>{hb_table_rows}</table>
{hb_img}

<h2>Ligand-residue contact occupancy</h2>
{contact_img}
<table><tr><th>Residue</th><th>Occupancy</th></tr>{contact_rows}</table>

<h2>MM/GBSA term decomposition</h2>
<table><tr><th>Term</th><th>&Delta; (kcal/mol)</th></tr>
{term_rows or '<tr><td colspan=2><i>n/a</i></td></tr>'}</table>

<p style="font-size:.8em;color:#888;border-top:1px solid #e0e0e0;padding-top:10px">
Generated by <code>combine_report.py</code> v3.
</p></body></html>"""

    out_path = A / "COMBINED_REPORT.html"
    out_path.write_text(html)
    return out_path


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: combine_report.py <ligand_workdir>")
    p = build(Path(sys.argv[1]).resolve())
    print(f"[combine_report] wrote {p}")
