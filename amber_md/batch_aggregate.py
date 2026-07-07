"""Aggregate DG_bind across all ligands in a batch directory (v2.5.1).

v2.5.1: locate analysis_kit and mmpbsa_report INSIDE the amber_md package
first (matching the current source layout), with a fallback to the legacy
parent.parent location.

v2.5.0: also runs analysis_kit/run_screen_analysis.sh and
build_screen_summary.py at the start.

Usage:
    python -m amber_md.batch_aggregate <batch_dir>
"""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import sys, re, html, base64, io, datetime, importlib.util, subprocess
from pathlib import Path


# ----------------------------- path discovery (v2.5.1) -----------------------------

def _find_analysis_kit():
    pkg_dir = Path(__file__).resolve().parent
    for cand in (pkg_dir / "analysis_kit", pkg_dir.parent / "analysis_kit"):
        if cand.is_dir():
            return cand
    return None


def _find_mmpbsa_report_py():
    pkg_dir = Path(__file__).resolve().parent
    for cand in (pkg_dir / "mmpbsa_report.py", pkg_dir.parent / "mmpbsa_report.py"):
        if cand.exists():
            return cand
    return None


# ----------------------------- parsing -----------------------------

def parse_mmpbsa_dat(path):
    if not path.exists(): return None
    text = path.read_text()
    m = re.search(r"DELTA TOTAL\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", text)
    if m:
        return {"delta_total": float(m.group(1)),
                "std":         float(m.group(2)),
                "sem":         float(m.group(3))}
    m = re.search(r"Delta\s+Total\s+Energy\s*[:=]?\s*(-?\d+\.\d+)", text, re.IGNORECASE)
    if m:
        return {"delta_total": float(m.group(1)), "std": None, "sem": None}
    return None


# ----------------------------- per-ligand report fallback -----------------------------

def _load_mmpbsa_report_module():
    """v2.5.1: package import first, then file discovery."""
    try:
        from . import mmpbsa_report  # type: ignore
        return mmpbsa_report
    except Exception:
        pass

    candidate = _find_mmpbsa_report_py()
    if candidate is None:
        return None
    try:
        spec = importlib.util.spec_from_file_location("mmpbsa_report", candidate)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["mmpbsa_report"] = mod
        spec.loader.exec_module(mod)
        return mod
    except Exception as e:
        sys.stderr.write(f"[batch_aggregate] cannot import mmpbsa_report.py: {e}\n")
        return None


def _ensure_per_ligand_report(dat_path, report_module):
    if report_module is None: return None
    if not dat_path.exists() or dat_path.stat().st_size == 0: return None
    html_path = dat_path.with_name("FINAL_RESULTS.report.html")
    if html_path.exists() and html_path.stat().st_size > 0:
        return html_path
    try:
        return report_module.generate_report(dat_path, html_path)
    except Exception as e:
        sys.stderr.write(f"[batch_aggregate] report failed for {dat_path}: {e}\n")
        return None


# ----------------------------- screen-wide analysis_kit -----------------------------

def _run_screen_analysis_kit(batch_dir):
    """v2.5.1: ak path via _find_analysis_kit() so in-package layout works."""
    ak = _find_analysis_kit()
    if ak is None:
        sys.stderr.write("[batch_aggregate] analysis_kit not found - "
                         "skipping screen analysis & summary.\n")
        return

    sh = ak / "run_screen_analysis.sh"
    py = ak / "build_screen_summary.py"

    if sh.exists():
        sys.stderr.write(f"[batch_aggregate] running screen analysis kit on {batch_dir}\n")
        try:
            subprocess.run(["bash", str(sh), str(batch_dir)],
                           cwd=str(ak), check=False, timeout=3600)
        except subprocess.TimeoutExpired:
            sys.stderr.write("[batch_aggregate] run_screen_analysis.sh timed out (3600 s)\n")
        except Exception as e:
            sys.stderr.write(f"[batch_aggregate] screen analysis failed: {e}\n")
    else:
        sys.stderr.write(f"[batch_aggregate] {sh} not found - skipping screen analysis.\n")

    if py.exists():
        sys.stderr.write("[batch_aggregate] building screen summary\n")
        try:
            subprocess.run([sys.executable, str(py), str(batch_dir)],
                           cwd=str(ak), check=False, timeout=600)
        except subprocess.TimeoutExpired:
            sys.stderr.write("[batch_aggregate] build_screen_summary.py timed out (600 s)\n")
        except Exception as e:
            sys.stderr.write(f"[batch_aggregate] screen summary failed: {e}\n")
    else:
        sys.stderr.write(f"[batch_aggregate] {py} not found - skipping screen summary.\n")


# ----------------------------- INDEX.html -----------------------------

INDEX_CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
       max-width: 1100px; margin: 24px auto; padding: 0 20px; color: #222; line-height: 1.5; }
h1 { border-bottom: 3px solid #1971c2; padding-bottom: 6px; }
h2 { color: #1971c2; margin-top: 28px; }
.meta { font-size: 0.88em; color: #666; }
.kpis { display: flex; gap: 18px; margin: 16px 0 22px; flex-wrap: wrap; }
.kpi { padding: 12px 18px; border: 1px solid #e0e0e0; border-radius: 8px;
       background: #f8f9fa; min-width: 140px; }
.kpi .label { font-size: 0.78em; color: #888; text-transform: uppercase;
              letter-spacing: 0.04em; }
.kpi .value { font-size: 1.5em; font-weight: 700; margin-top: 2px; }
table { border-collapse: collapse; width: 100%; margin: 10px 0 22px; font-size: 0.93em; }
th, td { padding: 7px 10px; border-bottom: 1px solid #e0e0e0; text-align: left;
         vertical-align: middle; }
th { background: #f5f7fa; font-weight: 600; cursor: pointer; user-select: none; }
th.sort-asc::after  { content: "  \u25b2"; color: #1971c2; }
th.sort-desc::after { content: "  \u25bc"; color: #1971c2; }
tr:nth-child(even) td { background: #fafbfc; }
td.num { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.rank { text-align: right; color: #888; font-variant-numeric: tabular-nums; }
.good    { color: #2b8a3e; font-weight: 600; }
.bad     { color: #c92a2a; font-weight: 600; }
.status-done    { color: #2b8a3e; font-weight: 600; }
.status-pending { color: #f08c00; font-weight: 600; }
.status-failed  { color: #c92a2a; font-weight: 600; }
a { color: #1971c2; text-decoration: none; }
a:hover { text-decoration: underline; }
.footer { font-size: 0.8em; color: #888; border-top: 1px solid #e0e0e0;
          padding-top: 10px; margin-top: 30px; }
img.chart { max-width: 100%; border: 1px solid #e0e0e0; border-radius: 6px;
            padding: 8px; background: white; }
"""

INDEX_JS = """
document.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('table.sortable th').forEach((th, colIdx) => {
    th.addEventListener('click', () => {
      const table = th.closest('table');
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.rows);
      const isNumeric = th.classList.contains('numeric');
      const dir = th.classList.contains('sort-asc') ? -1 : 1;
      table.querySelectorAll('th').forEach(h => h.classList.remove('sort-asc','sort-desc'));
      th.classList.add(dir === 1 ? 'sort-asc' : 'sort-desc');
      rows.sort((a, b) => {
        const av = a.cells[colIdx].dataset.sortKey ?? a.cells[colIdx].textContent.trim();
        const bv = b.cells[colIdx].dataset.sortKey ?? b.cells[colIdx].textContent.trim();
        const aMissing = (av === '' || av === '--');
        const bMissing = (bv === '' || bv === '--');
        if (aMissing && bMissing) return 0;
        if (aMissing) return 1;
        if (bMissing) return -1;
        if (isNumeric) return (parseFloat(av) - parseFloat(bv)) * dir;
        return av.localeCompare(bv) * dir;
      });
      rows.forEach(r => tbody.appendChild(r));
    });
  });
});
"""

def _make_ranking_chart_b64(done_rows):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    if not done_rows: return None
    rows = sorted(done_rows, key=lambda r: r["result"]["delta_total"])
    labels = [_short_ligand_label(r["workdir"]) for r in rows]
    vals   = [r["result"]["delta_total"] for r in rows]
    errs   = [r["result"].get("sem") or 0.0 for r in rows]
    h = max(2.6, 0.34 * len(rows) + 1.0)
    fig, ax = plt.subplots(figsize=(8.5, h))
    y = list(range(len(rows)))[::-1]
    colors = ["#2b8a3e" if v < 0 else "#c92a2a" for v in vals]
    ax.barh(y, vals, xerr=errs, color=colors, ecolor="#444",
            capsize=3, edgecolor="black", linewidth=0.5)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=9)
    ax.axvline(0, color="black", lw=0.7)
    ax.set_xlabel("\u0394G_bind  (kcal/mol)", fontsize=10)
    ax.set_title(f"Ranked binding free energies ({len(rows)} ligand(s))",
                 fontsize=11, weight="bold")
    ax.grid(axis="x", alpha=0.25, linestyle="--")
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=130, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _short_ligand_label(name, maxlen=42):
    if len(name) <= maxlen: return name
    return name[:maxlen-1] + "\u2026"


def _build_index_html(batch_dir, rows):
    done    = [r for r in rows if r["result"]]
    pending = [r for r in rows if not r["result"]]
    done.sort(key=lambda r: r["result"]["delta_total"])

    n_total   = len(rows)
    n_done    = len(done)
    n_pending = len(pending)

    best_dG = done[0]["result"]["delta_total"] if done else None
    mean_dG = (sum(r["result"]["delta_total"] for r in done) / n_done) if done else None

    chart_b64 = _make_ranking_chart_b64(done)
    if chart_b64:
        chart_html = (f'<img class="chart" src="data:image/png;base64,{chart_b64}" '
                      f'alt="Ranking chart">')
    else:
        chart_html = "<p class='meta'><i>(Chart unavailable: matplotlib missing or no DONE ligands.)</i></p>"

    screen_html = batch_dir / "screen_summary.html"
    screen_link = ""
    if screen_html.exists():
        screen_link = ('<p class="meta">Screen-wide cpptraj summary: '
                       '<a href="screen_summary.html">screen_summary.html</a></p>')

    ranked_rows_html = []
    for rank, r in enumerate(done, 1):
        wd = r["workdir"]
        res = r["result"]
        dG  = res["delta_total"]
        std = res.get("std"); sem = res.get("sem")
        rep = r.get("report_path")
        if rep:
            try: rel = rep.relative_to(batch_dir)
            except ValueError: rel = rep
            wd_cell = f'<a href="{html.escape(str(rel))}">{html.escape(wd)}</a>'
        else:
            wd_cell = html.escape(wd)
        std_txt = f"{std:.2f}" if std is not None else "\u2014"
        sem_txt = f"{sem:.2f}" if sem is not None else "\u2014"
        cls = "good" if dG < 0 else "bad"
        std_key = std if std is not None else 9e9
        sem_key = sem if sem is not None else 9e9
        ranked_rows_html.append(
            f"<tr>"
            f'<td class="rank" data-sort-key="{rank}">{rank}</td>'
            f"<td>{wd_cell}</td>"
            f'<td class="num {cls}" data-sort-key="{dG}">{dG:+.2f}</td>'
            f'<td class="num" data-sort-key="{std_key}">{std_txt}</td>'
            f'<td class="num" data-sort-key="{sem_key}">{sem_txt}</td>'
            f'<td><span class="status-done">DONE</span></td>'
            f"</tr>")

    pending_rows_html = []
    for r in pending:
        wd = r["workdir"]
        status = r["status"]
        sclass = "status-failed" if "fail" in status.lower() else "status-pending"
        pending_rows_html.append(
            f"<tr>"
            f'<td class="rank">\u2014</td>'
            f"<td>{html.escape(wd)}</td>"
            f'<td class="num" data-sort-key="9e9">\u2014</td>'
            f'<td class="num" data-sort-key="9e9">\u2014</td>'
            f'<td class="num" data-sort-key="9e9">\u2014</td>'
            f'<td><span class="{sclass}">{html.escape(status)}</span></td>'
            f"</tr>")

    kpi_best = f"{best_dG:+.2f}" if best_dG is not None else "\u2014"
    kpi_mean = f"{mean_dG:+.2f}" if mean_dG is not None else "\u2014"
    gen_time = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    name_esc = html.escape(batch_dir.name)
    path_esc = html.escape(str(batch_dir))
    ranked_joined  = "".join(ranked_rows_html)
    pending_joined = "".join(pending_rows_html)

    return f"""<!doctype html>
<html><head><meta charset="utf-8">
<title>Batch index - {name_esc}</title>
<style>{INDEX_CSS}</style>
<script>{INDEX_JS}</script>
</head><body>

<h1>Batch screen: <code>{name_esc}</code></h1>
<p class="meta">
  Path: <code>{path_esc}</code><br>
  Generated: {gen_time}
</p>

<div class="kpis">
  <div class="kpi"><div class="label">Ligands total</div><div class="value">{n_total}</div></div>
  <div class="kpi"><div class="label">Done</div><div class="value">{n_done}</div></div>
  <div class="kpi"><div class="label">Pending</div><div class="value">{n_pending}</div></div>
  <div class="kpi"><div class="label">Best \u0394G_bind</div><div class="value">{kpi_best}</div></div>
  <div class="kpi"><div class="label">Mean \u0394G_bind</div><div class="value">{kpi_mean}</div></div>
</div>

{screen_link}

<h2>Ranking</h2>
{chart_html}

<h2>Ligands</h2>
<p class="meta">Click any column header to sort. Click a ligand name to open its full per-ligand report.</p>
<table class="sortable">
  <thead><tr>
    <th class="numeric">Rank</th>
    <th>Ligand</th>
    <th class="numeric">\u0394G_bind (kcal/mol)</th>
    <th class="numeric">Std (kcal/mol)</th>
    <th class="numeric">SEM (kcal/mol)</th>
    <th>Status</th>
  </tr></thead>
  <tbody>
  {ranked_joined}
  {pending_joined}
  </tbody>
</table>

<div class="footer">
  Generated by <code>amber_md.batch_aggregate</code> v2.5.1.
  Per-ligand reports are at <code>lig_*/mmgbsa/FINAL_RESULTS.report.html</code>.
  \u0394G_bind = single-trajectory MM/GBSA, no entropy correction.
  Use for <b>relative ranking</b>, not absolute Kd estimates.
</div>

</body></html>
"""


# ----------------------------- driver -----------------------------

def main():
    if len(sys.argv) != 2:
        print(__doc__); sys.exit(1)
    batch_dir = Path(sys.argv[1]).resolve()
    if not batch_dir.is_dir():
        print(f"ERROR: not a directory: {batch_dir}"); sys.exit(1)

    _run_screen_analysis_kit(batch_dir)

    rows = []
    for wd in sorted(batch_dir.glob("lig_*")):
        if not wd.is_dir(): continue
        # final49: tolerate MMPBSA (pipeline) or MMGBSA (older/manual) name.
        dat = wd / "mmgbsa" / "FINAL_RESULTS_MMPBSA.dat"
        if not dat.exists():
            _alt = wd / "mmgbsa" / "FINAL_RESULTS_MMGBSA.dat"
            if _alt.exists():
                dat = _alt
        result = parse_mmpbsa_dat(dat)
        if result:
            status = "DONE"
        elif (wd / "jobs" / "prod.nc").exists():
            status = "MD done, MMGBSA pending"
        elif (wd / "build" / "complex.prmtop").exists():
            status = "BUILD done, MD pending"
        else:
            status = "Not started / failed"
        rows.append({"workdir": wd.name, "status": status,
                     "result": result, "path": dat})

    report_module = _load_mmpbsa_report_module()
    n_reports_built = 0
    n_reports_existing = 0
    for r in rows:
        if not r["result"]: continue
        existing = r["path"].with_name("FINAL_RESULTS.report.html")
        if existing.exists() and existing.stat().st_size > 0:
            r["report_path"] = existing
            n_reports_existing += 1
            continue
        rp = _ensure_per_ligand_report(r["path"], report_module)
        if rp:
            r["report_path"] = rp
            n_reports_built += 1
        else:
            r["report_path"] = None

    tsv = batch_dir / "binding_energies.tsv"
    with open(tsv, "w") as f:
        f.write("ligand\tstatus\tdelta_total_kcal_mol\tstd\tsem\n")
        for r in rows:
            if r["result"]:
                f.write(f"{r['workdir']}\t{r['status']}\t"
                        f"{r['result']['delta_total']:.3f}\t"
                        f"{r['result'].get('std') or ''}\t"
                        f"{r['result'].get('sem') or ''}\n")
            else:
                f.write(f"{r['workdir']}\t{r['status']}\t\t\t\n")

    done = [r for r in rows if r["result"]]
    done.sort(key=lambda r: r["result"]["delta_total"])
    md = batch_dir / "binding_energies_ranked.md"
    with open(md, "w") as f:
        f.write(f"# Binding free energies -- {batch_dir.name}\n\n")
        f.write(f"Aggregated from {len(done)}/{len(rows)} ligands (DONE).\n\n")
        f.write("| Rank | Ligand | DG_bind (kcal/mol) | std | SEM |\n")
        f.write("|------|--------|-------------------:|----:|----:|\n")
        for rank, r in enumerate(done, 1):
            res = r["result"]
            std = f"{res['std']:.2f}" if res.get("std") else "--"
            sem = f"{res['sem']:.2f}" if res.get("sem") else "--"
            f.write(f"| {rank} | `{r['workdir']}` | {res['delta_total']:+.2f} | {std} | {sem} |\n")
        pending = [r for r in rows if not r["result"]]
        if pending:
            f.write(f"\n## Still pending ({len(pending)})\n\n")
            for r in pending:
                f.write(f"- `{r['workdir']}` ({r['status']})\n")

    index_html = batch_dir / "INDEX.html"
    try:
        index_html.write_text(_build_index_html(batch_dir, rows))
        index_written = True
    except Exception as e:
        sys.stderr.write(f"[batch_aggregate] INDEX.html failed: {e}\n")
        index_written = False

    print(f"Wrote {tsv}")
    print(f"Wrote {md}")
    if index_written:
        print(f"Wrote {index_html}")
    if report_module is None:
        print("[note] mmpbsa_report not found; skipped per-ligand HTML reports.")
    else:
        print(f"Per-ligand HTML reports: {n_reports_existing} already present, "
              f"{n_reports_built} newly generated.")
    print(f"Done: {len(done)}/{len(rows)} ligands have MMGBSA results.")
    if done:
        print("\nTop 5 hits:")
        for rank, r in enumerate(done[:5], 1):
            print(f"  {rank}. {r['workdir']:50s}  "
                  f"{r['result']['delta_total']:+8.2f} kcal/mol")


if __name__ == "__main__":
    main()
