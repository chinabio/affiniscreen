#!/usr/bin/env python3
"""build_screen_summary.py
Aggregate every analysis/summary.json under a screen directory into
ONE sortable master HTML.

Usage:  build_screen_summary.py <screen_dir>
Writes: <screen_dir>/screen_summary.html
"""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
import json, sys
from pathlib import Path


def composite_score(s):
    """Lower-is-better composite score:
       composite = dG + 0.5 * mean_lig_RMSD + max(0, max_bb_RMSD - 3.0)
       Missing pieces are skipped/penalized as Inf."""
    if s.get("dG") is None:
        return float("inf")
    score = s["dG"]
    if s.get("rmsd_lig_mean") is not None:
        score += 0.5 * s["rmsd_lig_mean"]
    if s.get("rmsd_bb_max") is not None and s["rmsd_bb_max"] > 3.0:
        score += s["rmsd_bb_max"] - 3.0
    return score


def cls_dG(v):
    if v is None: return "neutral"
    if v <= -40:  return "good"
    if v <= -25:  return "warn"
    return "bad"


def cls_rmsd(v, soft=2.0, hard=4.0):
    if v is None: return "neutral"
    if v < soft:  return "good"
    if v < hard:  return "warn"
    return "bad"


def fmt(v, spec="{:+.2f}"):
    return "-" if v is None else spec.format(v)


def main(screen_dir: Path):
    rows = []
    for sj in sorted(screen_dir.glob("lig_*/analysis/summary.json")):
        try:
            rows.append(json.loads(sj.read_text()))
        except Exception as e:
            print(f"  skip {sj}: {e}", file=sys.stderr)

    if not rows:
        sys.exit(f"No analysis/summary.json found under {screen_dir}")

    rows.sort(key=composite_score)

    body = []
    for rank, s in enumerate(rows, 1):
        link = s.get("report_path", "")
        body.append(f"""<tr>
  <td class='num'>{rank}</td>
  <td><a href="{link}"><code>{s['ligand']}</code></a></td>
  <td class='num {cls_dG(s.get("dG"))}'>{fmt(s.get("dG"))}</td>
  <td class='num'>{fmt(s.get("sem"), "{:.2f}")}</td>
  <td class='num'>{fmt(s.get("VDWAALS"))}</td>
  <td class='num'>{fmt(s.get("EEL"))}</td>
  <td class='num'>{fmt(s.get("EGB"))}</td>
  <td class='num'>{fmt(s.get("ESURF"))}</td>
  <td class='num {cls_rmsd(s.get("rmsd_bb_mean"))}'>{fmt(s.get("rmsd_bb_mean"), "{:.2f}")}</td>
  <td class='num {cls_rmsd(s.get("rmsd_lig_mean"), 1.5, 3.0)}'>{fmt(s.get("rmsd_lig_mean"), "{:.2f}")}</td>
  <td class='num'>{s.get("n_hbonds_strong", 0)}</td>
  <td class='num'>{fmt(composite_score(s), "{:+.2f}")}</td>
</tr>""")

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Screen summary - {screen_dir.name}</title>
<style>
 body{{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:1300px;margin:24px auto;padding:0 20px;color:#222}}
 h1{{border-bottom:3px solid #1971c2;padding-bottom:6px}}
 .meta{{color:#666;font-size:.9em}}
 table{{border-collapse:collapse;width:100%;font-size:.92em;margin-top:14px}}
 th,td{{padding:6px 9px;border-bottom:1px solid #e6e6e6;text-align:left;vertical-align:top}}
 th{{background:#1971c2;color:white;cursor:pointer;user-select:none;position:sticky;top:0}}
 th:hover{{background:#1864ab}} th::after{{content:" \u21c5";opacity:.5}}
 tr:nth-child(even) td{{background:#fafbfc}}
 td.num{{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}}
 .good{{color:#2b8a3e;font-weight:600}}
 .warn{{color:#b08900;font-weight:600}}
 .bad {{color:#c92a2a;font-weight:600}}
 .legend{{font-size:.85em;color:#555;margin:10px 0 0}}
 code{{background:#f1f3f5;padding:1px 5px;border-radius:3px}}
 input#filter{{padding:6px 10px;border:1px solid #ccc;border-radius:5px;width:280px;margin-top:10px}}
</style></head><body>

<h1>Screen summary - <code>{screen_dir.name}</code></h1>
<p class="meta">{len(rows)} ligand(s) ranked by composite score (lower = better).
Click any header to sort. Type to filter.</p>

<input id="filter" placeholder="filter by ligand name...">

<table id="t">
 <thead><tr>
   <th data-type="num">#</th>
   <th>Ligand</th>
   <th data-type="num">&Delta;G (kcal/mol)</th>
   <th data-type="num">&plusmn;SEM</th>
   <th data-type="num">VDW</th>
   <th data-type="num">EEL</th>
   <th data-type="num">EGB</th>
   <th data-type="num">ESURF</th>
   <th data-type="num">&lang;BB RMSD&rang; (A)</th>
   <th data-type="num">&lang;Lig RMSD&rang; (A)</th>
   <th data-type="num">strong H-bonds</th>
   <th data-type="num">composite</th>
 </tr></thead>
 <tbody>
 {''.join(body)}
 </tbody>
</table>

<p class="legend">
 <span class="good">green</span> = good &middot; <span class="warn">amber</span> = caution &middot;
 <span class="bad">red</span> = poor &nbsp;|&nbsp;
 Composite = &Delta;G + 0.5&middot;&lang;Lig RMSD&rang; + max(0, BB RMSD<sub>max</sub> - 3 A).
 Strong H-bond = occupancy &ge; 0.5.
</p>

<script>
document.querySelectorAll('#t th').forEach((th, i) => {{
  let asc = true;
  th.addEventListener('click', () => {{
    const tbody = document.querySelector('#t tbody');
    const rows  = [...tbody.querySelectorAll('tr')];
    const num   = th.dataset.type === 'num';
    rows.sort((a,b) => {{
      let x = a.children[i].innerText.trim(), y = b.children[i].innerText.trim();
      if (num) {{ x = parseFloat(x); y = parseFloat(y);
                 if (isNaN(x)) x = Infinity; if (isNaN(y)) y = Infinity; }}
      return (x>y?1:x<y?-1:0) * (asc?1:-1);
    }});
    asc = !asc;
    rows.forEach(r => tbody.appendChild(r));
  }});
}});
document.getElementById('filter').addEventListener('input', e => {{
  const q = e.target.value.toLowerCase();
  document.querySelectorAll('#t tbody tr').forEach(r => {{
    r.style.display = r.children[1].innerText.toLowerCase().includes(q) ? '' : 'none';
  }});
}});
</script>
</body></html>"""

    out = screen_dir / "screen_summary.html"
    out.write_text(html)
    print(f"[screen_summary] wrote {out}  ({len(rows)} ligands)")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: build_screen_summary.py <screen_dir>")
    main(Path(sys.argv[1]).resolve())
