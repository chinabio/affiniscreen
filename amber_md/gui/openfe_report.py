"""amber_md.gui.openfe_report -- reusable OpenFE FEP/ABFE report generator.

Point :func:`build_campaign_report` at a results directory; it auto-discovers
``*_result.json`` files, groups them by ligand/transformation, treats every
matching file as one *replicate* (the replicate COUNT is whatever is on disk --
1, 2, 3, 5, ... -- never hard-coded), decodes the embedded MBAR diagnostics,
and emits one self-contained HTML per ligand plus a campaign index.html.
"""
from __future__ import annotations
import base64
import math, bz2, io, json, re, statistics, zipfile
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path


def _mag(x):
    return x.get("magnitude") if isinstance(x, dict) else x


def _decode_ndarray(node):
    if not isinstance(node, dict) or node.get("__class__") != "ndarray":
        return None
    import numpy as np
    try:
        dtype = node["dtype"]; shape = tuple(node["shape"]); b = node["bytes"]
        raw = b["latin-1"].encode("latin-1") if "latin-1" in b else base64.b64decode(b.get("base64", ""))
        try:
            import zstandard as zstd; raw = zstd.ZstdDecompressor().decompress(raw)
        except Exception:
            try: raw = bz2.decompress(raw)
            except Exception: pass
        return np.frombuffer(raw, dtype=dtype).reshape(shape)
    except Exception:
        return None


def _load_json(path):
    path = Path(path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path) as zf:
            jn = [n for n in zf.namelist() if n.endswith(".json")]
            if not jn: raise ValueError(f"no .json inside {path}")
            with zf.open(jn[0]) as fh: return json.load(fh)
    with open(path) as fh: return json.load(fh)


_REP_PATTERNS = [
    re.compile(r"^(?P<key>.+?)[_-]rep(?P<rep>\d+)$", re.I),
    re.compile(r"^(?P<key>.+?)[_-]r(?P<rep>\d+)$", re.I),
    re.compile(r"^(?P<key>.+?)[_-]replicate(?P<rep>\d+)$", re.I),
    re.compile(r"^(?P<key>.+?)[_-]repeat[_-]?(?P<rep>\d+)$", re.I),
]


def parse_ligand_replicate(stem):
    stem = re.sub(r"_result$", "", stem)
    stem = re.sub(r"^results?_", "", stem)
    for pat in _REP_PATTERNS:
        m = pat.match(stem)
        if m: return m.group("key"), f"rep{int(m.group('rep'))}"
    return stem, "rep0"


@dataclass
class LegResult:
    dG: object = None
    dGerr: object = None
    ssc: object = None
    prod: object = None
    equil: object = None
    overlap: object = None
    fr_frac: object = None
    forward_DGs: object = None
    reverse_DGs: object = None
    forward_dDGs: object = None
    reverse_dDGs: object = None
    min_adj_overlap: object = None


@dataclass
class ReplicateResult:
    label: str
    path: str
    estimate: object = None
    uncertainty: object = None
    n_ok: int = 0
    n_fail: int = 0
    legs: dict = field(default_factory=dict)
    error: str = ""


def parse_replicate(path, label):
    rr = ReplicateResult(label=label, path=str(path))
    try:
        d = _load_json(path)
    except Exception as e:
        rr.error = f"unreadable: {e}"; return rr
    rr.estimate = _mag(d.get("estimate")); rr.uncertainty = _mag(d.get("uncertainty"))
    ur = d.get("unit_results", {})
    items = list(ur.values()) if isinstance(ur, dict) else list(ur or [])
    for u in items:
        if not isinstance(u, dict): continue
        nm = u.get("name", "")
        if u.get("exception") or "Failure" in u.get("__qualname__", ""):
            rr.n_fail += 1; continue
        rr.n_ok += 1
        if "Analysis" not in nm: continue
        leg = "complex" if "complex" in nm else ("solvent" if "solvent" in nm else nm[:30])
        o = u.get("outputs", {})
        lr = LegResult(dG=_mag(o.get("unit_estimate")), dGerr=_mag(o.get("unit_estimate_error")),
                       ssc=_mag(o.get("standard_state_correction")),
                       prod=o.get("production_iterations"), equil=o.get("equilibration_iterations"))
        try: lr.overlap = _decode_ndarray((o.get("unit_mbar_overlap") or {}).get("matrix"))
        except Exception: lr.overlap = None
        if lr.overlap is not None and lr.overlap.shape[0] > 1:
            M = lr.overlap
            lr.min_adj_overlap = float(min(M[i, i+1] for i in range(M.shape[0]-1)))
        fr = o.get("forward_and_reverse_energies", {}) or {}
        lr.fr_frac = _decode_ndarray(fr.get("fractions"))
        for fk in ("forward_DGs","reverse_DGs","forward_dDGs","reverse_dDGs"):
            nd = fr.get(fk, {}); m = nd.get("magnitude") if isinstance(nd, dict) else None
            if isinstance(m, dict): setattr(lr, fk, _decode_ndarray(m))
            elif isinstance(m, list):
                import numpy as np; setattr(lr, fk, np.array(m))
        rr.legs[leg] = lr
    return rr


def discover_ligands(results_dir):
    """Return {ligand_key: [ReplicateResult,...]} under results_dir.

    Handles BOTH common OpenFE layouts:
      (a) flat files     : <root>/<ligand>_rep<N>_result.json
      (b) per-run dirs   : <root>/<ligand>_rep<N>/results.json
                           (replicate taken from the DIRECTORY name)
    Any number of replicates is supported. quickrun_cache is ignored.
    """
    root = Path(results_dir).expanduser()
    groups = {}
    if not root.exists():
        return groups

    seen_paths = set()

    # --- (a) *_result.json anywhere (replicate parsed from file stem) ---
    for rj in sorted(root.rglob("*_result.json")):
        if "quickrun_cache" in rj.parts:
            continue
        seen_paths.add(rj.resolve())
        key, label = parse_ligand_replicate(rj.stem)
        groups.setdefault(key, []).append(parse_replicate(rj, label))

    # --- (b) bare results.json inside a per-run directory ---
    #     ligand key + replicate come from the PARENT directory name.
    for rj in sorted(list(root.rglob("results.json")) + list(root.rglob("result.json"))):
        if "quickrun_cache" in rj.parts:
            continue
        if rj.resolve() in seen_paths:
            continue
        seen_paths.add(rj.resolve())
        key, label = parse_ligand_replicate(rj.parent.name)
        groups.setdefault(key, []).append(parse_replicate(rj, label))

    # de-duplicate replicate labels within a ligand
    for key, reps in groups.items():
        seen = {}
        for r in reps:
            if r.label in seen:
                seen[r.label] += 1
                r.label = f"{r.label}.{seen[r.label]}"
            else:
                seen[r.label] = 0
        reps.sort(key=lambda r: r.label)
    return groups


def replicate_stats(reps):
    ests = [r.estimate for r in reps if isinstance(r.estimate, (int, float))]
    n = len(ests); out = {"n": n, "values": ests, "mean": None, "sd": None, "sem": None}
    if n >= 1: out["mean"] = statistics.mean(ests)
    if n >= 2: out["sd"] = statistics.stdev(ests); out["sem"] = out["sd"]/(n**0.5)
    return out


def _fig_to_b64(fig):
    import matplotlib.pyplot as plt
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=110, bbox_inches="tight"); plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()

def _esc(x):
    import html as _h; return _h.escape(str(x))

def _fnum(x, n=2):
    return f"{x:.{n}f}" if isinstance(x, (int, float)) else "n/a"

def _fig_estimate(reps, st, reference=None):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    labels=[r.label for r in reps]
    vals=[r.estimate if isinstance(r.estimate,(int,float)) else float("nan") for r in reps]
    fig,ax=plt.subplots(figsize=(max(5,1.4*len(reps)),3.6))
    ax.bar(labels,vals,color="#4c78a8",width=0.55,label="replicate")
    if st["mean"] is not None:
        ax.axhline(st["mean"],color="#e45756",ls="--",lw=2,label=f"mean {st['mean']:.2f}")
        if st["sem"] is not None:
            ax.axhspan(st["mean"]-st["sem"],st["mean"]+st["sem"],color="#e45756",alpha=0.12,label=f"+/-SEM {st['sem']:.2f}")
    if reference is not None:
        ax.axhline(reference[1],color="#54a24b",ls=":",lw=2.2,label=f"{reference[0]} {reference[1]:.2f}")
    ax.set_ylabel("dG (kcal/mol)"); ax.set_title("Binding free energy by replicate")
    ax.legend(fontsize=8,loc="best"); ax.grid(axis="y",alpha=0.3)
    return _fig_to_b64(fig)

def _fig_legs(reps):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    legnames=sorted({lg for r in reps for lg in r.legs})
    if not legnames: return None
    fig,axs=plt.subplots(1,len(legnames),figsize=(4.2*len(legnames),3.3),squeeze=False)
    labels=[r.label for r in reps]
    for ax,leg in zip(axs[0],legnames):
        ys=[r.legs[leg].dG if leg in r.legs else None for r in reps]
        es=[r.legs[leg].dGerr if leg in r.legs else None for r in reps]
        ax.errorbar(labels,[y if y is not None else float("nan") for y in ys],
                    yerr=[e if e is not None else 0 for e in es],fmt="o-",capsize=4,color="#4c78a8")
        ax.set_title(f"{leg} leg dG"); ax.set_ylabel("kcal/mol"); ax.grid(alpha=0.3)
    fig.suptitle("Per-leg decoupling dG (+/- MBAR error)")
    return _fig_to_b64(fig)

def _fig_overlap(reps, leg="complex"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    have=[r for r in reps if leg in r.legs and r.legs[leg].overlap is not None]
    if not have: return None
    fig,axs=plt.subplots(1,len(have),figsize=(3.6*len(have),3.4),squeeze=False)
    for ax,r in zip(axs[0],have):
        M=r.legs[leg].overlap
        im=ax.imshow(M,cmap="viridis",vmin=0,vmax=max(0.5,float(M.max())))
        ax.set_title(f"{r.label}\nmin adj={_fnum(r.legs[leg].min_adj_overlap,3)}")
        ax.set_xlabel("lambda"); ax.set_ylabel("lambda"); fig.colorbar(im,ax=ax,fraction=0.046)
    fig.suptitle(f"MBAR lambda-overlap ({leg} leg)")
    return _fig_to_b64(fig)

def _fig_convergence(reps, leg="complex"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    have=[r for r in reps if leg in r.legs and r.legs[leg].fr_frac is not None and r.legs[leg].forward_DGs is not None]
    if not have: return None
    fig,axs=plt.subplots(1,len(have),figsize=(3.6*len(have),3.3),squeeze=False)
    for ax,r in zip(axs[0],have):
        e=r.legs[leg]
        ax.plot(e.fr_frac,e.forward_DGs,"o-",color="#4c78a8",ms=3,label="fwd")
        if e.forward_dDGs is not None: ax.fill_between(e.fr_frac,e.forward_DGs-e.forward_dDGs,e.forward_DGs+e.forward_dDGs,color="#4c78a8",alpha=.2)
        if e.reverse_DGs is not None:
            ax.plot(e.fr_frac,e.reverse_DGs,"s--",color="#e45756",ms=3,label="rev")
            if e.reverse_dDGs is not None: ax.fill_between(e.fr_frac,e.reverse_DGs-e.reverse_dDGs,e.reverse_DGs+e.reverse_dDGs,color="#e45756",alpha=.2)
        ax.set_title(r.label); ax.set_xlabel("sim fraction"); ax.set_ylabel("dG"); ax.legend(fontsize=7); ax.grid(alpha=0.3)
    fig.suptitle(f"Forward/reverse convergence ({leg} leg)")
    return _fig_to_b64(fig)

_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#1a1a1a;max-width:980px;margin:24px auto;padding:0 20px;line-height:1.45}
h1{font-size:22px;border-bottom:3px solid #4c78a8;padding-bottom:6px}
h2{font-size:17px;color:#2c3e50;margin-top:28px;border-left:4px solid #4c78a8;padding-left:8px}
table{border-collapse:collapse;width:100%;margin:12px 0;font-size:13px}
th,td{border:1px solid #ccc;padding:6px 9px;text-align:center}
th{background:#4c78a8;color:#fff} tr:nth-child(even){background:#f4f7fb}
.kpi{display:flex;gap:14px;flex-wrap:wrap;margin:14px 0}
.card{flex:1;min-width:140px;background:#f4f7fb;border:1px solid #d6e0ec;border-radius:8px;padding:12px;text-align:center}
.card .v{font-size:23px;font-weight:700;color:#2c3e50}.card .l{font-size:11px;color:#666}
.good{color:#2a8a4a;font-weight:600}.warn{color:#c47f17;font-weight:600}.bad{color:#c0392b;font-weight:600}
img{max-width:100%;border:1px solid #e3e3e3;border-radius:6px;margin:6px 0}
.note{background:#fff8e6;border-left:4px solid #e0a91b;padding:8px 12px;font-size:13px;border-radius:4px}
.foot{margin-top:28px;font-size:11px;color:#888;border-top:1px solid #ddd;padding-top:8px} a{color:#2b6cb0}
"""



# --------------------------------------------------------------------------- #
# OpenFE-produced diagnostic PNGs (embedded, not regenerated)
# --------------------------------------------------------------------------- #
_PNG_NAMES = ("forward_reverse_convergence", "mbar_overlap_matrix",
              "replica_exchange_matrix", "replica_state_timeseries")


def find_openfe_pngs(result_path):
    """Return {leg: {plotname: Path}} for PNGs OpenFE wrote next to results.json.

    Looks for shared_ABFE{Complex,Solvent}AnalysisUnit-*/<plot>.png in the
    directory containing the result file (the run folder is assumed intact).
    """
    base = Path(result_path).parent
    found = {}
    for unit in sorted(base.glob("shared_ABFE*AnalysisUnit-*")):
        nm = unit.name.lower()
        leg = "complex" if "complex" in nm else ("solvent" if "solvent" in nm else "other")
        for png in _PNG_NAMES:
            p = unit / f"{png}.png"
            if p.exists():
                found.setdefault(leg, {})[png] = p
    return found


def _png_to_b64(path):
    try:
        return "data:image/png;base64," + base64.b64encode(Path(path).read_bytes()).decode()
    except Exception:
        return None


def render_ligand_html(ligand_key, reps, reference=None):
    st=replicate_stats(reps); n=st["n"]
    legnames=sorted({lg for r in reps for lg in r.legs})
    if n>=2:
        unc_txt=f"{st['mean']:.2f} +/- {st['sem']:.2f} kcal/mol (mean +/- SEM, n={n})"
    elif n==1:
        err=None
        for lg in reps[0].legs.values():
            if lg.dGerr is not None: err=lg.dGerr
        unc_txt=f"{st['mean']:.2f} kcal/mol (n=1 -- no replicate spread; single-run MBAR error ~ {_fnum(err)})"
    else:
        unc_txt="no numeric estimate found"
    fig1=_fig_estimate(reps,st,reference) if n>=1 else None
    fig2=_fig_legs(reps)
    fig3=_fig_overlap(reps,"complex") or (_fig_overlap(reps,legnames[0]) if legnames else None)
    fig4=_fig_convergence(reps,"complex")
    rows=""
    for r in reps:
        cl=r.legs.get("complex"); sv=r.legs.get("solvent")
        status=(f"<span class='good'>{r.n_ok}/{r.n_ok+r.n_fail} ok</span>" if r.n_fail==0 and not r.error
                else f"<span class='bad'>{r.error or (str(r.n_fail)+' of '+str(r.n_ok+r.n_fail)+' failed')}</span>")
        # simulation length from the complex leg (falls back to solvent) -- real
        # info from the JSON, replacing the previously-confusing bare unit count.
        _leg = cl or sv
        if _leg and (_leg.prod is not None or _leg.equil is not None):
            simlen = f"{int(_leg.prod) if _leg.prod else '-'} / {int(_leg.equil) if _leg.equil else '-'}"
        else:
            simlen = "-"
        rows+=(f"<tr><td><b>{_esc(r.label)}</b></td><td>{_fnum(r.estimate)}</td>"
               f"<td>{(_fnum(cl.dG,3)+' +/- '+_fnum(cl.dGerr,2)) if cl else '-'}</td>"
               f"<td>{(_fnum(sv.dG,3)+' +/- '+_fnum(sv.dGerr,2)) if sv else '-'}</td>"
               f"<td>{simlen}</td><td>{status}</td></tr>")
    summary_row=""
    if n>=2:
        summary_row=(f"<tr style='background:#eaf3ea;font-weight:700'><td>Mean +/- SEM</td>"
                     f"<td colspan='5'>{st['mean']:.2f} +/- {st['sem']:.2f} kcal/mol (SD {st['sd']:.2f}, n={n})</td></tr>")
    ref_html=""
    # Cross-check section is OMITTED entirely unless an external reference
    # (e.g. FEP+) was explicitly supplied for THIS ligand. No reference is the
    # normal case, so by default no cross-check appears anywhere in the report.
    if reference and st["mean"] is not None:
        diff=st["mean"]-reference[1]; verdict="good" if abs(diff)<=1.5 else "warn"
        ref_html=(f"<h2>Cross-check vs {_esc(reference[0])}</h2><p>OpenFE mean <b>{st['mean']:.2f}</b> vs "
                  f"{_esc(reference[0])} <b>{reference[1]:.2f}</b> -> difference <b class='{verdict}'>{diff:+.2f} kcal/mol</b>. "
                  f"<span class='note'>Only the final dG is comparable across engines; per-leg conventions differ.</span></p>")
    overlaps=[r.legs["complex"].min_adj_overlap for r in reps if "complex" in r.legs and r.legs["complex"].min_adj_overlap is not None]
    min_ov=min(overlaps) if overlaps else None
    ov_verdict=("All >= 0.03 -- lambda-spacing adequate." if (min_ov or 0)>=0.03 else "Some adjacent overlaps < 0.03 -- consider more lambda-windows.")
    # ---- Embed OpenFE's own diagnostic PNGs (per replicate, per leg) ----
    png_sections = []
    # gather: {plotname: [(replabel, leg, b64), ...]}
    by_plot = {}
    for r in reps:
        pngs = find_openfe_pngs(r.path)
        for leg, plots in pngs.items():
            for plotname, p in plots.items():
                b64 = _png_to_b64(p)
                if b64:
                    by_plot.setdefault(plotname, []).append((r.label, leg, b64))

    _PLOT_TITLES = {
        "forward_reverse_convergence": "Forward / Reverse Convergence",
        "mbar_overlap_matrix": "MBAR Overlap Matrix (OpenFE)",
        "replica_exchange_matrix": "Replica-Exchange Transition Matrix",
        "replica_state_timeseries": "Replica State Time-series",
    }
    for plotname in ("forward_reverse_convergence", "mbar_overlap_matrix",
                     "replica_exchange_matrix", "replica_state_timeseries"):
        entries = by_plot.get(plotname)
        if not entries:
            continue
        title = _PLOT_TITLES[plotname]
        legs_present = sorted({lg for _, lg, _ in entries})
        caveat = ""
        if plotname == "forward_reverse_convergence" and legs_present == ["solvent"]:
            caveat = ("<p class='note'>OpenFE produced forward/reverse "
                      "convergence for the <b>solvent leg only</b> (the complex "
                      "leg plot was not generated). This reflects ligand "
                      "decoupling in water, not the complex/binding leg.</p>")
        imgs = "".join(
            f"<div style='display:inline-block;margin:6px;text-align:center'>"
            f"<div style='font-size:12px;color:#555'>{_esc(lbl)} · {leg} leg</div>"
            f"<img src='{b64}' style='max-width:340px'></div>"
            for lbl, leg, b64 in entries)
        png_sections.append(f"<h2>{title}</h2>{caveat}{imgs}")

    conv_html = "\n".join(png_sections) if png_sections else (
        "<h2>Convergence diagnostics</h2><p class='note'>No OpenFE diagnostic "
        "PNGs (forward/reverse, overlap, replica-exchange) were found next to "
        "the result files, and the forward/reverse arrays in the JSON are null, "
        "so these plots are omitted rather than fabricated.</p>")
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>OpenFE report -- {_esc(ligand_key)}</title><style>{_CSS}</style></head><body>
<h1>OpenFE FEP Report -- {_esc(ligand_key)}</h1>
<p><b>Replicates discovered:</b> {n} &nbsp;|&nbsp; <b>Generated:</b> {datetime.now():%Y-%m-%d %H:%M}</p>
<div class="kpi">
 <div class="card"><div class="v">{_fnum(st['mean'])}</div><div class="l">Mean dG (kcal/mol)</div></div>
 <div class="card"><div class="v">{('+/-'+_fnum(st['sem'])) if st['sem'] is not None else 'n/a'}</div><div class="l">SEM</div></div>
 <div class="card"><div class="v">{_fnum(st['sd']) if st['sd'] is not None else 'n/a'}</div><div class="l">Std dev</div></div>
 <div class="card"><div class="v">{n}</div><div class="l">Replicates</div></div></div>
<p><b>Result:</b> {unc_txt}</p>
<h2>Replicate Results</h2>
<table><tr><th>Replicate</th><th>dG (kcal/mol)</th><th>Complex leg</th><th>Solvent leg</th><th>Sim length<br><span style='font-weight:400;font-size:11px'>prod/equil iters</span></th><th>Protocol units</th></tr>
{rows}{summary_row}</table>
{('<h2>dG by replicate</h2><img src="'+fig1+'">') if fig1 else ''}
{ref_html}
{('<h2>Per-leg free energies</h2><img src="'+fig2+'">') if fig2 else ''}
{('<h2>MBAR lambda-overlap</h2><img src="'+fig3+'"><p>Minimum adjacent overlap: <b>'+_fnum(min_ov,3)+'</b>. '+ov_verdict+'</p>') if fig3 else ''}
{conv_html}
<div class="foot">Auto-generated by amber_md.gui.openfe_report. Replicate count is whatever was found on disk ({n}).
Headline uncertainty is the between-replicate SEM when n>=2, else the single-run MBAR error.</div></body></html>"""

def render_campaign_index(report_rows):
    head=(f"<!doctype html><html><head><meta charset='utf-8'><title>OpenFE campaign summary</title><style>{_CSS}</style></head><body>"
          f"<h1>OpenFE Campaign Summary</h1><p>{len(report_rows)} ligand(s)/transformation(s)</p>"
          f"<table><tr><th>Ligand / transformation</th><th>n reps</th><th>Mean dG</th><th>SEM</th><th>SD</th><th>Status</th><th>Report</th></tr>")
    body=""
    for row in report_rows:
        st=row["stats"]; status=row["status"]; cls={"ok":"good","partial":"warn","failed":"bad"}.get(status,"")
        body+=(f"<tr><td style='text-align:left'>{_esc(row['key'])}</td><td>{st['n']}</td><td>{_fnum(st['mean'])}</td>"
               f"<td>{('+/-'+_fnum(st['sem'])) if st['sem'] is not None else '-'}</td><td>{_fnum(st['sd']) if st['sd'] is not None else '-'}</td>"
               f"<td class='{cls}'>{status}</td><td><a href='{_esc(row['html_name'])}'>open</a></td></tr>")
    return head+body+"</table></body></html>"

def build_campaign_report(results_dir, out_dir=None, reference=None, experimental=None):
    """Discover replicates per ligand and emit per-ligand + campaign HTML.

    reference is OPTIONAL and defaults to None. When omitted (the normal case --
    there is usually no FEP+ data) NO cross-check section is rendered. Pass it
    only for the rare ligand that has an external value, e.g.
    reference={"12944901": ("FEP+", -7.56)}.
    """
    results_dir=Path(results_dir).expanduser()
    out_dir=Path(out_dir).expanduser() if out_dir else results_dir/"_report"
    out_dir.mkdir(parents=True, exist_ok=True)
    reference=reference or {}

    # --- RBFE auto-detection: if this directory holds rbfe_* edge files,
    #     produce the relative-binding (ddG = complex - solvent) report instead
    #     of the per-ligand ABFE layout. ---
    if is_rbfe_dir(results_dir):
        # Experiment comparison is OPTIONAL and OFF by default. The report's
        # core output (per-edge ddG = complex - solvent) NEVER depends on it.
        # It is used only when the caller explicitly passes `experimental`
        # (a dict {ligand: dG_kcal_mol} or a path to a ligands.yml). We do NOT
        # auto-search the filesystem for any experimental file.
        exp = experimental
        if isinstance(exp, (str, Path)):
            exp = parse_experimental_yaml(exp)   # explicit path only
        rres = render_rbfe_report(results_dir, out_dir, experimental=exp or {})
        return {"index": rres["index"], "ligands": {}, "rows": rres["rows"],
                "mode": "rbfe", "stats": rres.get("stats", {})}

    groups=discover_ligands(results_dir)
    out={"index":None,"ligands":{},"rows":[]}; rows=[]
    for key,reps in sorted(groups.items()):
        st=replicate_stats(reps)
        any_fail=any(r.n_fail>0 or r.error for r in reps)
        status="failed" if st["mean"] is None else ("partial" if any_fail else "ok")
        ref=None
        for rk,rv in reference.items():
            if rk in key: ref=rv; break
        safe=re.sub(r"[^A-Za-z0-9._-]","_",key); html_name=f"ligand_{safe}.html"
        (out_dir/html_name).write_text(render_ligand_html(key,reps,ref))
        out["ligands"][key]=str(out_dir/html_name)
        rows.append({"key":key,"stats":st,"status":status,"html_name":html_name})
    out["rows"]=rows
    index=out_dir/"index.html"; index.write_text(render_campaign_index(rows)); out["index"]=str(index)
    return out


# =========================================================================== #
# RBFE (relative) campaign support
#   Edges are named:
#     rbfe_<ligA>_complex_<ligB>_complex_rep<N>_result.json   (complex leg)
#     rbfe_<ligA>_solvent_<ligB>_solvent_rep<N>_result.json   (solvent leg)
#   Predicted ddG(A->B) = dG_complex - dG_solvent  (per replicate, then averaged)
# =========================================================================== #
import re as _re_rbfe

_EDGE_RE = _re_rbfe.compile(
    r"rbfe_(lig_.+?)_(complex|solvent)_(lig_.+?)_(complex|solvent)_rep(\d+)_result")

_RT_KCAL = 0.001987204258 * 298.15   # kcal/mol at 298.15 K


def is_rbfe_dir(results_dir):
    """True if the directory contains rbfe_* edge result files."""
    root = Path(results_dir).expanduser()
    if not root.exists():
        return False
    for rj in root.rglob("rbfe_*_result.json"):
        if "quickrun_cache" not in rj.parts:
            return True
    return False


def parse_experimental_yaml(yaml_path):
    """Parse a ligands.yml -> {ligand_name: dG_exp_kcal_mol}.

    OPTIONAL helper. Returns {} (never raises) if pyyaml is unavailable, the
    path is missing/unreadable, or no usable measurements are present. The
    report does not require this; it only enriches the RBFE table when supplied.
    """
    try:
        import yaml  # noqa
    except Exception:
        return {}
    try:
        p = Path(yaml_path)
        if not p.exists():
            return {}
        doc = yaml.safe_load(p.read_text())
    except Exception:
        return {}
    scale = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9, "pM": 1e-12}
    out = {}
    for lig, d in (doc or {}).items():
        m = (d or {}).get("measurement", {})
        if m.get("type", "").lower() in ("ki", "kd") and m.get("value") and m.get("unit") in scale:
            ki_M = float(m["value"]) * scale[m["unit"]]
            if ki_M > 0:
                out[lig] = _RT_KCAL * math.log(ki_M)
    return out


def discover_rbfe_edges(results_dir):
    """Return {(ligA,ligB): {'complex': {rep: dG}, 'solvent': {rep: dG}}}."""
    root = Path(results_dir).expanduser()
    edges = {}
    seen = set()
    for rj in sorted(root.rglob("rbfe_*_result.json")):
        if "quickrun_cache" in rj.parts or rj.resolve() in seen:
            continue
        seen.add(rj.resolve())
        m = _EDGE_RE.match(rj.stem)
        if not m:
            continue
        ligA, legA, ligB, legB, rep = m.groups()
        if legA != legB:
            continue
        rep = int(rep)
        rec = parse_replicate(rj, f"rep{rep}")   # reuse robust JSON parser
        dG = rec.estimate
        edges.setdefault((ligA, ligB), {"complex": {}, "solvent": {}})
        edges[(ligA, ligB)][legA][rep] = dG
    return edges


def _mean_sem(vals):
    vals = [v for v in vals if isinstance(v, (int, float))]
    n = len(vals)
    if n == 0:
        return None, None, None, 0
    mean = sum(vals) / n
    if n == 1:
        return mean, None, None, 1
    sd = (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5
    return mean, sd / (n ** 0.5), sd, n


def compute_rbfe_results(edges):
    """Per edge: predicted ddG = dG_complex - dG_solvent (paired by replicate
    where possible, else mean-difference). Returns list of dicts."""
    rows = []
    for (ligA, ligB), legs in sorted(edges.items()):
        cx, sv = legs["complex"], legs["solvent"]
        # pair by replicate index
        shared = sorted(set(cx) & set(sv))
        if shared:
            per_rep = [cx[r] - sv[r] for r in shared
                       if isinstance(cx[r], (int, float)) and isinstance(sv[r], (int, float))]
            ddg, sem, sd, n = _mean_sem(per_rep)
        else:
            mc, _, _, _ = _mean_sem(list(cx.values()))
            ms, _, _, _ = _mean_sem(list(sv.values()))
            ddg = (mc - ms) if (mc is not None and ms is not None) else None
            sem = sd = None
            n = 0
        cmean = _mean_sem(list(cx.values()))[0]
        smean = _mean_sem(list(sv.values()))[0]
        rows.append({"edge": (ligA, ligB),
                     "ligA": ligA, "ligB": ligB,
                     "ddg_pred": ddg, "sem": sem, "sd": sd,
                     "n_reps": n if n else max(len(cx), len(sv)),
                     "dG_complex": cmean, "dG_solvent": smean,
                     "n_complex": len(cx), "n_solvent": len(sv)})
    return rows


def _stats_vs_exp(pairs):
    """pairs = [(pred, exp), ...] -> MUE, RMSE, Pearson R, Kendall tau."""
    pairs = [(p, e) for p, e in pairs if p is not None and e is not None]
    n = len(pairs)
    if n == 0:
        return {}
    pr = [p for p, _ in pairs]; ex = [e for _, e in pairs]
    mue = sum(abs(p - e) for p, e in pairs) / n
    rmse = (sum((p - e) ** 2 for p, e in pairs) / n) ** 0.5
    # Pearson
    mp = sum(pr) / n; me = sum(ex) / n
    cov = sum((p - mp) * (e - me) for p, e in pairs)
    vp = sum((p - mp) ** 2 for p in pr); ve = sum((e - me) ** 2 for e in ex)
    R = cov / ((vp * ve) ** 0.5) if vp > 0 and ve > 0 else None
    # Kendall tau
    conc = disc = 0
    for i in range(n):
        for j in range(i + 1, n):
            a = (pr[i] - pr[j]); b = (ex[i] - ex[j])
            s = a * b
            if s > 0: conc += 1
            elif s < 0: disc += 1
    tau = (conc - disc) / (0.5 * n * (n - 1)) if n > 1 else None
    return {"n": n, "mue": mue, "rmse": rmse, "pearson_r": R, "kendall_tau": tau}



def compute_network_mle(rows):
    """Weighted least-squares MLE: integrate per-edge ddG into per-ligand
    relative dG (mean-centered), with per-ligand SEM from the fit covariance.

    rows: list of dicts with ligA, ligB, ddg_pred, sem (from compute_rbfe_results).
    Returns {ligand: {"dG": float, "sem": float}} or {} if it cannot be solved.
    Pure NumPy; degrades to {} (never raises) if numpy is missing or the system
    is singular.
    """
    usable = [r for r in rows
              if r.get("ddg_pred") is not None and r.get("ligA") and r.get("ligB")]
    if not usable:
        return {}
    try:
        import numpy as np
    except Exception:
        return {}
    nodes = sorted({n for r in usable for n in (r["ligA"], r["ligB"])})
    idx = {n: i for i, n in enumerate(nodes)}
    N, E = len(nodes), len(usable)
    A = np.zeros((E + 1, N))
    y = np.zeros(E + 1)
    w = np.zeros(E + 1)
    for e, r in enumerate(usable):
        A[e, idx[r["ligB"]]] += 1.0
        A[e, idx[r["ligA"]]] -= 1.0
        y[e] = r["ddg_pred"]
        s = r.get("sem") or 0.0
        w[e] = 1.0 / (max(s, 1e-3) ** 2)
    # gauge: sum(dG) = 0, heavily weighted
    A[E, :] = 1.0
    y[E] = 0.0
    w[E] = 1e6
    W = np.diag(w)
    try:
        AtW = A.T @ W
        M = AtW @ A
        dG = np.linalg.solve(M, AtW @ y)
        cov = np.linalg.inv(M)
        sem = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    except Exception:
        return {}
    return {nodes[i]: {"dG": float(dG[i]), "sem": float(sem[i])} for i in range(N)}


def compute_cycle_closure(rows):
    """Find independent cycles in the edge network and report hysteresis
    (sum of signed ddG around each loop, which should be ~0).

    Returns {"n_nodes","n_edges","n_independent_cycles","loops":[{...}],
             "mean_abs","max_abs"}. No experimental data needed.
    """
    usable = [r for r in rows
              if r.get("ddg_pred") is not None and r.get("ligA") and r.get("ligB")]
    nodes = sorted({n for r in usable for n in (r["ligA"], r["ligB"])})
    info = {"n_nodes": len(nodes), "n_edges": len(usable),
            "n_independent_cycles": 0, "loops": [], "mean_abs": None,
            "max_abs": None}
    if not usable:
        return info
    import collections
    adj = collections.defaultdict(list)
    for e, r in enumerate(usable):
        adj[r["ligA"]].append((r["ligB"], e, +1))
        adj[r["ligB"]].append((r["ligA"], e, -1))
    # spanning forest via BFS; build node potentials from tree edges
    parent, parent_edge, pot = {}, {}, {}
    tree = set()
    for root in nodes:
        if root in pot:
            continue
        pot[root] = 0.0
        stack = [root]
        while stack:
            u = stack.pop()
            for v, e, sign in adj[u]:
                if v not in pot:
                    pot[v] = pot[u] + (usable[e]["ddg_pred"] * sign)
                    parent[v] = u
                    parent_edge[v] = e
                    tree.add(e)
                    stack.append(v)
    loops = []
    for e, r in enumerate(usable):
        if e in tree:
            continue
        closure = r["ddg_pred"] - (pot[r["ligB"]] - pot[r["ligA"]])
        loops.append({"edge": f"{r['ligA']} -> {r['ligB']}",
                      "hysteresis": float(closure)})
    info["n_independent_cycles"] = len(loops)
    info["loops"] = loops
    if loops:
        mags = [abs(l["hysteresis"]) for l in loops]
        info["mean_abs"] = sum(mags) / len(mags)
        info["max_abs"] = max(mags)
    return info


def render_rbfe_report(results_dir, out_dir, experimental=None):
    """Build a single-page RBFE campaign HTML (edge table + experiment stats)."""
    edges = discover_rbfe_edges(results_dir)
    rows = compute_rbfe_results(edges)
    exp = experimental or {}

    for r in rows:
        r["ddg_exp"] = (exp[r["ligB"]] - exp[r["ligA"]]
                        if r["ligA"] in exp and r["ligB"] in exp else None)
        r["abs_err"] = (abs(r["ddg_pred"] - r["ddg_exp"])
                        if r["ddg_pred"] is not None and r["ddg_exp"] is not None else None)

    stats = _stats_vs_exp([(r["ddg_pred"], r["ddg_exp"]) for r in rows]) if exp else {}

    # ---- network analysis (no experimental data required) ----
    mle = compute_network_mle(rows)            # {lig: {dG, sem}}
    closure = compute_cycle_closure(rows)      # hysteresis around loops

    # per-ligand stats vs experiment (only if experiment supplied AND MLE solved)
    lig_stats = {}
    lig_rank_rows = []
    if mle:
        # center experiment on the same set of ligands present in the MLE
        exp_common = {k: exp[k] for k in mle if k in exp}
        exp_mean = (sum(exp_common.values()) / len(exp_common)) if exp_common else None
        for lig, d in mle.items():
            e_c = (exp[lig] - exp_mean) if (exp_mean is not None and lig in exp) else None
            err = abs(d["dG"] - e_c) if e_c is not None else None
            lig_rank_rows.append({"lig": lig, "dG": d["dG"], "sem": d["sem"],
                                  "exp": e_c, "abs_err": err})
        if exp_common and len(exp_common) >= 2:
            lig_stats = _stats_vs_exp(
                [(r["dG"], r["exp"]) for r in lig_rank_rows if r["exp"] is not None])

    # ---- table ----
    trs = ""
    for r in sorted(rows, key=lambda x: (x["ligA"], x["ligB"])):
        edge = f"{r['ligA']} &rarr; {r['ligB']}"
        pred = (f"{r['ddg_pred']:+.2f}" + (f" +/- {r['sem']:.2f}" if r['sem'] else "")
                if r["ddg_pred"] is not None else "n/a")
        expv = f"{r['ddg_exp']:+.2f}" if r["ddg_exp"] is not None else "-"
        err = (f"{r['abs_err']:.2f}" if r["abs_err"] is not None else "-")
        ecl = ""
        if r["abs_err"] is not None:
            ecl = "good" if r["abs_err"] <= 1.0 else ("warn" if r["abs_err"] <= 2.0 else "bad")
        trs += (f"<tr><td style='text-align:left'>{edge}</td>"
                f"<td>{r['n_reps']}</td>"
                f"<td>{_fnum(r['dG_complex'])}</td>"
                f"<td>{_fnum(r['dG_solvent'])}</td>"
                f"<td><b>{pred}</b></td><td>{expv}</td>"
                f"<td class='{ecl}'>{err}</td></tr>")

    stat_cards = ""
    if stats:
        def card(v, l, n=2):
            return (f"<div class='card'><div class='v'>"
                    f"{v:.{n}f}" if isinstance(v, (int, float)) else
                    f"<div class='card'><div class='v'>n/a") + f"</div><div class='l'>{l}</div></div>"
        stat_cards = (
            "<div class='kpi'>"
            f"<div class='card'><div class='v'>{stats['n']}</div><div class='l'>edges vs exp</div></div>"
            f"<div class='card'><div class='v'>{stats['mue']:.2f}</div><div class='l'>MUE (kcal/mol)</div></div>"
            f"<div class='card'><div class='v'>{stats['rmse']:.2f}</div><div class='l'>RMSE (kcal/mol)</div></div>"
            f"<div class='card'><div class='v'>{(('%.2f'%stats['pearson_r']) if stats['pearson_r'] is not None else 'n/a')}</div><div class='l'>Pearson R</div></div>"
            f"<div class='card'><div class='v'>{(('%.2f'%stats['kendall_tau']) if stats['kendall_tau'] is not None else 'n/a')}</div><div class='l'>Kendall tau</div></div>"
            "</div>")

    # ---- scatter plot pred vs exp ----
    scatter = ""
    if exp:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        pts = [(r["ddg_exp"], r["ddg_pred"]) for r in rows
               if r["ddg_exp"] is not None and r["ddg_pred"] is not None]
        if pts:
            xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
            fig, ax = plt.subplots(figsize=(5, 5))
            lim = max(abs(min(xs+ys)), abs(max(xs+ys))) + 0.5
            ax.plot([-lim, lim], [-lim, lim], color="#888", lw=1)
            ax.fill_between([-lim, lim], [-lim-1, lim-1], [-lim+1, lim+1],
                            color="#4c78a8", alpha=0.10, label="+/-1 kcal/mol")
            ax.scatter(xs, ys, c="#e45756", s=42, edgecolor="k", zorder=3)
            ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
            ax.set_xlabel("Experimental ddG (kcal/mol)")
            ax.set_ylabel("Predicted ddG (kcal/mol)")
            ax.set_title("RBFE predicted vs experimental ddG")
            ax.legend(fontsize=8); ax.grid(alpha=0.3)
            scatter = "<h2>Predicted vs Experimental</h2><img src='" + _fig_to_b64(fig) + "'>"

    # ---- per-ligand MLE ranking table ----
    mle_section = ""
    if lig_rank_rows:
        rk = sorted(lig_rank_rows, key=lambda r: r["dG"])
        body = ""
        for i, r in enumerate(rk, 1):
            ev = f"{r['exp']:+.2f}" if r["exp"] is not None else "-"
            er = f"{r['abs_err']:.2f}" if r["abs_err"] is not None else "-"
            ecl = ""
            if r["abs_err"] is not None:
                ecl = "good" if r["abs_err"] <= 1.0 else ("warn" if r["abs_err"] <= 2.0 else "bad")
            body += (f"<tr><td>{i}</td><td style='text-align:left'>{r['lig']}</td>"
                     f"<td><b>{r['dG']:+.2f}</b></td><td>{r['sem']:.2f}</td>"
                     f"<td>{ev}</td><td class='{ecl}'>{er}</td></tr>")
        ls_cards = ""
        if lig_stats:
            ls_cards = (
                "<div class='kpi'>"
                f"<div class='card'><div class='v'>{lig_stats['n']}</div><div class='l'>ligands vs exp</div></div>"
                f"<div class='card'><div class='v'>{lig_stats['mue']:.2f}</div><div class='l'>MUE (kcal/mol)</div></div>"
                f"<div class='card'><div class='v'>{lig_stats['rmse']:.2f}</div><div class='l'>RMSE (kcal/mol)</div></div>"
                f"<div class='card'><div class='v'>{(('%.2f'%lig_stats['pearson_r']) if lig_stats['pearson_r'] is not None else 'n/a')}</div><div class='l'>Pearson R</div></div>"
                f"<div class='card'><div class='v'>{(('%.2f'%lig_stats['kendall_tau']) if lig_stats['kendall_tau'] is not None else 'n/a')}</div><div class='l'>Kendall tau</div></div>"
                "</div>")
        mle_section = (
            "<h2>Per-ligand ranking (network MLE)</h2>"
            "<p class='note'>Per-edge ddG values integrated into relative per-ligand "
            "dG by weighted least squares (weights = 1/SEM&sup2;, gauge: mean dG = 0). "
            "This ranking, not the raw edge ddG, is what to use for prioritisation.</p>"
            + ls_cards +
            "<table><tr><th>Rank</th><th>Ligand</th><th>dG (rel)</th><th>SEM</th>"
            "<th>exp dG (rel)</th><th>|error|</th></tr>" + body + "</table>")

    # ---- cycle-closure section ----
    closure_section = ""
    if closure["n_edges"]:
        if closure["n_independent_cycles"] == 0:
            closure_section = (
                "<h2>Cycle closure</h2><p class='note'>The network is a "
                f"<b>spanning tree</b> ({closure['n_nodes']} ligands, "
                f"{closure['n_edges']} edges, 0 independent cycles). There is no "
                "redundancy, so cycle-closure cannot be assessed and each ligand's "
                "dG rests on a single path. Consider adding closing edges for "
                "internal error estimates.</p>")
        else:
            lb = "".join(
                f"<tr><td style='text-align:left'>{l['edge']}</td>"
                f"<td class='{('good' if abs(l['hysteresis'])<=0.5 else ('warn' if abs(l['hysteresis'])<=1.0 else 'bad'))}'>"
                f"{l['hysteresis']:+.2f}</td></tr>" for l in closure["loops"])
            closure_section = (
                "<h2>Cycle closure</h2>"
                f"<p class='note'>{closure['n_independent_cycles']} independent "
                f"cycle(s). Mean |hysteresis| = {closure['mean_abs']:.2f}, "
                f"max = {closure['max_abs']:.2f} kcal/mol. Values near 0 indicate "
                "internally consistent edges (no experiment needed).</p>"
                "<table><tr><th>Loop (closed by edge)</th><th>hysteresis (kcal/mol)</th></tr>"
                + lb + "</table>")

    html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>OpenFE RBFE campaign</title><style>{_CSS}</style></head><body>
<h1>OpenFE RBFE Campaign Report</h1>
<p><b>{len(rows)} transformations (edges)</b> &nbsp;|&nbsp; Generated: {datetime.now():%Y-%m-%d %H:%M}</p>
<p class="note">Each edge's result is <b>ddG = dG(complex) - dG(solvent)</b>, averaged over
replicates (paired by replicate index when available). Per-leg dG values alone are not
binding free energies.</p>
{stat_cards}
{mle_section}
{closure_section}
{scatter}
<h2>Edge results</h2>
<table><tr><th>Transformation (A &rarr; B)</th><th>n reps</th>
<th>dG complex</th><th>dG solvent</th><th>ddG pred</th><th>ddG exp</th><th>|error|</th></tr>
{trs}</table>
<div class="foot">Auto-generated by amber_md.gui.openfe_report (RBFE mode).
Experimental ddG from Ki via dG = RT ln(Ki), RT={_RT_KCAL:.3f} kcal/mol at 298.15 K.
|error| coloring: green &le;1, amber &le;2, red &gt;2 kcal/mol.</div>
</body></html>"""

    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    idx = out_dir / "rbfe_index.html"
    idx.write_text(html)
    return {"index": str(idx), "rows": rows, "stats": stats,
            "mle": mle, "cycle_closure": closure, "ligand_stats": lig_stats,
            "ligand_ranking": sorted(lig_rank_rows, key=lambda r: r["dG"])}
