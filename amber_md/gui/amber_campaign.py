# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Part of AffiniScreen.
"""Amber-TI RBFE campaign adapter (final61). Mirrors openfe_campaign's API so
the shared FEP Campaign page can drive either engine behind an engine selector.

On-disk layout (Setup & Launch -> RBFE/Amber + fep_driver --mode rbfe):
    <wd>/rbfe_map/{edges.csv,mapping.json,diagnostics.txt,map.graphml}
    <wd>/RBFE_CAMPAIGN.json
    <wd>/edges/<A~B>/fep/complex/summary.json  -> dG_kcal_mol
    <wd>/edges/<A~B>/fep/solvent/summary.json  -> dG_kcal_mol
    ddG(A->B) = dG_complex - dG_solvent
"""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import sys

from .openfe_campaign import read_edges_table, mapping_network_mermaid
from .openfe_common import count_my_jobs

JOB_PREFIX = "amb_rbfe"


@dataclass
class CampaignLayout:
    work_dir: Path
    network_dir: Path
    results_dir: Path
    log_dir: Path
    dg_tsv: Path

    @classmethod
    def at(cls, work_dir: Path) -> "CampaignLayout":
        wd = Path(work_dir).expanduser()
        return cls(work_dir=wd, network_dir=wd / "rbfe_map",
                   results_dir=wd / "edges", log_dir=wd / "logs",
                   dg_tsv=wd / "amber_rbfe_ddg.tsv")


def _edge_name(row: dict) -> str:
    return row.get("edge") or f"{row.get('lig_a','?')}~{row.get('lig_b','?')}"


def _edge_rows(layout: CampaignLayout) -> list[dict]:
    csv_path = layout.network_dir / "edges.csv"
    if csv_path.exists():
        try:
            return read_edges_table(csv_path)
        except Exception:
            pass
    rows = []
    if layout.results_dir.is_dir():
        for dd in sorted(p for p in layout.results_dir.iterdir() if p.is_dir()):
            a, _, b = dd.name.partition("~")
            rows.append({"edge": dd.name, "lig_a": a, "lig_b": b})
    return rows


def edges(layout: CampaignLayout) -> list[Path]:
    return [layout.results_dir / _edge_name(r) for r in _edge_rows(layout)]


def task_list(layout: CampaignLayout, n_replicates: int):
    return [(e, r) for e in edges(layout) for r in range(max(1, n_replicates))]


def _leg_dg(edge_dir: Path, leg: str):
    sj = edge_dir / "fep" / leg / "summary.json"
    if not sj.exists():
        return None
    try:
        data = json.loads(sj.read_text())
    except Exception:
        return None
    if isinstance(data, dict) and "dG_kcal_mol" in data:
        return data.get("dG_kcal_mol")
    if isinstance(data, dict) and leg in data and isinstance(data[leg], dict):
        return data[leg].get("dG_kcal_mol")
    return None


def edge_ddg(edge_dir: Path):
    c = _leg_dg(edge_dir, "complex")
    s = _leg_dg(edge_dir, "solvent")
    if c is None or s is None:
        return None, None
    return (c - s), 0.5


def edge_status(edge_dir: Path) -> str:
    ddg, _ = edge_ddg(edge_dir)
    if ddg is not None:
        return "DONE"
    if (edge_dir / "fep" / "complex").is_dir() or (edge_dir / "fep" / "solvent").is_dir():
        return "running"
    if edge_dir.exists():
        return "analyzing"
    return "not started"


def edge_ddg_table(layout: CampaignLayout):
    out = []
    for row in _edge_rows(layout):
        name = _edge_name(row)
        ed = layout.results_dir / name
        ddg, err = edge_ddg(ed)
        if ddg is None:
            continue
        a = row.get("lig_a") or name.partition("~")[0]
        b = row.get("lig_b") or name.partition("~")[2]
        out.append((str(a), str(b), float(ddg), float(err)))
    return out


def campaign_status(layout: CampaignLayout, n_replicates: int) -> dict:
    rows, buckets = [], {"DONE": 0, "analyzing": 0, "running": 0, "not started": 0}
    eds = edges(layout)
    for ed in eds:
        ddg, _ = edge_ddg(ed)
        stt = edge_status(ed)
        buckets[stt] = buckets.get(stt, 0) + 1
        rows.append({"edge": ed.name, "status": stt,
                     "ddG (kcal/mol)": (round(ddg, 2) if ddg is not None else None)})
    total = len(eds); done = buckets.get("DONE", 0)
    return {"pct": (done / total) if total else 0.0, "done": done,
            "total": total, "buckets": buckets, "rows": rows}


def submit_campaign(layout: CampaignLayout, settings, n_replicates: int,
                    progress_cb=None, throttle_n: int = 0,
                    protein: Path = None, ligands: Path = None,
                    python_bin: str = "python", only_edges=None):
    if protein is None or ligands is None:
        return 0, 0, ["Amber RBFE submit needs protein + ligands (the files the "
                      "map was planned from)."]
    cmd = [python_bin, "-m", "amber_md.fep_driver", "--mode", "rbfe",
           "--work-dir", str(layout.work_dir),
           "--protein-pdb", str(Path(protein)),
           "--ligands", str(Path(ligands)),
           "--rbfe-map-dir", str(layout.network_dir),
           "--analyze", "--submit"]
    if getattr(settings, "queue", None):
        cmd += ["--queue", str(settings.queue)]
    if getattr(settings, "walltime", None):
        cmd += ["--walltime", str(settings.walltime)]
    if getattr(settings, "project", None):
        cmd += ["--project", str(settings.project)]
    if only_edges:
        cmd += ["--only-edges", *list(only_edges)]
    layout.log_dir.mkdir(parents=True, exist_ok=True)
    log = layout.log_dir / "amber_rbfe_campaign.log"
    try:
        from .common import spawn_detached
        spawn_detached(cmd, log, cwd=str(layout.work_dir))
    except Exception as e:  # noqa: BLE001
        return 0, 0, [f"spawn failed: {e}"]
    n = len(edges(layout))
    if progress_cb:
        progress_cb(n, n, n, 0)
    return n, 0, []


def submit_preview_cmd(layout, settings, protein="<protein.pdb>",
                       ligands="<ligands.sdf>", python_bin="python"):
    cmd = [python_bin, "-m", "amber_md.fep_driver", "--mode", "rbfe",
           "--work-dir", str(layout.work_dir),
           "--protein-pdb", str(protein), "--ligands", str(ligands),
           "--rbfe-map-dir", str(layout.network_dir), "--analyze", "--submit"]
    if getattr(settings, "queue", None):
        cmd += ["--queue", str(settings.queue)]
    if getattr(settings, "walltime", None):
        cmd += ["--walltime", str(settings.walltime)]
    if getattr(settings, "project", None):
        cmd += ["--project", str(settings.project)]
    return cmd


def gather_command(layout: CampaignLayout, report: str = "ddg") -> list[str]:
    return [sys.executable or "python", "-m", "amber_md.gui.amber_campaign",
            "--gather", str(layout.work_dir)]


def gather_now(layout: CampaignLayout) -> Path:
    rows = edge_ddg_table(layout)
    layout.dg_tsv.parent.mkdir(parents=True, exist_ok=True)
    with open(layout.dg_tsv, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["ligand_a", "ligand_b", "ddG(kcal/mol)", "uncertainty"])
        for a, b, ddg, err in rows:
            w.writerow([a, b, f"{ddg:.4f}", f"{err:.4f}"])
    return layout.dg_tsv


def find_mapping_artifacts(layout: CampaignLayout) -> dict:
    found, nd = {}, layout.network_dir
    if (nd / "diagnostics.txt").exists(): found["diagnostics"] = nd / "diagnostics.txt"
    if (nd / "edges.csv").exists():       found["edges"] = nd / "edges.csv"
    if (nd / "mapping.json").exists():    found["mapping"] = nd / "mapping.json"
    if (nd / "map.graphml").exists():     found["graphml"] = nd / "map.graphml"
    if edges(layout):                     found["transformations_dir"] = nd
    return found


def preview_script(layout: CampaignLayout, settings) -> str:
    cmd = submit_preview_cmd(layout, settings)
    return ("# Amber-TI RBFE campaign (driver builds dual-topology + bsubs all\n"
            "# per-edge legs; --analyze writes each edge's fep/*/summary.json).\n"
            + " ".join(str(c) for c in cmd) + "\n")


def _main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description="Amber RBFE campaign helper.")
    ap.add_argument("--gather", type=Path, help="work_dir to gather ddG from")
    a = ap.parse_args(argv)
    if a.gather:
        print(f"Wrote {gather_now(CampaignLayout.at(a.gather))}")
        return 0
    ap.print_help(); return 1


if __name__ == "__main__":
    raise SystemExit(_main())
