"""Snapshot generation: PyMOL → NGLview → MDAnalysis (fallback)."""

# -*- coding: utf-8 -*-
# Author:    AffiniScreen contributors
# Maintainer: AffiniScreen team
# Contact:   (see repository)
# Part of AffiniScreen.

from __future__ import annotations
import re
from pathlib import Path
from typing import Optional
from .utils import ensure_dir
from .logger import get_logger
log = get_logger()


# nglview 4.x writes "model_module_version": "4.0" into the embedded widget
# state, but the matching nglview-js-widgets JS bundle was never published to
# npm/jsDelivr beyond 3.1.5.  Without a require.config hint *and* a version
# downgrade, browsers fail with:
#   "Failed to load model class 'ColormakerRegistryModel' from module
#    'nglview-js-widgets' — Script error for 'nglview-js-widgets'"
# See https://github.com/nglviewer/nglview/issues/1172
_NGL_JS_VERSION = "3.1.5"  # latest version actually on the CDN
_NGL_REQUIRE_CONFIG = f"""
<script>
  // Patch added by Visualizer.nglview_html(): tell require.js where to find
  // the nglview-js-widgets AMD bundle (not auto-resolved by embed-amd.js).
  require.config({{
    paths: {{
      'nglview-js-widgets':
        'https://cdn.jsdelivr.net/npm/nglview-js-widgets@{_NGL_JS_VERSION}/dist/index'
    }}
  }});
</script>
"""
# Matches "model_module"/"view_module":"nglview-js-widgets","..._version":"X"
_NGL_VER_RE_A = re.compile(
    r'("(?:model|view)_module"\s*:\s*"nglview-js-widgets"\s*,\s*'
    r'"(?:model|view)_module_version"\s*:\s*")([^"]+)(")'
)
# Same pair in reverse order
_NGL_VER_RE_B = re.compile(
    r'("(?:model|view)_module_version"\s*:\s*")([^"]+)'
    r'("\s*,\s*"(?:model|view)_module"\s*:\s*"nglview-js-widgets")'
)
_EMBED_AMD_TAG = (
    '<script src="https://cdn.jsdelivr.net/npm/@jupyter-widgets/'
    'html-manager@^1.0.1/dist/embed-amd.js" crossorigin="anonymous"></script>'
)


def _patch_nglview_html(path: Path) -> None:
    """Make an nglview.write_html() output actually load in a browser.

    1. Rewrites any nglview-js-widgets module_version that's newer than what's
       on the CDN down to the latest published one (currently 3.1.5).
    2. Injects a require.config block before embed-amd.js so require.js can
       resolve the 'nglview-js-widgets' module.
    """
    try:
        html = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return

    pin = lambda m: m.group(1) + _NGL_JS_VERSION + m.group(3)
    html, n_a = _NGL_VER_RE_A.subn(pin, html)
    html, n_b = _NGL_VER_RE_B.subn(pin, html)

    if _EMBED_AMD_TAG in html and "nglview-js-widgets" not in html.split(
            _EMBED_AMD_TAG)[0][-400:]:
        html = html.replace(_EMBED_AMD_TAG,
                            _NGL_REQUIRE_CONFIG + _EMBED_AMD_TAG)
        injected = True
    else:
        injected = False

    path.write_text(html, encoding="utf-8")
    log.info("Patched %s (version rewrites: %d, require.config injected: %s)",
             path.name, n_a + n_b, injected)


class Visualizer:
    def __init__(self, work_dir):
        self.dir = ensure_dir(Path(work_dir)/"viz")
    def pymol_snapshot(self, structure, ligand_resname="LIG",
                       out_png: Optional[Path]=None,
                       width=1200, height=900, ray=True):
        try:
            from pymol import cmd
        except ImportError:
            log.warning("PyMOL not available — skipping PyMOL snapshot.")
            return None
        out_png = out_png or self.dir/f"{Path(structure).stem}_pymol.png"
        cmd.reinitialize(); cmd.load(str(structure), "sys")
        cmd.hide("everything", "sys"); cmd.show("cartoon", "sys and polymer")
        cmd.color("gray80", "sys and polymer")
        cmd.show("sticks", f"sys and resn {ligand_resname}")
        cmd.color("magenta", f"sys and resn {ligand_resname} and elem C")
        cmd.select("pocket", f"byres (polymer and (sys within 5 of resn {ligand_resname}))")
        cmd.show("sticks", "pocket and not (name C+N+O)")
        cmd.color("cyan", "pocket and elem C")
        cmd.bg_color("white"); cmd.orient(f"resn {ligand_resname}")
        cmd.zoom(f"resn {ligand_resname}", 8)
        if ray: cmd.ray(width, height)
        cmd.png(str(out_png), width=width, height=height, dpi=150, ray=int(ray))
        return out_png
    def nglview_html(self, topology, trajectory=None, ligand_resname="LIG",
                     out_html: Optional[Path]=None):
        try:
            import nglview as nv
            import mdtraj as md
        except ImportError:
            log.warning("nglview/mdtraj not available — skipping NGLview HTML.")
            return None
        out_html = out_html or self.dir/f"{Path(topology).stem}_view.html"
        if trajectory:
            traj = md.load(str(trajectory), top=str(topology))
            view = nv.show_mdtraj(traj)
        else:
            view = nv.show_file(str(topology))
        view.clear_representations()
        view.add_cartoon("protein", color="lightgrey")
        view.add_ball_and_stick(f":{ligand_resname}", color="magenta")
        try:
            nv.write_html(str(out_html), [view])
            # Post-process so the file actually renders in a plain browser:
            # nglview 4.x ships a widget-state version that has no matching
            # JS bundle on any CDN, and embed-amd.js doesn't know where to
            # fetch nglview-js-widgets from. See _patch_nglview_html docstring.
            _patch_nglview_html(out_html)
        except Exception as exc:
            log.warning("nv.write_html failed (%s); writing placeholder.", exc)
            out_html.write_text(
                "<html><body><p>Open with nglview in Jupyter.</p></body></html>"
            )
        return out_html
    def mda_snapshot(self, topology, trajectory=None, ligand_resname="LIG",
                     out_png: Optional[Path]=None):
        try:
            import MDAnalysis as mda
            import matplotlib; matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D  # noqa
        except ImportError:
            log.warning("MDAnalysis/matplotlib not available — skipping fallback snapshot.")
            return None
        out_png = out_png or self.dir/f"{Path(topology).stem}_mda.png"
        u = (mda.Universe(str(topology), str(trajectory))
             if trajectory else mda.Universe(str(topology)))
        prot = u.select_atoms("protein and name CA")
        lig  = u.select_atoms(f"resname {ligand_resname}")
        fig = plt.figure(figsize=(7,7))
        ax  = fig.add_subplot(111, projection="3d")
        if len(prot):
            ax.scatter(prot.positions[:,0], prot.positions[:,1], prot.positions[:,2],
                       c="lightgray", s=8, label=f"protein Cα (n={len(prot)})")
        if len(lig):
            ax.scatter(lig.positions[:,0], lig.positions[:,1], lig.positions[:,2],
                       c="magenta", s=30, label=f"ligand {ligand_resname}")
        ax.set_xlabel("X (Å)"); ax.set_ylabel("Y (Å)"); ax.set_zlabel("Z (Å)")
        ax.legend(loc="upper right")
        plt.tight_layout(); plt.savefig(out_png, dpi=140); plt.close()
        return out_png
    def render_all(self, topology, trajectory=None, ligand_resname="LIG"):
        struct = topology
        if Path(topology).suffix == ".prmtop":
            cand = Path(topology).with_name("complex_solv.pdb")
            if cand.exists(): struct = cand
        out = {"pymol":   self.pymol_snapshot(struct, ligand_resname),
               "nglview": self.nglview_html(topology, trajectory, ligand_resname),
               "mda":     self.mda_snapshot(topology, trajectory, ligand_resname)}
        return {k:v for k,v in out.items() if v}
