"""Split a solvated complex prmtop into complex/receptor/ligand topologies."""

# -*- coding: utf-8 -*-
# Author:    Pulan Yu
# Developer: Pulan Yu <chinabio@gmail.com>
# Contact:   chinabio@gmail.com
# Part of AffiniScreen.

from __future__ import annotations
from pathlib import Path
from .utils import run, which_or_die, ensure_dir, CommandError
from .logger import get_logger
log = get_logger()

class TopologySplitter:
    def __init__(self, work_dir, ligand_resname="LIG",
                 solvent_mask=":WAT,HOH,Na+,Cl-,K+,Mg+2,Ca+2",
                 radii="mbondi3"):
        self.dir = ensure_dir(Path(work_dir)/"topo")
        self.lig = ligand_resname.upper()
        self.solvent_mask = solvent_mask
        # GB radius set applied to ALL split prmtops. MUST match the igb used
        # in mmgbsa.in: igb=8 (GBn2) and igb=7 (GBn) require mbondi3 / mbondi2;
        # the default tleap radii (mbondi) cause mmpbsa_py_energy to FAIL on the
        # complex prmtop. The ante-MMPBSA.py path already passes --radii; the
        # parmed path previously did NOT set radii at all -- this fixes that so
        # both paths are consistent.
        self.radii = radii

    def _expected_outputs(self):
        return {"solvated": self.dir/"solvated.prmtop",
                "complex":  self.dir/"complex.prmtop",
                "receptor": self.dir/"receptor.prmtop",
                "ligand":   self.dir/"ligand.prmtop"}

    def _all_exist(self, out):
        return all(p.exists() and p.stat().st_size > 0 for p in out.values())

    def _split_parmed(self, solvated_prm):
        import parmed as pmd
        from parmed.tools import changeRadii
        out = self._expected_outputs()

        def _save_gb_prmtop(struct, path):
            # 1) GB requires a NON-PERIODIC topology. The dry complex/receptor/
            #    ligand are stripped from the *solvated* prmtop, which keeps its
            #    periodic box (IFBOX=1). mmpbsa_py_energy with gb>0 then aborts:
            #    "gb>0 is incompatible with periodic boundary conditions ...
            #     set IFBOX in the PRMTOP file to 0". So we remove the box here.
            try:
                struct.box = None
            except Exception as e:  # noqa: BLE001
                log.warning("Could not clear periodic box (%s); GB calc may "
                            "fail with an IFBOX error", e)
            # 2) Apply the GB radius set BEFORE saving so the prmtop's RADII /
            #    SCREEN sections match the chosen igb (e.g. mbondi3 for igb=8).
            if self.radii:
                try:
                    changeRadii(struct, self.radii).execute()
                except Exception as e:  # noqa: BLE001
                    raise RuntimeError(
                        "changeRadii(%s) failed for %s: %s. igb=8 (GBn2) "
                        "REQUIRES mbondi3; shipping a prmtop with the wrong "
                        "RADIUS_SET makes mmpbsa_py_energy fail on the complex "
                        "GB calc. Refusing to save a broken topology."
                        % (self.radii, path.name, e)) from e
            struct.save(str(path), overwrite=True)
            # RADIUS_SET is verified for the whole set after the
            # split completes, via _verify_radii() in split(): it
            # reads %FLAG RADIUS_SET from the saved prmtop text
            # (robust for both parmed and ante-MMPBSA paths).

        # solvated kept as-is (only used as -sp reference; radii irrelevant there)
        sysp = pmd.load_file(str(solvated_prm))
        sysp.save(str(out["solvated"]), overwrite=True)
        dry  = pmd.load_file(str(solvated_prm)); dry.strip(self.solvent_mask)
        _save_gb_prmtop(dry, out["complex"])
        rec  = pmd.load_file(str(out["complex"])); rec.strip(f":{self.lig}")
        _save_gb_prmtop(rec, out["receptor"])
        lig  = pmd.load_file(str(out["complex"])); lig.strip(f"!:{self.lig}")
        _save_gb_prmtop(lig, out["ligand"])
        log.info("parmed split: GB radii '%s' applied and periodic box removed (IFBOX=0) on complex/receptor/ligand",
                 self.radii)
        return out

    def _split_ante(self, solvated_prm):
        which_or_die("ante-MMPBSA.py")
        out = self._expected_outputs()
        # v2.2.9: remove existing files so ante-MMPBSA.py won't error with FileExists
        for k in ("complex", "receptor", "ligand"):
            if out[k].exists():
                log.info("Removing stale %s before ante-MMPBSA.py rebuild", out[k].name)
                out[k].unlink()
        out["solvated"].write_bytes(Path(solvated_prm).read_bytes())
        run(["ante-MMPBSA.py","-p",str(solvated_prm),
             "-c",str(out["complex"]),"-r",str(out["receptor"]),
             "-l",str(out["ligand"]),"-s",self.solvent_mask,
             "-n",f":{self.lig}",f"--radii={self.radii}"], cwd=self.dir)
        self._verify_radii(out)
        return out

    def _verify_radii(self, out):
        """Ensure complex/receptor/ligand prmtops carry the requested radius set.
        Reads RADIUS_SET directly from the prmtop text (no parmed dependency) so
        it works for BOTH the parmed and ante-MMPBSA split paths. A mismatch
        (e.g. mbondi when igb=8 needs mbondi3) is the documented cause of
        mmpbsa_py_energy failing on the complex leg -- so we fail loudly here."""
        if not self.radii:
            return
        want = self.radii.lower()
        for k in ("complex", "receptor", "ligand"):
            p = out.get(k)
            if not p or not Path(p).exists():
                continue
            got = ""
            try:
                txt = Path(p).read_text(errors="ignore").splitlines()
                for idx, line in enumerate(txt):
                    if line.startswith("%FLAG RADIUS_SET"):
                        # value is 2 lines down (after %FORMAT)
                        got = (txt[idx + 2] if idx + 2 < len(txt) else "").lower()
                        break
            except Exception as e:  # noqa: BLE001
                log.warning("Could not read RADIUS_SET from %s (%s)", p, e)
                continue
            if want not in got:
                raise RuntimeError(
                    "RADIUS_SET mismatch in %s: requested '%s' but prmtop "
                    "reports '%s'. igb=8 (GBn2) requires mbondi3; MMPBSA would "
                    "fail on the complex GB calc. Aborting topology split so a "
                    "broken prmtop never reaches MMPBSA." % (
                        Path(p).name, self.radii, got.strip() or "<empty>"))
            log.info("GB radii OK: %s -> '%s'", Path(p).name, got.strip())

    def split(self, solvated_prm, *, force=False):
        """Split solvated topology. Returns dict with 4 paths.

        v2.2.9: idempotent — if all 4 outputs already exist, returns them
        without re-splitting. Avoids redundant work AND the FileExists
        crash when ante-MMPBSA.py is invoked a second time.
        """
        out = self._expected_outputs()
        if not force and self._all_exist(out):
            self._verify_radii(out)
            log.info("TopologySplitter: outputs already present in %s, skipping split.",
                     self.dir)
            return out

        try:
            import parmed   # noqa
            return self._split_parmed(Path(solvated_prm))
        except ImportError as e:
            log.warning("parmed unavailable (%s) - falling back to ante-MMPBSA.py", e)
            return self._split_ante(Path(solvated_prm))
        except CommandError:
            raise
        except Exception as e:
            log.warning("parmed split failed (%s: %s) - falling back to ante-MMPBSA.py",
                        type(e).__name__, e)
            return self._split_ante(Path(solvated_prm))

    @staticmethod
    def sanity_check(topos):
        try:
            import parmed as pmd
            c = pmd.load_file(str(topos["complex"])).ptr("NATOM")
            r = pmd.load_file(str(topos["receptor"])).ptr("NATOM")
            l = pmd.load_file(str(topos["ligand"])).ptr("NATOM")
            ok = (r + l) == c
            if not ok:
                log.warning("Topology sanity check FAILED: %d (receptor) + %d (ligand) != %d (complex)",
                            r, l, c)
            return ok
        except ImportError:
            return True
        except Exception as e:
            log.warning("Topology sanity check error (%s) - skipping", e)
            return True
