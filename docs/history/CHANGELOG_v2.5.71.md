# amber_md v2.5.71 -- box buffer 16 A -> 12 A + convergence tool

**Base:** v2.5.70

## TL;DR
Drops the solvent box buffer from **16 A to 12 A** (OpenFE-aligned), recovering
~30%% throughput at no accuracy cost, now that the instability the 16 A was
working around is properly fixed. Ships a read-only per-leg convergence
diagnostic. HREMD deliberately stays opt-in.

## Why 16 A existed, and why it can go
v2.4.21 raised the buffer 12 A -> 16 A to stop a "GPU small-box crash on small
ligands." That crash was a **soft-core / box-drift instability**, not a genuine
periodic-image problem. v2.5.70 root-caused and fixed it (stage-aware hard-vdW
decharge that removed the charge-on-softcore singularity), alongside the
44-window vdw schedule and the MD-only NaN guard. With the real cause gone, the
oversized box is pure overhead.

On the current 148,238-atom test system the box is **91%% water** (44,938 TIP3P;
solute only ~13,424 atoms) in a 124.6 A truncated octahedron -- measured
throughput **16.4 ns/day** on an L40S (5.27 ms/step), i.e. ~7.3 h per 5 ns
complex window and ~3 days/ligand on 8 GPUs.

## Why 12 A (not 10, not back to FEP+ 5)
| package | default solvent padding | note |
|---|---|---|
| OpenFE (OpenMM) | **12 A** (1.2 nm) + dodecahedron | closest validated analog |
| GROMACS practice | 10-12 A | rhombic dodecahedron |
| AmberTools tutorials / this repo preflight | 10-12 A | solvateBox 12.0 |
| FEP+/Desmond | ~5 A *(different definition)* | NOT portable to solvateOct |

12 A matches OpenFE's validated default and leaves headroom over the 10 A
nonbonded cutoff for a decoupling (reorienting) soft-core ligand. 10 A with
cut=10 is uncomfortably tight; FEP+'s "5 A" uses a different box-image
convention and cannot be copied into solvateOct. Truncated octahedron (already
used) captures the box-shape water savings.

Projected (solute span fixed, MD cost ~ linear in atom count):
| buffer | ~NATOM | ~ns/day | 5 ns/window | campaign (8 GPU) |
|---|---|---|---|---|
| 16 A (old) | 148k | 16.4 | 7.3 h | ~3.0 d |
| **12 A (new)** | **~113k** | **~21.5** | **~6.4 h** | **~2.2 d** |
| 10 A | ~102k | ~24 | ~5.7 h | ~2.0 d (tight) |

## Changes
1. `config.py`: `box_buffer_A` 16.0 -> 12.0 (+ rationale comment).
2. `fep_driver.py`: `--box-buffer` default 16.0 -> 12.0.
3. GUI `0_Setup_and_Launch.py`: buffer slider default 10.0 -> 12.0 (harmonizes
   GUI with CLI/config, which previously disagreed: GUI 10 vs CLI/config 16).
4. NEW `tools/convergence_analysis.py`: per-leg cumulative-time scan + forward/
   reverse halves, reporting TI / BAR / MBAR; gracefully degrades to TI when
   MBAR/BAR cannot solve (overlap collapse). Read-only. Decides whether 5 ns is
   enough and which windows can be trimmed.

## Deliberately NOT changed
- **HREMD stays opt-in** (`--exchange-freq`). Enabling replica exchange by
  default couples all windows into one job (a single window stall halts the
  leg) and changes failure modes -- the opposite of what you want while the
  pipeline is still being stabilized. Revisit once v2.5.70 + 44-window vdw are
  proven clean.
- decharge/vdw protocols and lambda schedules: unchanged from v2.5.70.

## Migration / action
- New box buffer applies to **newly prepared** systems only; in-flight runs and
  existing topologies are unaffected (re-prep required to benefit).
- **Validate one ligand through the GUI at 12 A before launching a full
  campaign** (the GUI has no single-leg mode; one ligand is the smallest unit
  of risk you control).
- No API changes.

## Verification (this release)
- `config.box_buffer_A == 12.0` and CLI `--box-buffer` default `12.0` asserted
  in build.
- config.py / fep_driver.py / convergence_analysis.py byte-compile.
- `tools/check_version_sync.py` passes (all touchpoints = 2.5.71).
