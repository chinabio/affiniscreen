# v2.5.3 build final51: unify MM-GBSA production MD to 10 ns for BOTH engines

## Request
The MM-GBSA production length differed by engine (Amber 10 ns, OpenMM 5 ns) and
an old batch path still defaulted to 50 ns. Make every default the same so
Amber and OpenMM MM-GBSA scores are directly comparable. Decision: 10 ns prod
(equil stays 1 ns).

## All prod-ns defaults set to 10.0 ns (audited - 6 real sites):
| File | Site | Before | After |
|------|------|--------|-------|
| run_amber.py:94                 | Amber CLI argparse        | 10  | 10 (kept) |
| amber_md/mmgbsa_openmm.py:~329  | OpenMM fn default         | 5   | 10 |
| amber_md/mmgbsa_openmm.py:760   | OpenMM CLI argparse       | 5   | 10 |
| amber_md/batch.py:368           | batch CLI argparse        | 50  | 10 |
| gui/pages/0_Setup_and_Launch.py | _def_prod (field default) | 10/5 by engine | 10 (both) |
| gui/pages/0_Setup_and_Launch.py | submit fallbacks (538/588)| 5/10 | 10/10 |

Also: OpenMM docstring example (--prod-ns 5 -> 10); GUI help/comment text now
says "10 ns for both engines"; engine-conditional _is_amber_mm branch removed
for the default (still used elsewhere).

## Not changed
* equil_ns default stays 1.0 ns everywhere.
* The compute-estimate line params.get("prod_ns", 0) keeps 0 as a guard (only
  hit if the widget value is missing, which never happens) - intentionally not
  forced to 10.
* MMPBSA.py frame interval default (2) unchanged.

## Net effect
Whether a user runs Amber or OpenMM MM-GBSA, via GUI or CLI, single or batch,
production MD now defaults to 10 ns. Users can still override via the GUI field
or --prod-ns.

## Carried forward
final46 (driver hardening), v2.5.3 (submit.py activate fix), final47 (anti-storm
guard), final49 (report kit: exec bits + MM[PG]BSA filename tolerance +
standalone report), final50 (unified Generate-HTML-report button).

## Verified
* changed modules compile; full package audit shows no prod default != 10
  (other than the intentional estimate guard).
