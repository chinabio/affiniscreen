# v2.5.48 - FEPAnalyzer._window_complete hotfix

Date: 2026-06-21

Run abfe_20260620_222307 lig_12944901: prod COMPLETED on windows (v2.5.44 fix works; recovery.log: 'prod completed OK (attempt 1)', prod.out has Final Performance Info). But the analyze job failed:
  AttributeError: 'FEPAnalyzer' object has no attribute '_window_complete'  (fep.py run() line ~2154)

Root cause: _window_complete() and _PROD_DONE_MARKERS were defined on FEPSetup (class @200), but FEPAnalyzer (class @1821) calls self._window_complete() in _collect_dvdl() and run(). Latent since the method was added -- only reachable now that prod.out finally exists (pre-2.5.44 the analyzer bailed at 'no dHdl' before the call).

Fix: define _PROD_DONE_MARKERS and _window_complete on FEPAnalyzer (self-contained copy).
