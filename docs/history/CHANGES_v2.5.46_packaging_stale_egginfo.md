# v2.5.46 - remove stale egg-info / pyc (fixes 'installed 2.5.31')

Date: 2026-06-21

Symptom: `pip install -e .` for v2.5.45 ended with 'Successfully installed amber_md-2.5.31.post7'.

Root cause: the source tree carried a committed amber_md.egg-info/ (PKG-INFO Version: 2.5.31.post8) and 66 __pycache__/*.pyc files, copied forward since 2.5.31. setuptools' editable backend reused the cached egg-info metadata instead of re-evaluating the dynamic version (attr=amber_md.__pep440_version__), so pip recorded the OLD version. The actual code/__version__ was correct (2.5.45).

Fix: deleted amber_md.egg-info/ and all __pycache__/. Added .gitignore and MANIFEST.in to keep them out. Extended tools/check_version_sync.py to fail if egg-info/PKG-INFO drifts or any .pyc is committed.

User action on the login node: after unzip, ensure a clean install:
  rm -rf amber_md.egg-info **/__pycache__
  pip install -e . --no-build-isolation --force-reinstall
  python -c 'import amber_md; print(amber_md.__version__)'   # must print 2.5.46
