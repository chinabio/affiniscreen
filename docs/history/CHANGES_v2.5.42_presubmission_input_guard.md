# v2.5.42 - pre-submission restraint-input guard

Date: 2026-06-20

After three Option A startup bugs (heat mdin, prod recovery wrapper) each from a place that still assumed TI, added a build-time guard. FEPSetup._validate_restraint_inputs() scans lambda_*/*.in in the restraint leg for icfe=1/timask/scmask/clambda/ifmbar/mbar_lambda/mbar_states/crgmask and raises RuntimeError listing offenders; called from setup_leg when stage=='restraint'. Standalone CLI tools/validate_restraint_leg_inputs.py mirrors the check (exit 0/1/2).
