# v2.5.41 - prod recovery wrapper Option A hotfix

Date: 2026-06-20

min/heat/dens/eq PASS but no prod.out on any window: run_prod_with_restart() aborted (return 74) because the Option A prod.in has no clambda. Fix: clambda guard now icfe=1-only; recovery mdin emits TI keywords only when icfe=1.
