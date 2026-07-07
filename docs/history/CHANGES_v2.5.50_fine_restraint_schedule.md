# v2.5.50 - finer restraint lambda schedule (opt-in)

Date: 2026-06-21

Fallback for clash-prone restraint windows (lig_12944901 lambda=0.150). Default grid jumps 0.0->0.15->0.30; fine grid inserts 0.05/0.10/0.20 (16->19 windows). Enable with use_fine_restraint_lambdas=True. v2.5.49 recovery escalation is tried first.
