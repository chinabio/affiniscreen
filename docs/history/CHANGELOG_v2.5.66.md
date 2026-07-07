# amber_md v2.5.66

Added tools/analyze_campaign.py: campaign-wide ABFE re-analyzer.
- Combines legs: dG_complex=c_decharge+c_vdw+c_restraint(+Boresch once); dG_solvent=s_decharge+s_vdw; dG_bind=-(dG_complex-dG_solvent)+dG_charge.
- Option-A restraint leg already carries Boresch in dG_kcal_mol -> NOT added twice.
- Resume-safe (reuse summary.json unless --force); --recurse for many edges.
- Same trust gating + exit codes as the cycle-closer.
