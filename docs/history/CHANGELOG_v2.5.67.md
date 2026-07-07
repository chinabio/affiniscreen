# amber_md v2.5.67

Added CSV export to tools/analyze_campaign.py.
- --csv <path> writes one row per edge; with --recurse and no --csv, writes campaign_summary.csv at the campaign root.
- Rows ranked trusted-first then ascending dG_bind (tightest binder first).
- Columns: edge, edge_dir, dG_bind, trusted, complex/solvent totals, each leg dG, Boresch, and incomplete/unreliable/mbar_failed leg lists.
