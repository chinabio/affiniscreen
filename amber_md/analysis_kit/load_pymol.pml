# load_pymol.pml -- auto-generated
# Loads the solvated complex + trajectory, hides water, highlights ligand.
# Usage:  pymol load_pymol.pml
#         (or)  open with PyMOL

# --- Paths are RELATIVE to the analysis/ folder ------------------------
load    ../complex_solv.prmtop, complex
load_traj ../md/prod.nc, complex, state=1, interval=10   # every 10th frame for speed

# --- Cosmetic ----------------------------------------------------------
hide everything
bg_color white
show cartoon, polymer.protein
color grey80, polymer.protein

# Ligand: try residue name LIG first, then any non-polymer hetero
select lig, resn LIG or (hetatm and not resn HOH+WAT+NA+CL+K+MG+ZN+CA)
show sticks, lig
color magenta, lig
util.cnc lig

# Hide water + ions
hide everything, resn HOH+WAT+NA+CL+K+MG+ZN
remove resn HOH+WAT          # comment this out if you want to *see* waters

# Binding-site pocket: residues within 5 Å of ligand, lines only
select pocket, byres (polymer.protein within 5 of lig)
show lines, pocket
color cyan, pocket and elem C

zoom lig, 8
set ray_shadow, 0
set cartoon_transparency, 0.15
mplay
