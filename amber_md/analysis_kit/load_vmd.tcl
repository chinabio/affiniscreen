# load_vmd.tcl -- auto-generated
# Usage:  vmd -e load_vmd.tcl
#         (or, inside VMD console:  source load_vmd.tcl)

# --- Load topology + trajectory ----------------------------------------
set top  "../complex_solv.prmtop"
set traj "../md/prod.nc"

mol new     $top  type parm7  waitfor all
mol addfile $traj type netcdf waitfor all step 1

# --- Representations ---------------------------------------------------
mol delrep 0 top
# protein cartoon
mol representation NewCartoon
mol selection      "protein"
mol color          ColorID 6
mol addrep         top
# ligand sticks (resname LIG, else any non-water hetero)
mol representation Licorice 0.25 12 12
mol selection      "resname LIG or (not protein and not water and not ion)"
mol color          Name
mol addrep         top
# binding-site residues within 5 Å
mol representation Lines
mol selection      "protein and (same residue as (within 5 of (resname LIG)))"
mol color          ColorID 4
mol addrep         top

# --- View setup --------------------------------------------------------
display projection Orthographic
display depthcue   off
color   Display    Background white
axes    location   Off

# Center on ligand
set sel [atomselect top "resname LIG"]
if {[$sel num] > 0} { set c [measure center $sel]; molinfo top set center [list $c] }
$sel delete

animate goto 0
