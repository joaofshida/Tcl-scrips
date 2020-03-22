# Original script by Joao F. Shida
# 21 march 2020 - Coronavirus Quarentine
#
# Measures the average diameter of a single beta barrel from the residues that belonged
# to a beta barrel for a T percentage of the entire trajectory. Barrel and or relevant structure must be aligned to the
# z-axis, from which the radius will be calculated.
# Usage: source Bbarrel-radius.tcl
# barrelRadius T outfile
# Outputs the average radius of the beta barrel every trajectory step

proc barrelRadius {T outfile} {
  set outfile [open $outfile w]
	set mol top
  set nSteps [molinfo $mol get numframes]

  # Initialize list of beta occurence
  set sel [atomselect top "name CA"]
  set struct [$sel get structure]
  set len [llength $struct]
  set betaList {}
  for {set i 0} {$i <= $len} {incr i} {
      lappend betaList 0
  }

  # For every Frame
  for {set i 0} {$i <= $nSteps} {incr i} {
    # Get structure list for this frame
    animate top goto $i
    mol ssrecalc top
    set thisstruc [$sel get structure]
    set len [llength $thisstruc]
    # For every residue, check if beta sheet. if so, then add to the betaList
    for {set i 0} {$i <= $len} {incr i} {
      set struct [lindex $thisstruc $i]
      if {$struct == "B"} {
      set previousCount [llindex betaList $i]
      lset betaList $i [expr $previousCount + (1/$nSteps)]
      } else {
      }
    }


  set betaIndexList {}
  for {set i 0} {$i <= $nSteps} {incr i} {
    set ammountOfBeta [lindex betaList $i]
    if {$ammountOfBeta >= %T}{
        lappend betaIndexList $i
    } else {}
  }

  set lengthOfResidues [llength betaIndexList]
  # For every Frame
  for {set i 0} {$i <= $nSteps} {incr i} {
    set frameSum 0
    for {set j 0} {$j <= lengthOfResidues} {incr j} {
        set resNum [lindex $betaIndexList $j]
        set residueSel [atomselect top "backbone and resid $resNum"]
        set pos [center $residueSel]
        set x [lindex $pos 0]
        set y [lindex $pos 1]
        set x2 [expr $x*$x]
        set y2 [expr $y*$y]
        # Assumption: z axis is central axis of barrel
        set dOverSteps [expr { sqrt($x2 + $y2)/$nSteps }]
        set frameSum [expr $frameSum + $dOverSteps]
    puts $outfile "$i\t$frameSum"
    }
  }
  }
	close $outfile
}
