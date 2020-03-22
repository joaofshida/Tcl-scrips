# Original script by Joao F. Shida
# 21 march 2020 - Coronavirus Quarentine
#
# Measures the average diameter of a single beta barrel from the residues that belonged
# to a beta barrel for a T percentage of the entire trajectory. Barrel and or relevant structure must be aligned to the
# z-axis, from which the radius will be calculated.
# Usage: source Bbarrel-radius.tcl
# barrelRadius T outfile max_resid
# Outputs the average radius of the beta barrel every trajectory step

proc barrelRadius {T outfile max_resid} {
  set outfile [open $outfile w]
	set mol top
  set nSteps [molinfo $mol get numframes]
  puts $nSteps
  # Initialize list of beta occurence
  set sel [atomselect top "name CA"]
  set struct [$sel get structure]
  set len [llength $struct]
  set betaList {}
  for {set i 0} {$i < $len} {incr i} {
      lappend betaList 0
  }

  # For every Frame
  for {set i 0} {$i < $nSteps} {incr i} {
    # Get structure list for this frame
    set sel [atomselect top "name CA"]
    $sel frame $i
    set thisstruc [$sel get structure]
    # puts $thisstruc
    set len [llength $thisstruc]
    # For every residue, check if beta sheet. if so, then add to the betaList
    for {set j 0} {$j < $len} {incr j} {
      set struct [lindex $thisstruc $j]
      if {$struct == "E"} {
      set previousCount [lindex $betaList $j]
      lset betaList $j [expr $previousCount + (double(1)/$nSteps)]
      } else {
      }
    }
  }

  set betaIndexList {}
  set j 0
  for {set i 1} {$i <= $max_resid} {incr i} {
    set originali $i
    set residueSel [atomselect top "resid $i"]
    set atoms [$residueSel num]
    # puts "$atoms in resid $i"
    while {$atoms == 0} {
        set i [expr $i + 1]
        set residueSel [atomselect top "resid $i"]
        set atoms [$residueSel num]
	if { $i > $max_resid} {
    	    break
	    # puts "Something Went Wrong"
  	}
	} 
    # puts "$i $originali" 
    set ammountOfBeta [lindex $betaList $j]
    if {$ammountOfBeta >= $T} {
        lappend betaIndexList $i 
        # puts $ammountOfBeta
    } else {}
    #puts "$i"
    set j [expr $j + 1]
  }
  # puts $betaIndexList
  set lengthOfResidues [llength $betaIndexList]
  # puts $lengthOfResidues 
  # For every Frame
  for {set i 0} {$i < $nSteps} {incr i} {
    set frameSum 0
    for {set j 0} {$j < $lengthOfResidues} {incr j} {
        set resNum [expr [lindex $betaIndexList $j] + 1]
        set residueSel [atomselect top "resid $resNum"]
        $residueSel frame $i
	# puts $resNum
        set pos [measure center $residueSel weight mass]
        set x [lindex $pos 0]
        set y [lindex $pos 1]
        set x2 [expr $x*$x]
        set y2 [expr $y*$y]
        # Assumption: z axis is central axis of barrel
        set dOverSteps [expr { double(sqrt(double($x2 + $y2)))/$lengthOfResidues}]
        # puts $dOverSteps 
        set frameSum [expr $frameSum + $dOverSteps]
    }
  puts $outfile "$i\t$frameSum"
  }

	close $outfile
}
