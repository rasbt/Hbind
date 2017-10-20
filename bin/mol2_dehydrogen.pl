#!perl
#
# mol2dehydrogen -- remove hydrogens from a mol2 file, including necessary
#   changes to atom and bond counts
#
# 23-Jan-01: This does NOT handle multiple substructure moelcules correctly 
#   in that the SUBSTRUCTURE section is left unchanged. - PCS

@ARGV == 1 || die "Usage: $0 <mol2file>\n";
$fn = $ARGV[0];
open IN, $fn or die "Unable to open input: $fn\n";
$atomcount = $bondcount = 0;
$numheavy = $numheavybonds = 0;
$numH = $numHbonds = 0;

while (<IN>) {
    if ($_ =~ /\<TRIPOS\>MOLECULE/) { # We've reached the molecule tag
	$_ = <IN>; # next line is the molname
	$_ = <IN>; # next line is the atom count line
	s/^\s+//g;
	($numatoms, $numbonds, $numsubstrs) = split(/\s+/,$_);
#	print "$numatoms, $numbonds, $numsubstrs\n";
    }
    if ($_ =~ /\<TRIPOS\>ATOM/) {
	$atomflag = 1;
    }
    if ($_ =~ /\<TRIPOS\>BOND/) {
	$atomflag = 0;
	$bondflag = 1;
    }
    if ($_ =~ /\<TRIPOS\>SUBSTRUCTURE/) {
	$bondflag = 0;
    }
    if ($atomflag && ($_ !~ /ATOM/)) {
	s/^\s+//g;
	($atomnum, $atomname, $x, $y, $z, $orbital, $res, $resname, 
	 $potential) = split(/\s+/,$_);
#	print "$atomnum, $atomname, $x, $y, $z, $orbital, $res, $resname, $potential\n";
	$atomcount ++;
	if ($orbital =~ /H/) {
#	    print "$atomnum, $atomname, $x, $y, $z, $orbital, $res, $resname, $potential\n";
	    $Hnums[$atomnum] = 1;
	    $numH++;	    
	} else {
	    $heavyatoms[$numheavy] = $atomnum;
	    $atomnums[$atomnum] = $numheavy+1;
	    $atomnames[$atomnum] = $atomname;
	    $xs[$atomnum] = $x;
	    $ys[$atomnum] = $y;
	    $zs[$atomnum] = $z;
	    $orbitals[$atomnum] = $orbital;
	    $ress[$atomnum] = $res;
	    $resnames[$atomnum] = $resname;
	    $potentials[$atomnum] = $potential;	    
	    $numheavy++;
	}
    }
    if ($bondflag && ($_ !~ /BOND/)) {
	s/^\s+//g;
	($bondnum, $bond1, $bond2, $bondorder) = split(/\s+/,$_);
#	print "$bondnum, $bond1, $bond2, $bondorder\n";
	$bondcount ++;
	if ($Hnums[$bond1] == 1 || $Hnums[$bond2] == 1) {
#	    print "$bondnum, $bond1, $bond2, $bondorder\n";
	    $numHbonds++;
	} else {
	    $heavybonds[$numheavybonds] = $bondnum;
	    $bondnums[$bondnum] = $numheavybonds+1;
	    $bond1s[$bondnum] = $bond1;
	    $bond2s[$bondnum] = $bond2;
	    $bondorders[$bondnum] = $bondorder;
	    $numheavybonds++;
	}
    }
}
#print "$numatoms, $numbonds\n";
if ($atomcount != $numatoms) {
    print "ERROR: Number of atoms (count: $atomcount, entry: $numatoms)\n";
    exit;
}
if ($bondcount != $numbonds) {
    print "ERROR: Number of bonds (count: $bondcount, entry: $numbonds)\n";
    exit;
}
if (($atomcount-$numH) != $numheavy) {
    print "ERROR: Num H (Total-H: ".($atomcount-$numH).", Count: $numheavy)\n";
    exit;
}
if (($bondcount-$numHbonds) != $numheavybonds) {
    print "ERROR: Num H (Total-H: ".($bondcount-$numHbonds).
	", Count: $numheavybonds)\n";
    exit;
}

if (0) {
    for ($i=0; $i<$numheavy;$i++) {
	$atomnum = $heavyatoms[$i];
	print "$atomnum:  ";
	print $atomnums[$atomnum].", ".
	  $atomnames[$atomnum].", ".
	  $xs[$atomnum].", ".
	  $ys[$atomnum].", ".
	  $zs[$atomnum].", ".
	  $orbitals[$atomnum].", ".
	  $ress[$atomnum].", ".
	  $resnames[$atomnum].", ".
	  $potentials[$atomnum]."\n";
    }
    for ($i=0; $i<$numheavybonds;$i++) {
	$bondnum=$heavybonds[$i];
	print "$bondnum: ";
	print $bondnums[$bondnum].", ".
	  $bond1s[$bondnum].", ".
	  $bond2s[$bondnum].", ".
	  $bondorders[$bondnum]." :: ";
	print $bondnums[$bondnum].", ".
	  $atomnums[$bond1s[$bondnum]].", ".
	  $atomnums[$bond2s[$bondnum]].", ".
	  $bondorders[$bondnum]."\n";
    }
}
# Close and reopen the input file to do the reading with the correct output
close IN;
open IN, $fn;

while (<IN>) {
    if ($_ =~ /MOLECULE/) {
	print $_; # print the MOLECULE line
	$_ = <IN>;
	print $_; # print the molecule name line
	$_ = <IN>; # we are now at the atoms/bonds/substrs lines
	print "$numheavy   $numheavybonds   $numsubstrs\n";
	next;
	
    }
    if ($_ =~ /\<TRIPOS\>ATOM/) {
	print $_;
	$atomflag = 1;
    }
    if ($_ =~ /\<TRIPOS\>BOND/) { # we have reached the end of the atoms -- 
                                  # output the stored atoms
	$atomflag = 0;
	$bondflag = 1;
	for ($i=0; $i<$numheavy;$i++) {
	    $atomnum = $heavyatoms[$i];
	    printf("%5d %-8s %9.4f %9.4f %9.4f %-9s %2d %5s %7.4f ****\n",
		   $atomnums[$atomnum], 
		   $atomnames[$atomnum], 
		   $xs[$atomnum], 
		   $ys[$atomnum], 
		   $zs[$atomnum], 
		   $orbitals[$atomnum], 
		   $ress[$atomnum], 
		   $resnames[$atomnum], 
		   $potentials[$atomnum]);
	}
	print $_;
	next;
    }
    if ($_ =~ /\<TRIPOS\>/ && $bondflag) { # we have reached the end of the bonds -- 
	                        # output the stored bonds
	$bondflag = 0;
	for ($i=0; $i<$numheavybonds;$i++) {
	    $bondnum=$heavybonds[$i];
	    printf ("%5d %5d %5d %s\n",
		    $bondnums[$bondnum], 
		    $atomnums[$bond1s[$bondnum]], 
		    $atomnums[$bond2s[$bondnum]], 
		    $bondorders[$bondnum]);
	}
    }	
    if (!$atomflag && !$bondflag) { # print the rest of the mol2 file
	print $_;
    }
}
