#!perl

# Usage: check_number_of_atoms_mol2.pl <directory> <min_atoms> <max_atoms>
#
# Volker Schnecke   Tue Feb 24 18:09:00 EST 1998
#
# This script checks each mol2-files in a directory and checks its
# number of atoms. If there are less than a minimal number or more
# than a maximal number of atoms, the file is deleted.
# 20-Dec-2001 PCS: Modified to examine only heavy atoms

if ( $#ARGV != 2 )
{
  die "Usage: $0 <directory> <min_atoms> <max_atoms>\n";
}

$dirname = $ARGV[0];
$min_atoms = $ARGV[1];
$max_atoms = $ARGV[2];

opendir DIR, $dirname;
@files = grep ( /.mol2$/, grep ( !/^.\.?$/, readdir ( DIR ) ) );
close DIR;

foreach $file ( @files )
  {
    next unless open IN, "$dirname/$file";
    print "$file: ";
    $numheavy = $numhydrogen = 0;
    while (<IN>) {
      chop;
      if ($_ =~ /\<TRIPOS\>ATOM/) { # we've reached the atom section
	$atomflag = 1;
      }
      if (($_ =~ /\<TRIPOS\>/) && ($_ !~ /ATOM/)) {
	# we've reached the next section
	$atomflag = 0;
      }
      if ($atomflag && ($_ !~ /ATOM/)) { # read the atoms and check for H's
	s/^\s+//g;
	($atomnum, $atomname, $x, $y, $z, $orbital, $res, $resname, 
	 $potential) = split(/\s+/,$_);
	if ($orbital ne "H") {
	  $numheavy++;
	} else {
	  $numhydrogen++;
	}
      }
    }
    $numtotal = $numheavy + $numhydrogen;
    print "nonH: $numheavy, H: $numhydrogen, Total: $numtotal\n";
    close IN;
    if ( $numheavy < $min_atoms || $numheavy > $max_atoms )
      {
	unlink $dirname."/".$file;
	print "deleting $file ($numheavy nonH atoms)\n";
      }
  }
