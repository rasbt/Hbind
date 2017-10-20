#!perl
#
# Time-stamp: <Paul Sanschagrin -- Thu Jul 19 14:34:03 EDT 2001>
# merge_pdb_target.pl: replace the atoms coords in a base pdb file (i.e.
#   that used as a base for a hbind run) with the coords for side-chains
#   rotated during the run
# Usage: merge_pdb_target.pl [-n] <rotation file> <base PDB file>
#       -n:  no chain IDs used for the HBIND run (older versions)

use File::Basename;
use Getopt::Std;

sub usage {
  $this = basename ($0);
  print "Usage:    $this [-n] <rotation file> <base PDB file>\n".
    "Options:  -n:  no chain IDs used for the HBIND run ".
      "(older versions of HBIND\n";
    exit 1;
}


# Usage check and Parse CL args
if (getopts('n') == false) {
  usage();
}
usage() unless @ARGV == 2;


$rotfn = $ARGV[0];
$pdbfn = $ARGV[1];

# Process the rotation file

open ROT, $rotfn or die "Unable to open rotation file: $rotfn\n";
while (<ROT>) {
  next if $_ =~ /^REMARK/; # skip comment lines
  if ($opt_n) {
    # older verisons of  HBIND do _NOT_ output chain info for the target files
    $newcoords{substr($_,12,9)." ".substr($_,22,8)} = substr($_,30,24);
  } else {
    $newcoords{substr($_,12,18)} = substr($_,30,24);
  }
}

close ROT;

# Process the base PDB file
open PDB, $pdbfn or die "Unable to open base PDB file: $pdbfn\n";
while (<PDB>) {
  if ($opt_n) { # Pre-ChainID output versions
    $key = substr($_,12,9)." ".substr($_,22,8);
  } else { # new versions with ChainID (& AltLoc, InsCode)
    $key = substr($_,12,18);
  }
  if (exists ($newcoords{$key})) {
    print substr($_,0,30).$newcoords{$key}.substr($_,54,27);
  } else {
    print $_;
  }
}
