#!perl
# Usage: pdb_to_template.pl <pdb file> > <template file>
#
# This script reads the PDB-file, which has water molecules
# at the positions of template points and outputs a template-file,
# the types of the template points in that file are based on the 
# temperature factor listed for each water and the key designation
# (key or nonkey) is based on the occupancy

@ARGV == 1 || die "Usage: <pdb file> > <template file>\n";

$fn = $ARGV[0];
open ( IN, $fn )
    or die "couldn't open file `$fn` for reading";

%TypeHash = (100 => "H", # Hydrophobic point
	     50 => "D",  # Donor point
	     0 => "A",  # Acceptor point
	     25 => "N");  # Doneptor point

%KeyHash = (1.00 => "*",  # '*' indicates a key point (occ = 1.0)
	    0.50 => " "); # ' ' indicates a nonkey point (occ = 0.5)

$count = 0;
while ( <IN> )
{
    /.{30}(.{8})(.{8})(.{8})(.{6})(.{6})/;
    printf "%1s%1s  %8.3f %8.3f %8.3f\n", 
    $TypeHash{$5+0},
    $KeyHash{$4+0},
    $1, $2, $3;


}

close IN;
