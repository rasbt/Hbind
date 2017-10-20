#!perl
#
# renum_pdb_atoms.pl       Volker Schnecke    Wed Sep 15 16:18:39 EDT 1999
#
# Usage: renum_pdb_atoms pdb_file type
#
# This script renumbers the atoms in the PDB-file <pdb_file> and is
# especially used for preparing the binding site and ligand file for the
# web interface.  The command-line parameter <type> has to be one of the
# keywords "TARGET" or "LIGAND", when renumbering the atoms for the
# binding site or the ligand, respectively.  The main purpose is to
# assign a chain ID to these files so that they can be easily
# distinguished in Rasmol when using the web interface.

if ( $#ARGV != 1 )
{
    die "Usage: $0 <pdb file> <LIGAND|TARGET>\n";
}

$number = 0;
open IN, $ARGV[0];
while ( <IN> )
{
    if ( /^ATOM/ || /^HETATM/ )	
    {
	# assign a new atom number
	substr ( $_, 7, 4 ) = sprintf "%4d", $number;
	if ( $ARGV[1] eq "LIGAND" )
	{
	    # if we are renumbering a ligand file, then set the chain ID to "K"
	    # and the temperature factor to 100.00
	    substr ( $_, 21, 1 ) = "K";
	    substr ( $_, 60, 6 ) = "100.00";
	}
	else
	{
	    # if we are renumbering a binding site file, then set chain ID to
	    # "T" and temperature factor to 0.00
	    substr ( $_, 21, 1 ) = "T";
	    substr ( $_, 60, 6 ) = "  0.00";
	}
	$number++;
    }
    print $_;
}
close IN;
