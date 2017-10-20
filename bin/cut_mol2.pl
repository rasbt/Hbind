#!perl
#
# cut_mol2.pl
# Time-stamp: <Paul Sanschagrin -- Fri May 18 09:46:45 EDT 2001>
# Usage: cut_mol2.pl file
#
# This script cuts a mol2 file into into separate files containing 
# individual molecules. The files are named 'MOLNAME_#.mol2'.


if ( $#ARGV != 0 )
{
    die "Usage: $0 <mol2-file>\n";
}

$num_molecules = 0;
open IN, "$ARGV[0]";
while ( <IN> )
{
    if ( /@<TRIPOS>MOLECULE/ )
    {
	$_ = <IN>;
	chop;
	if ($num_molecules < 10) {
	  $numstring = "00".$num_molecules;
	} elsif ($num_molecules < 100) {
	  $numstring = "0".$num_molecules;
	} else {
	  $numstring = $num_molecules;
	}
	$filename = $_."_".$numstring;
	$num_molecules++;
	# we do not explicitely have to close OUT, since close is
	# automatcally done when we open an existing filehandle again
	open OUT, ">$filename.mol2";
	print OUT "@<TRIPOS>MOLECULE\n";
	print OUT "$filename\n";
    }
    else
    {
	print OUT "$_";
    }
}
close OUT;
close IN;



