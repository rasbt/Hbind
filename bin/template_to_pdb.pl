#!perl
#
# Usage: template_to_pdb.pl <template file> > <pdb file>
#
# Converts <template file> to PDB format with BValues indicating type:
#     A (acceptor) = 0
#     D (donor) = 50
#     H (hydrophobic) = 100
#     N (doneptor) = 25
# and Occupancies indicating key/nonkey:
#     Key = 1.00
#     Nonkey = 0.50

@ARGV == 1 || die "Usage: $0 <template file> > <pdb file>\n";

$fn = $ARGV[0];
open ( IN, $fn )
    or die "couldn't open file `$fn` for reading\n";

%BValHash = ("A" => 0,
	     "D" => 50,
	     "H" => 100,
	     "N" => 25 );
%OccHash = ("*" => 1.0,
	    "" => 0.5 );

$count = 1;
while ( <IN> )
{
    next if $_ =~ /\#/;
    @fields = split;
    printf "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  %4.2f%6.2f\n", 
    $count,
    $count,
    $fields[1],
    $fields[2],
    $fields[3],
    $OccHash{substr($fields[0],1,1)},
    $BValHash{substr($fields[0],0,1)};
    $count++;
}

close IN;

