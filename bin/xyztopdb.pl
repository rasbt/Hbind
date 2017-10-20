#!perl
#
# Usage: xyztopdb <XYZ file> > <pdb file>
#
# Converts <XYZ file> to PDB format of waters

@ARGV == 1 || die "Usage: $0 <coordinate file> ".
    "(use '-' for stdin)\n";

if ($ARGV[0] eq "-") {
    open (IN, "-");
} else {
    open (IN, $ARGV[0]) or die "Unable to open data file ".$ARGV[0]."\n";
}
$count = 1;
while ( <IN> )
{
    s/^\s+//g;
    @fields = split;
    printf "HETATM %4d  O   XXX %5d     %7.3f %7.3f %7.3f  %4.2f%6.2f\n",
    $count,
    $count,
    $fields[0],
    $fields[1],
    $fields[2],
    1.0,
    50.0;
    $count++;
}

close IN;
