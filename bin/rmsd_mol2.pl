#!perl
                                                                                
if ( $#ARGV != 1 )
{
    die "Usage: $0 <Candidate mol2-file1> <Reference mol2-file2>\n";
}
                                                                                
sub dist
{
    my (  $a, @points ) = @_;
                                                                                
    $distance = sqrt ( ( $points[$a][0] - $points[$b][0] )
                        * ( $points[$a][0] - $points[$b][0] )
                      + ( $points[$a][1] - $points[$b][1] )
                        * ( $points[$a][1] - $points[$b][1] )
                      + ( $points[$a][2] - $points[$b][2] )
                         * ( $points[$a][2] - $points[$b][2] ) );
}
                                                                                
$count1 = $count = 0;
                                                                                
open IN, $ARGV[0] or die "Unable to open <mol2-file1>: $ARGV[0]\n";
while ( <IN> )
{
    last if /^@<TRIPOS>ATOM/;
}
while ( <IN> )
{
    last if /@<TRIPOS>BOND/;
    s/^\s+//;
    @line = split;
    next if $line[5] eq "H";
    $atoms1[$count1][0] = $line[2];
    $atoms1[$count1][1] = $line[3];
    $atoms1[$count1][2] = $line[4];
    $count1++;
}
close IN;
                                                                                
open IN, $ARGV[1] or die "Unable to open <mol2-file2>: $ARGV[1]\n";
while ( <IN> )
{
    last if /^@<TRIPOS>ATOM/;
}
while ( <IN> )
{
    last if /@<TRIPOS>BOND/;
    s/^\s+//;
    @line = split;
    next if $line[5] eq "H";
    $atoms2[$count2][0] = $line[2];
    $atoms2[$count2][1] = $line[3];
    $atoms2[$count2][2] = $line[4];
    $count2++;
}
close IN;
                                                                                
if ( $count1 != $count2 )
{
    die "$0: Different number of atoms in both files $ARGV[0]<->$ARGV[1] ($count1<->$count2) \n";
}
                                                                                
$err = 0.0;
for ( $i = 0; $i < $count1; $i++)
{
    $distance = sqrt ( ( $atoms1[$i][0] - $atoms2[$i][0] )
                        * ( $atoms1[$i][0] - $atoms2[$i][0] )
                      + ( $atoms1[$i][1] - $atoms2[$i][1] )
                        * ( $atoms1[$i][1] - $atoms2[$i][1] )
                      + ( $atoms1[$i][2] - $atoms2[$i][2] )
                         * ( $atoms1[$i][2] - $atoms2[$i][2] ) );
    #print "$i  $distance\n";
    $err += $distance * $distance;
}
$err /= $count1;
$err = sqrt ( $err );
                                                                                
$outname = $ARGV[0];
if ($outname =~ /\.mol2/) {
  $outname =~ s/\.mol2//g;
}
 print "$outname  $err\n";

