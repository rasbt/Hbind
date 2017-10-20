#!perl

if ( $#ARGV < 2 )
{
    die "Usage: $0 <max ligand-template distance> <template file> <pdb file> <ligand mol2 file #1> ... <ligand mol2 file #n>\n";
    }

$threshold = $ARGV[0];
$extended_threshold = $ARGV[0] * 2;
$template_extend_threshold = $threshold;
$buried_thresh = 4.5;
$neighbor_thresh = 10;
$pass = 1;
$fail = 0;

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

$count1 = $count2 = $count3 = 0;

open TEMPLATE, $ARGV[1] or die "Unable to open <template file>: $ARGV[1]\n";
while ( <TEMPLATE> )
{
    s/^\s+//;
    @line = split;
    if ($line[1] == '#')
    {
	printf $_;
    }
    $atoms1[$count1][0] = $line[1];
    $atoms1[$count1][1] = $line[2];
    $atoms1[$count1][2] = $line[3];
    $atoms1[$count1][3] = $line[0];
    $count1++;
}
close TEMPLATE;

open PDB, $ARGV[2] or die "Unable to open <pdb-file>: $ARGV[2]\n";
#    print "opening $lig_done $ARGV[$lig_done]\n";
while ( <PDB> )
{
    $record_name = substr $_, 0, 6;  # extract the record name field
    $atom_num = substr $_, 6, 5;  # extract the atom number field
    $atom_name = substr $_, 12, 4;  # extract the atom name field
    $altLoc = substr $_, 16, 1;  # extract the alternate location field
    $resName = substr $_, 17, 3;  # extract the residue name field
    $chainID = substr $_, 21, 1;  # extract the chain ID field
    $resNum = substr $_, 22, 4;  # extract the residue number field
    $iCode = substr $_, 26, 1;  # extract the code for insertion of residues field
    $X_coord = substr $_, 30, 8;  # extract the X coordinate field
    $Y_coord = substr $_, 38, 8;  # extract the Y coordinate field
    $Z_coord = substr $_, 46, 8;  # extract the Z coordinate field
    $occupancy = substr $_, 54, 6;  # extract the Occupancy field
    $tempFactor = substr $_, 60, 6;  # extract the Temperature factor field
    $segID = substr $_, 72, 4;  # extract the segment ID (left justified) field
    $element = substr $_, 76, 2;  # extract the element symbol (right justified) field
    $charge = substr $_, 78, 2;  # extract the charge on the atom field
    
    if (($record_name eq 'ATOM  ') || ($record_name eq 'HETATM'))
    {
	$atoms3[$count3][0] = $X_coord;
	$atoms3[$count3][1] = $Y_coord;
	$atoms3[$count3][2] = $Z_coord;
	$atoms3[$count3][3] = $atom_name;
	$count3++;
    }
}
close PDB;

$lig_done = 3;
while ($lig_done <=  $#ARGV)
{
    open LIGAND, $ARGV[$lig_done] or die "Unable to open <mol2-file>: $ARGV[$lig_done]\n";
#    print "opening $lig_done $ARGV[$lig_done]\n";
    while ( <LIGAND> )
    {
	last if /^@<TRIPOS>ATOM/;
    }
    while ( <LIGAND> )
    {
	last if /@<TRIPOS>BOND/;
	s/^\s+//;
	@line = split;
	next if ($line[5] eq "H") || ($line[5] eq "H_ADD");
	$atoms2[$count2][0] = $line[2];
	$atoms2[$count2][1] = $line[3];
	$atoms2[$count2][2] = $line[4];
	$atoms2[$count2][3] = $line[0];
	$count2++;
    }
    close LIGAND;
    $lig_done++;
}

$err = 0.0;
$max = 0.0;
$min = 100.0;
$closest = 9999.0;
$k = $l = 0;
$count_new = $count_far = 0;

# compare ligand-template distance and keep templates within <threshold> distance of ligand atoms (label new_temp).  Label the remaining template points far_temp 
for ( $i = 0; $i < $count1; $i++)
{
    $closest = 9999.0;
    for ( $j = 0; $j < $count2; $j++)
    {
	$distance = sqrt ( ( $atoms1[$i][0] - $atoms2[$j][0] )
			   * ( $atoms1[$i][0] - $atoms2[$j][0] )
			   + ( $atoms1[$i][1] - $atoms2[$j][1] )
			   * ( $atoms1[$i][1] - $atoms2[$j][1] )
			   + ( $atoms1[$i][2] - $atoms2[$j][2] )
			   * ( $atoms1[$i][2] - $atoms2[$j][2] ) );
	#print "$i: $distance\n";
	#$err += $distance * $distance;  
#	printf ("%d %s -- %d %s : %0.3f %0.3f\n", $i+1, $atoms1[$i][3], $j+1, $atoms2[$j][3], $distance, $closest);
	
	if ($distance < $closest){
	    $closest = $distance;}
	if ($distance > $max){
	    $max = $distance;}
	if ($distance < $min){
	    $min = $distance;}
    }
    if ($closest <= $threshold)
    {
	$new_temp[$k][0] = $atoms1[$i][0];
	$new_temp[$k][1] = $atoms1[$i][1];
	$new_temp[$k][2] = $atoms1[$i][2];
	$new_temp[$k][3] = $atoms1[$i][3];
	printf("%s %0.3f %0.3f %0.3f\n", $atoms1[$i][3], $atoms1[$i][0], $atoms1[$i][1], $atoms1[$i][2]);
	$k++;
	$count_new++;
    }
    if (($closest > $threshold) && ($closest <= $extended_threshold))
    {
	$far_temp[$l][0] = $atoms1[$i][0];
	$far_temp[$l][1] = $atoms1[$i][1];
	$far_temp[$l][2] = $atoms1[$i][2];
	$far_temp[$l][3] = $atoms1[$i][3];
#	printf("%s %0.3f %0.3f %0.3f\n", $atoms1[$i][3], $atoms1[$i][0], $atoms1[$i][1], $atoms1[$i][2]);
	$l++;
	$count_far++;
    }
}

#this part of the code finds pockets near the known ligands in the active site that may be possible to exploit in ligand docking

$added = $count_new;
$k = 0;
while ($added != 0) #stop iterating when no new template points found
{
    $added = 0;
    for ( $i = 0; $i < $count_far; $i++)
    {
	$neighbors[$i] = 0;
	for ( $j = 0; $j < $count3; $j++)
	{
	    $distance = sqrt ( ( $far_temp[$i][0] - $atoms3[$j][0] )
			       * ( $far_temp[$i][0] - $atoms3[$j][0] )
			       + ( $far_temp[$i][1] - $atoms3[$j][1] )
			       * ( $far_temp[$i][1] - $atoms3[$j][1] )
			       + ( $far_temp[$i][2] - $atoms3[$j][2] )
			       * ( $far_temp[$i][2] - $atoms3[$j][2] ) );
	    if ($distance <= $buried_thresh)
	    {
		$neighbors[$i]++; # find how many protein atoms are within <buried_thresh> of the far_temp points
	    }
	}
	if ($neighbors[$i] >= $neighbor_thresh)
	{
	    $smallest_gap = 9999.0;
	    $status = $fail;
#	    print "cn $count_new\n";
	    for ( $l=0; $l < $count_new; $l++)
	    {   # for far_temp points with enough close protein points, find distance between the far_temp and the new_temp points
		$distance_a = sqrt ( ( $far_temp[$i][0] - $new_temp[$l][0] )
				     * ( $far_temp[$i][0] - $new_temp[$l][0] )
				     + ( $far_temp[$i][1] - $new_temp[$l][1] )
				     * ( $far_temp[$i][1] - $new_temp[$l][1] )
				     + ( $far_temp[$i][2] - $new_temp[$l][2] )
				     * ( $far_temp[$i][2] - $new_temp[$l][2] ) );
#		print "a= $distancs_a\n";
		$smallest_gap = (($distance_a*$distance_a)/2.0 + 2.0*(1.5 * 1.5)); # this ensures that the closest pdb atom between the new and far template points is at least 1.5 A perpendicular to the line connecting new_temp and far_temp points (basically that the 2 points have a "line of sight" to each other)
#		print "a= $distance_a, smallest_gap = $smallest_gap\n";
#		print "status = $status\n";
		if ($distance_a <= $template_extend_threshold)
		{
		    for ( $m=0; $m < $count3; $m++)
		    {
			$distance_b = sqrt ( ( $far_temp[$i][0] - $atoms3[$m][0] )
					     * ( $far_temp[$i][0] - $atoms3[$m][0] )
					     + ( $far_temp[$i][1] - $atoms3[$m][1] )
					     * ( $far_temp[$i][1] - $atoms3[$m][1] )
					     + ( $far_temp[$i][2] - $atoms3[$m][2] )
					     * ( $far_temp[$i][2] - $atoms3[$m][2] ) );
			#	print "b= $distancs_b\n";
			$distance_c = sqrt ( ( $new_temp[$l][0] - $atoms3[$m][0] )
					     * ( $new_temp[$l][0] - $atoms3[$m][0] )
					     + ( $new_temp[$l][1] - $atoms3[$m][1] )
					     * ( $new_temp[$l][1] - $atoms3[$m][1] )
					     + ( $new_temp[$l][2] - $atoms3[$m][2] )
					     * ( $new_temp[$l][2] - $atoms3[$m][2] ) );
			#print "c= $distancs_c\n";
			if ($distance_b * $distance_b + $distance_c * $distance_c >= $smallest_gap) # see prior comment about smallest_gap
			{
#			print "$distance_a\n";
			    if ($distance_a <= $template_extend_threshold) # this makes sure that only far_temp points that are within reasonable distance are accepted as new_temp points
			    {
				$status = $pass;
				next;
			    }
			}
		    }
		}
	    }
	    if ($status == $pass)
	    {
		$count_new++;
		$new_temp[$count_new][0] = $far_temp[$i][0];
		$new_temp[$count_new][1] = $far_temp[$i][1];
		$new_temp[$count_new][2] = $far_temp[$i][2];
		$new_temp[$count_new][3] = $far_temp[$i][3];
		$far_temp[$i][0]*=10000;
		$far_temp[$i][1]*=10000;
		$far_temp[$i][2]*=10000;
		printf("%s %0.3f %0.3f %0.3f\n", $new_temp[$count_new][3], $new_temp[$count_new][0], $new_temp[$count_new][1], $new_temp[$count_new][2]);
		$k++;
#		$count_new++;
#		printf("new %s %0.3f %0.3f %0.3f\n", $far_temp[$i][3], $far_temp[$i][0], $far_temp[$i][1], $far_temp[$i][2]);
		$added++;
	    }
	}
	
#    if ($neighbors[$i] >= $neighbor_thresh)
#    {
#	$buried_temp[$k][0] = $far_temp[$i][0];
#	$buried_temp[$k][1] = $far_temp[$i][1];
#	$buried_temp[$k][2] = $far_temp[$i][2];
#	$buried_temp[$k][3] = $far_temp[$i][3];
#	printf("%s %0.3f %0.3f %0.3f\n", $far_temp[$i][3], $far_temp[$i][0], $far_temp[$i][1], $far_temp[$i][2]);
#	$k++;
#	$count_buried++;
#    }
    }
}

for ( $i = 0; $i < $count1; $i++)
{
    $neighbors[$i]=0;
    $closest = 9999.0;
    for ( $j = 0; $j < $count3; $j++)
    {
	$distance = sqrt ( ( $atoms1[$i][0] - $atoms3[$j][0] )
			   * ( $atoms1[$i][0] - $atoms3[$j][0] )
			   + ( $atoms1[$i][1] - $atoms3[$j][1] )
			   * ( $atoms1[$i][1] - $atoms3[$j][1] )
			   + ( $atoms1[$i][2] - $atoms3[$j][2] )
			   * ( $atoms1[$i][2] - $atoms3[$j][2] ) );
	if ($distance <= $buried_thresh)
	{
	    $neighbors[$i]++;
	}
    }
#3    printf("template point %d has %d pdb neighbors within %.2f Angstroms.\n",$i+1, $neighbors[$i], $buried_thresh );
}

