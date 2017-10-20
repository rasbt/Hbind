#!perl
#
# binding_site_residues.pl    Matt Tonero      Thu Feb 08, 2007
#                             Volker Schnecke  Wed Sep 15 15:04:03 EDT 1999
#
# Usage: binding_site_residues.pl target template distance
#
# This script reads the PDB-file <target>.pdb, which should be
# located in the base directory for this target protein (which is
# $HBIND_DATA_DIR/<target>/) and the corresponding template file, and 
# outputs a PDB file <target>.rad in the input directory for this 
# template (which is $HBIND_DATA_DIR/<target>/<template>/in/)
# that includes only those residues of the target protein that have 
# at least one atom within the distance <distance> from any template point.

if ( $#ARGV != 2 )
{
    die "Usage: $0 <target> <template> <distance>\n";
}

$target = $ARGV[0];
$template = $ARGV[1];
$distance = $ARGV[2];

$marker = 0;
$marker2 = 1;
$TARGET_BASE = $ENV{HBIND_DATA_DIR}."/".$target."/";
$TEMPLATE_BASE = $TARGET_BASE.$template."/";

$target_file = $TARGET_BASE.$target.".pdb";
$template_file = $TEMPLATE_BASE."in/"."template";
$binding_site = $TEMPLATE_BASE."in/".$target.".rad";

# get all template points
open IN, $template_file;
$number_of_template_points = 0;
while ( <IN> )
{
    if ( /^\#/ )
    {			
	# skip comment lines
	next;
    }
    @line = split;
    $t_x[$number_of_template_points] = $line[1];
    $t_y[$number_of_template_points] = $line[2];
    $t_z[$number_of_template_points] = $line[3];
    $number_of_template_points++;
}
close IN;

# now read through the target PDB file and compute the distance of its
# atoms to the template points
open IN, $target_file;
open OUT, ">$binding_site";

WHILE:
while ( <IN> )		       
{
    if ( /^ATOM/ || ( /HETATM/ && ! /^HETATM.{11}HOH/ && ! /^HETATM.{11}WAT/ ) )
    {
	$record_name = substr $_, 0, 6;  # extract the record name field
	$atom_num = substr $_, 6, 5;  # extract the atom number field
	$atom_name = substr $_, 12, 4;  # extract the atom name field
	$altLoc = substr $_, 16, 1;  # extract the alternate location field
	$resName = substr $_, 17, 3;  # extract the residue name field
	$chainID = substr $_, 21, 1;  # extract the chain ID field
	$resNum = substr $_, 22, 4;  # extract the residue number field
	if ($resNum < 0)
	{
	    $resNum = 1000000 + $resNum;
	}
	$iCode = substr $_, 26, 1;  # extract the code for insertion of residues field
	$X_coord = substr $_, 30, 8;  # extract the X coordinate field
	$Y_coord = substr $_, 38, 8;  # extract the Y coordinate field
	$Z_coord = substr $_, 46, 8;  # extract the Z coordinate field
	$occupancy = substr $_, 54, 6;  # extract the Occupancy field
	$tempFactor = substr $_, 60, 6;  # extract the Temperature factor field
	$segID = substr $_, 72, 4;  # extract the segment ID (left justified) field
	$element = substr $_, 76, 2;  # extract the element symbol (right justified) field
	$charge = substr $_, 78, 2;  # extract the charge on the atom field
	if (($resmark[$chainID][$resNum] == 1 ) && ($oldRes == $resNum) && ($oldChain == $chainID))
	{
#	    print "$chainID $resNum\n";
	    next;
	}
	for ( $i = 0; $i < $number_of_template_points; $i++ )
	{
	    $distance2 = ( sqrt ( ( $X_coord - $t_x[$i] ) * ( $X_coord - $t_x[$i] )
				  + ( $Y_coord - $t_y[$i] ) * ( $Y_coord - $t_y[$i] )
				  + ( $Z_coord - $t_z[$i] ) * ( $Z_coord - $t_z[$i] ) ));
	    if ($distance2 < $distance)
	    {
		$resmark[$chainID][$resNum] = 1;
		$oldRes = $resNum;
		$oldChain = $chainID;
		# reopen the target PDB file and output all residues that have at least
# one atom close to any template point				
		open IN2, $target_file;
	      WHILE:
		while ( <IN2> )		       
		{
		    if ( /ATOM/ || ( /HETATM/ && ! /^HETATM.{11}HOH/ && ! /^HETATM.{11}WAT/ ) )
		    {	
			$record_name2 = substr $_, 0, 6;  # extract the record name field
			$atom_num2 = substr $_, 6, 5;  # extract the atom number field
			$atom_name2 = substr $_, 12, 4;  # extract the atom name field
			$altLoc2 = substr $_, 16, 1;  # extract the alternate location field
			$resName2 = substr $_, 17, 3;  # extract the residue name field
			$chainID2 = substr $_, 21, 1;  # extract the chain ID field
			$resNum2 = substr $_, 22, 4;  # extract the residue number field
			if ($resNum2 < 0)
			{
			    $resNum2 = 1000000 + $resNum2;
			}
			$iCode2 = substr $_, 26, 1;  # extract the code for insertion of residues field
			$X_coord2 = substr $_, 30, 8;  # extract the X coordinate field
			$Y_coord2 = substr $_, 38, 8;  # extract the Y coordinate field
			$Z_coord2 = substr $_, 46, 8;  # extract the Z coordinate field
			$occupancy2 = substr $_, 54, 6;  # extract the Occupancy field
			$tempFactor2 = substr $_, 60, 6;  # extract the Temperature factor field
			$segID2 = substr $_, 72, 4;  # extract the segment ID (left justified) field
			$element2 = substr $_, 76, 2;  # extract the element symbol (right justified) field
			$charge2 = substr $_, 78, 2;  # extract the charge on the atom field
			# get the residue number
#	/.{22}(.{4})/;	
			if (($chainID2 eq $chainID) && ($resNum2 == $resNum))
			{
			    print OUT $_;
			}
		    }
		}
		last;
	    }

	}
    }
    
#    if ( /HETATM/ && ! /^HETATM.{11}HOH/ && ! /^HETATM.{11}WAT/ )
 #   {	
	# if the target PDB file includes heteroatoms, which are no waters,
	# we assume that they are located in the binding site and directly
	# copy the line after removing the chain ID
	#s/(.{21})./$1 /;
	#print OUT $_;
    #}
}


close IN;
close OUT;
close IN2;
