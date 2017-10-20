#!perl
#
# binding_site_waters.pl    Volker Schnecke  Wed Sep 15 15:26:34 EDT 1999
#
# Usage:  binding_site_waters.pl target template distance
#
# This script extracts all water molecules from the PDB file 
# $HBIND_DATA_DIR/<target>/<target>.pdb that are within <distance> 
# Angstrom from any of the template points specified in 
# $HBIND_DATA_DIR/<target>/<template>/in/. These waters are
# all HETATM entries with residue name HOH or WAT.
# The water molecules are then stored in the file 
# $HBIND_DATA_DIR/<target>/<template>/in/waters.pdb
# and read by HBIND during screening.

if ( $#ARGV != 2 )
{
    die "Usage: $0 <target> <template> <distance>\n";
}

$target = $ARGV[0];
$template = $ARGV[1];
$distance = $ARGV[2];

$TARGET_BASE = $ENV{HBIND_DATA_DIR}."/".$target."/";
$INPUT_DIRECTORY = $TARGET_BASE.$template."/in/";

$pdb_file = $TARGET_BASE.$target.".pdb";
$template_file = $INPUT_DIRECTORY."template";
$water_file = $INPUT_DIRECTORY."waters.pdb";

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

# now open the target PDB file and find all waters close to template points
open IN, $pdb_file;
open OUT, ">$water_file";
WHILE:
while ( <IN> )		       
{
    if ( /^HETATM.{11}HOH/ || /^HETATM.{11}WAT/ )
    {
	# all waters need as the residue specifier either HOH or WAT
	/.{30}(.{8})(.{8})(.{8})/;		
	$x = $1;
	$y = $2;
	$z = $3;
	for ( $i = 0; $i < $number_of_template_points; $i++ )
	{
	    if ( sqrt ( ( $x - $t_x[$i] ) * ( $x - $t_x[$i] )
		       + ( $y - $t_y[$i] ) * ( $y - $t_y[$i] )
		       + ( $z - $t_z[$i] ) * ( $z - $t_z[$i] ) ) < $distance )
	    {
		# this water is close to the current template point, so
		# write the atom line into the output file and read the
		# next water from the PDB file
		print OUT "$_";
		next WHILE;
	    }
	}
    }
}
close IN;
close OUT;
