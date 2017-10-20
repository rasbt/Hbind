#!perl
#
# filter_ligands.pl        Volker Schnecke   Wed Sep 15 16:13:00 EDT 1999
#
# Usage: filter_ligands.pl target template db  cutoff
#
# This script checks HBIND's score for all ligands found in database
# <db> when screening for target <target> with template <template>.
# The command-line parameter <cutoff> defines the lower bound for
# the score, all ligands with the corresponding target and water files
# in the directories <db>_targets and <db>_waters, respectively,
# that have a lower score than <cutoff> are deleted.

if ( $#ARGV != 3 )
{
    die "Usage: $0 <target> <template> <db> <score_cutoff>\n";
}

$base_name = $ENV{HBIND_DATA_DIR}."/".$ARGV[0]."/".$ARGV[1]."/".$ARGV[2];

$ligand_dir = $base_name."_ligands";
$target_dir = $base_name."_targets";
$water_dir = $base_name."_waters";

opendir DIR, $ligand_dir;
@files = grep ( /.mol2$/, grep ( !/^.\.?$/, readdir ( DIR ) ) );
close DIR;

for ( $i = 0; $i <= $#files; $i++ )
{
    open IN, "$ligand_dir/$files[$i]";
    while ( <IN> )
    {
	if ( /^# affiscore/ )
	    {		
		chop;
		/\s([\-0-9.]*)$/;
		$score = $1;
		if ( $score > $ARGV[3] )
		{
		    print "deleting $files[$i] (score $score)\n";
		    unlink "$ligand_dir/$files[$i]";
		    $files[$i] =~ s/.mol2$/.pdb/;
		    unlink "$target_dir/$files[$i]";
		    unlink "$water_dir/$files[$i]";
		}
		last;
	    }
    }
    close IN;
}


