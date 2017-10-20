#!perl
#
# show_ligands.pl        Volker Schnecke   Wed Sep 15 15:59:56 EDT 1999
#
# Usage: show_ligands.pl target template db
#
# This script generates a Rasmol script that can be used to browse
# through the ligands that HBIND found in database <db> for protein 
# <target> using binding-site template <template>.
# The script reads in all corresponding mol2 files based on their ranking,
# shows the ligands in rasmol, and outputs the score and other relevant
# information in the terminal window in which Rasmol was started.

# This variable has to point to the directory including the *.data files
# which lists name, formula, and publication information for all entries
# in the screening database. These files were created using the script
# 'create_csd_set_data.pl'
$CSD_SET_DATA_DIR = "/psa/db/csd-select/names";

if ( $#ARGV != 2 )
{
    die "Usage: $0 <target> <template> <db>\n";
}

sub sort_by_best_orientscore {
    $best_scores{$a} <=> $best_scores{$b};
}

sub sort_by_best_affiscore {
    $best_affiscore{$a} <=> $best_affiscore{$b};
}
                                                                                
$ligand_dir = $ENV{HBIND_DATA_DIR}."/".$ARGV[0]."/".$ARGV[1]."/".$ARGV[2]."_ligands";

opendir DIR, $ligand_dir;
@files = grep ( /.mol2$/, grep ( !/^.\.?$/, readdir ( DIR ) ) );
close DIR;

for ( $i = 0; $i <= $#files; $i++ )
{
    $base[$i] = $files[$i];
    $best_score{$base[$i]} = 0.0;
    $binding_modes{$base[$i]} = 0;
}
                                                                               
for ( $i = 0; $i <= $#files; $i++ )
{
    open IN, "$ligand_dir/$files[$i]";
    while ( <IN> )
    {
        if ( /^# Affinity score/ )
        {
            chop;
            /\s([\-0-9.]*)$/;
            $affiscore = $1;
        }

        if ( /^# buried protein hydrophobic term/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $buried_hphob = $1;
        }
        if ( /^# hydrophobic complementarity term/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $hphob_comp = $1;
        }
        if ( /^# polar component/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $polar_component = $1;
        }
        if ( /^# affinity constant/ )
        {
            chop;
            /\s([\-0-9.]*)$/;
            $affinity_const = $1;
        }
        if ( /^# protein-ligand hydrophobic contacts/ )
        {
            chop;
            /\s([0-9]*)$/;
            $pl_hphob = $1;
        }
        if ( /^# protein-ligand H-bond count/ )
        {
            chop;
            /\s([0-9]*)$/;
            $pl_hbond = $1;
        }
        if ( /^# protein-ligand salt-bridge count/ )
        {
            chop;
            /\s([0-9]*)$/;
            $pl_saltbridge = $1;
        }
        if ( /^# metal-ligand interactions count/ )
        {
            chop;
            /\s([0-9]*)$/;
            $ml_interactions = $1;
        }
        if ( /^# unsatisfied interfacial polar atom count/ )
        {
            chop;
            /\s([0-9]*)$/;
            $unsat_interfacial_polar = $1;
        }
        if ( /^# unsatisfied interfacial charged atom count/ )
        {
            chop;
            /\s([0-9]*)$/;
            $unsat_interfacial_charged = $1;
        }
        if ( /^# buried carbons/ )
        {
            chop;
            /\s([0-9.]*)$/;
            $buried = $1;
        }
        if ( /^# remaining vdW collisions/ )
        {
            chop;
            /\s([0-9]*)$/;
            $final_bumps = $1;
        }
        if ( /^# total vdW overlap/ )
        {
            chop;
            /\s([0-9.]*)$/;
            $final_overlap = $1;
        }
        if ( /^# anchor fragment translations/ )
        {
            chop;
            /\s([0-9]*)$/;
            $anchor_translate = $1;
        }
        if ( /^# side-chain mean-field iterations/ )
        {
            chop;
            /\s([0-9]*)$/;
            $iter = $1;
        }
        if ( /^# ligand side-chain rotations/ )
        {
            chop;
            /\s([0-9]*)$/;
            $lig_side_chain = $1;
        }
        if ( /^# protein side-chain rotations/ )
        {
            chop;
            /\s([0-9]*)$/;
            $prot_side_chain = $1;
        }
        if ( /^# number of water translations/ )
        {
            chop;
            /\s([0-9]*)$/;
            $wat_trans = $1;
        }
        if ( /^# number of conserved waters/ )
        {
            chop;
            /\s([0-9]*)$/;
            $wat_cons = $1;
        }
        if ( /^# number of displaced waters/ )
        {
            chop;
            /\s([0-9]*)$/;
            $wat_disp = $1;
        }
        if ( /^# number of polar displaced waters/ )
        {
            chop;
            /\s([0-9]*)$/;
            $wat_pol_disp = $1;
        }
#       if ( /^# conserved waters/ && $score == $best_scores{$base[$i]} )
#       {
#           @line = split;
#           $best_conserved_waters{$base[$i]} = $#line - 3;
#       }
#       if ( /^# water translations/ && $score == $best_scores{$base[$i]} )
#       {
#           chop;
#           /\s([0-9]*)$/;
#           $best_water_translations{$base[$i]} = $1;
#       }
#       if ( /^# displaced waters/ && $score == $best_scores{$base[$i]} )
#       {
#           @line = split;
#           $best_displaced_waters{$base[$i]} = $#line - 3;
#       }
#       if ( /^# polar displaced waters/ && $score == $best_scores{$base[$i]} )
#       {
#           @line = split;
#           $best_polar_displaced_waters{$base[$i]} = $#line - 4;
#           last;
#       }
        if ( /^# affiscore/ )
        {
            chop;
            /\s([\-0-9.]*)$/;
            $orientscore = $1;
        }
        if ( /^# affiscore/ )
            {
                chop;
                /\s([\-0-9.]*)$/;
                $score = $1;
                $binding_modes{$base[$i]}++;
                if ( $score < $best_scores{$base[$i]} )
                {
                    $best_scores{$base[$i]} = $score;
                    $best_affiscore{$base[$i]} = $affiscore;
                    $best_buried_hphob{$base[$i]} = $buried_hphob;
                    $best_hphob_comp{$base[$i]} = $hphob_comp;
                    $best_polar_component{$base[$i]} = $polar_component;
                    $best_affinity_const{$base[$i]} = $affinity_const;
                    $best_pl_hphob{$base[$i]} = $pl_hphob;
                    $best_pl_hbond{$base[$i]} = $pl_hbond;
                    $best_pl_saltbridge{$base[$i]} = $pl_saltbridge;
                    $best_ml_interactions{$base[$i]} = $ml_interactions;
                    $best_unsat_interfacial_polar{$base[$i]} = $unsat_interfacial_polar;
                    $best_unsat_interfacial_charged{$base[$i]} = $unsat_interfacial_charged;
                    $best_buried{$base[$i]} = $buried;
                    $best_final_bumps{$base[$i]} = $final_bumps;
                    $best_final_overlap{$base[$i]} = $final_overlap;
                    $best_orientscore{$base[$i]} = $orientscore;
                    $best_anchor_translate{$base[$i]} = $anchor_translate;
                    $best_iter{$base[$i]} = $iter;
                    $best_lig_side_chain{$base[$i]} = $lig_side_chain;
                    $best_prot_side_chain{$base[$i]} = $prot_side_chain;
                    $best_wat_trans{$base[$i]} = $wat_trans;
                    $best_wat_cons{$base[$i]} = $wat_cons;
                    $best_wat_disp{$base[$i]} = $wat_disp;
                    $best_wat_pol_disp{$base[$i]} = $wat_pol_disp;
                    $best_file{$base[$i]} = $files[$i];
                    $best_file{$base[$i]} =~ s/.mol2$//;
                }
            }
    }
    close IN;
}
	     

$rank = 1;
print "set mouse insight\n";

foreach $key ( sort sort_by_best_affiscore keys %best_affiscore )
{
    $code = $key;
    if ( $ARGV[2] eq "csd" )	
    {
	$code =~ s/_.*$//;
	$code =~ /^(.)/;
	$filename = $CSD_SET_DATA_DIR."/".$1.".data";
	open IN, "$filename";
	$name = "unknown";
	while ( <IN> )
	{
	    if ( /^CSD-code: $code/ )
	    {	    
		$_ = <IN>;
		chop;
		s/name: //;
		$name = $_;
		last;
		close IN;
	    }
	}
    }
    print "zap\n";		 
    print "load mol2 ".$ligand_dir."/".$best_file{$key}."\n";
    print "wireframe 0.1\n";
    print "echo \"\"\n";
    print "echo \"File:             $best_file{$key}\"\n";
    print "echo \"Rank:             $rank\"\n";
    if ( $ARGV[2] eq "csd" )	
    {
	print "echo \"CSD-code:         $code\"\n";
	$len = length ( $name );
	print "echo \"Compound:         "; 
	while ( $len > 40 )
	{
	    $name =~ s/^\s*(.{40}\S*)//;
	    print "$1\"\necho \"                  "; 
	    $name =~ s/^\s+//;
	    $len = length ( $name );
	}
	print "$name\"\n"; 
    }
    print "echo \"Affinity Score:                      $best_affiscore{$key}\"\n";
    print "echo \"Buried Protein Hydrophobic Term:     $best_buried_hphob{$key}\"\n";
    print "echo \"Hydrophobic Complementarity Term:    $best_hphob_comp{$key}\"\n";
    print "echo \"Polar Component Term:                $best_polar_component{$key}\"\n";
    print "echo \"Protein-ligand Hydrophobic Contact:  $best_pl_hbond{$key}\"\n";
    print "echo \"Protein-ligand H-bonds:              $best_pl_hbond{$key}\"\n";
    print "echo \"Protein-ligand Salt-bridges:         $best_pl_saltbridge{$key}\"\n";
    print "echo \"Metal-ligand Interactions:           $best_ml_interactions{$key}\"\n";
    print "echo \"Unsatisfied Interface Polar Atoms:   $best_unsat_interfacial_polar{$key}\"\n";
    print "echo \"Unsatisfied Interface Charged Atoms: $best_unsat_interfacial_charged{$key}\"\n";
    print "echo \"Buried Carbons:                      $best_buried{$key}\"\n";
    print "echo \"Remaining van der Waals Bumps:       $best_final_bumps{$key}\"\n";
    print "echo \"van der Waals Overlap:               $best_final_overlap{$key}\"\n";
    print "echo \"Ligand Side-chain Rotations:         $best_lig_side_chain{$key}\"\n";
    print "echo \"Protein Side-chain Rotations:        $best_prot_side_chain{$key}\"\n";
    print "pause\n";
    $rank++;
}
