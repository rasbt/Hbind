#!perl
#
# results_table.pl
#
# modified version of the below cited perl scripts
#                            Maria Zavodszky Thu Jan 27 2005
#
# create_global_genericDB_results_table.pl
#                            Paul Sanschagrin Wed Mar 15 2000
# modified version of create_global_csd_results_table.pl (written by)
#                            Volker Schnecke  Wed Sep 15 15:46:25 EDT 1999
#
# Usage: create_global_genericDB_results_table.pl target template database number
#
# This script creates a table listing all important data for the top potential
# ligands that HBIND has identified.
                                                                                
                                                                             
                                                                                if ( $#ARGV < 3 )
{
#    die "Usage: $0 <target> <template> <database> <# displayed results> [-confs/-confb/-bumps]\n\nOption:\n[-confs] Groups all conformations of a ligand, separating bumped & unbumped\n[-confb] Groups all conformations without separating bumped & unbumped\n[-bumps] Doesn't group conformations together, but groups bumped & unbumped\n";
    die "Usage: $0 <target> <template> <database> <number> [-conf]\n";
}
                                                                                
sub sort_by_best_orientscore {
    $best_scores{$a} <=> $best_scores{$b};
}

sub sort_by_best_affiscore {
    $best_affiscore{$a} <=> $best_affiscore{$b};
}
                                                                                
$ligand_dir = $ENV{HBIND_DATA_DIR}."/".$ARGV[0]."/".$ARGV[1]."/".$ARGV[2]."_ligands";
$param_file = $ENV{HBIND_DATA_DIR}."/".$ARGV[0]."/".$ARGV[1]."/in/hbind.parameters";
$conformers = $ARGV[4];
                                                                                
opendir DIR, $ligand_dir;
@files = grep ( /.mol2$/, grep ( !/^.\.?$/, readdir ( DIR ) ) );
close DIR;

for ( $i = 0; $i <= $#files; $i++ )
{
    $base[$i] = $files[$i];
    $base_a[$i] = $files[$i];
    $base_b[$i] = $files[$i];
    
    if ($conformers =~ /^-confs/){
	$base_a[$i] =~ s/_\d+_\D?_\d+.mol2$//; 
	$base_b[$i] =~ s/_\d+.mol2$//;
	$length_bb = length($base_b[$i]);
	$base_c[$i] = substr($base_b[$i], $length_bb-1, 1);
	$base[$i] = "$base_a[$i]_$base_c[$i]";
#    print "ba = $base_a[$i]\n";
#    print "bb = $base_b[$i]\n";
#    print "bc = $base_c[$i]\n";
#    print "basei = $base[$i]\n";
    }
#    elsif ($conformers =~ /^-confb/){
    elsif ($conformers =~ /^-conf/){
	$base[$i] =~ s/_\d+_\d+.mol2$//;  # for: ligandname_<conformation>b_<orientation>.mol2
    }
    elsif ($conformers =~ /^-bumps/){
	$base[$i] =~ s/_\D?_\d+.mol2$//; 
    }
    else {
	$base[$i] =~ s/_\d+.mol2$//;      # for: ligandname_<orientation>.mol2
    }
#    print "$base[$i]\n";
    $best_score{$base[$i]} = 0.0;
    $binding_modes{$base[$i]} = 0;
}

for ( $i = 0; $i <= $#files; $i++ )
{
    open IN, "$ligand_dir/$files[$i]";
    while ( <IN> )
    {
        if ( /^# Affiscore \(kcal/ )
        {
            chop;
            /\s([\-0-9.]*)$/;
            $affiscore = $1;
        }

        if ( /^# hydrophobic complementarity term/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $hphob_comp = $1;
        }
        if ( /^# polar term/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $polar_component = $1;
        }
        if ( /^# unsatisfied polar term/ )
        {
            chop;
            /\s([\-\0-9.]*)$/;
            $unsat_polar_component = $1;
        }

        if ( /^# constant term/ )
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
        if ( /^# Orientscore/ )
        {
            chop;
            /\s([\-0-9.]*)$/;
            $orientscore = $1;
        }
        if ( /^# Orientscore/ )
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
                    $best_unsat_polar_component{$base[$i]} = $unsat_polar_component;
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

                                                                                
$total_ligands = scalar ( keys %best_affiscore );
$max_ligands = $ARGV[3];
$max_ligands = $total_ligands if $total_ligands < $max_ligands;
$title = "TOP ".$max_ligands." LIGANDS (OUT OF ".$total_ligands.") FOR ".$ARGV[0]." (TEMPLATE ".$ARGV[1].")";
             print "$title\n\n";
	     if ($conformers =~ /^-confs/){
#		 print "You grouped ligands by conformation, keeping bumped and unbumbed dockings separate.\n"
		 }
	     elsif ($conformers =~ /^-conf/){
#		 print "You grouped ligands by conformation.\n"
		 }
	     else
	     {
#		 print "You treated each ligand separately.\n"
	     }
 #            print "  hbind.parameters\n(those currently in your in/ directory, not necessarily the ones used for the run)\n\n";
                                                                                
#open IN, "$param_file";
#             while  ( <IN> ){
#                 print "$_";
#             }
#close IN;
	     print "\nFields:   1  Rank\n";
             print "          2  Ligand name\n";
             print "          3  Affinity score\n";
             print "          4  Hydrophobic complementarity term\n";
             print "          5  Polar component term\n";
             print "          6  Unsatisfied polar component term\n";
             print "          7  Protein-ligand hydrophobic contact count\n";
             print "          8  Protein-ligand H-bond count\n";
             print "          9  Protein-ligand salt-bridge count\n";
             print "         10  Metal-Ligand interaction count\n";
             print "         11  Unsatisfied Interfacial Polar Atom count\n";
             print "         12  Unsatisfied Interfacial Charged Atom count\n";
             print "         13  Buried Carbons\n";
             print "         14  Remaining van der Waals collisions\n";
             print "         15  Remaining van der Waals overlap\n";
             print "         16  Ligand side-chain rotations\n";
             print "         17  Protein side-chain rotations\n";
             print "         18  Top scoring orientation\n";
             print "         19  Binding modes\n\n";
             print "Affinity Constant = -5.218\n
                  Affinity Score = (4) +  (5) + (6) + (Affinity Constant)\n\n";
	     
#23456789012345678901234567890123456789012345678901234567890123456789012345678901234567890
             print "
  1      2                  3       4      5      6      7  8  9 10 11 12   13   14   15  16 17       18             19\n\n";
                                                                                
                                                                                
$rank = 1;
foreach $key ( sort sort_by_best_affiscore keys %best_affiscore )
{
    $code = $key;
    printf "%-3d ", $rank;
    printf "%-20s ", $code;
#   printf "%6.3f ", $best_orientscore{$key};
    printf "%6.3f ", $best_affiscore{$key};
    printf "[%6.3f ", $best_hphob_comp{$key};
    printf "%6.3f ", $best_polar_component{$key};
    printf "%6.3f] ", $best_unsat_polar_component{$key};
    printf "%3d ", $best_pl_hphob{$key};
    printf "%2d ", $best_pl_hbond{$key};
    printf "%2d ", $best_pl_saltbridge{$key};
    printf "%2d ", $best_ml_interactions{$key};
    printf "%2d ", $best_unsat_interfacial_polar{$key};
    printf "%2d ", $best_unsat_interfacial_charged{$key};
    printf "%6.3f ", $best_buried{$key};
    printf "%2d ", $best_final_bumps{$key};
    printf "%5.3f ", $best_final_overlap{$key};
    printf "%2d ", $best_lig_side_chain{$key};
    printf "%2d ", $best_prot_side_chain{$key};
#   printf "%2d ", $best_wat_trans{$key};
#   printf "%2d ", $best_wat_cons{$key};
#   printf "%2d ", $best_wat_disp{$key};
#   printf "%3d ", $best_wat_pol_disp{$key};
    printf "%-20s ", $best_file{$key};
    printf "%2d\n", $binding_modes{$key};
    $rank++;
    last if $rank > $max_ligands;
}
                                                                          
$rank = 1;
foreach $key ( sort sort_by_best_affiscore keys %best_affiscore ) {
    if ($rank >$max_ligands) {
#       print "deleting orientations of $key (rank $rank score $best_scores{$key})\n";
        unlink <$ligand_dir/$key*.mol2>;
        unlink <$target_dir/$key*.pdb>;
        unlink <$water_dir/$key*.pdb>;
    }
    $rank++;
}

