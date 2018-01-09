unlink "./src/hbind/main.o";
unlink "./src/hbind/screen_single_compounds.o";
unlink "./src/hbind/read_pts_file.o";
unlink "./src/hbind/mean_field_minimization.o";
unlink "./src/hbind/read_parameter_file.o";
unlink "./src/hbind/match_triangles.o";
unlink "./src/hbind/distance_array.o";
unlink "./src/hbind/check_connectivity.o";
unlink "./src/hbind/find_flexible_bonds.o";
unlink "./src/hbind/read_rule.o";
unlink "./src/hbind/read_waters.o";
unlink "./src/hbind/print_interaction.o";
unlink "./src/hbind/check_complementarity.o";
unlink "./src/hbind/quicksort.o";
unlink "./src/hbind/docking_features.o";
unlink "./src/hbind/transform_molecule.o";
unlink "./src/hbind/intra_hbonds.o";
unlink "./src/hbind/find_hyd_atoms.o";
unlink "./src/hbind/check_rule.o";
unlink "./src/hbind/rotate_unbump_bonds.o";
unlink "./src/hbind/unbump_anchor.o";
unlink "./src/hbind/read_flex_defn.o";
unlink "./src/hbind/least_square_fit.o";
unlink "./src/hbind/read_template.o";
unlink "./src/hbind/intra_bump_check.o";
unlink "./src/hbind/insertion_sort.o";
unlink "./src/hbind/unbump_translate.o";
unlink "./src/hbind/write_log_file.o";
unlink "./src/hbind/read_pdb.o";
unlink "./src/hbind/compute_target_angles.o";
unlink "./src/hbind/main_score.o";
unlink "./src/hbind/hbind_score.o";
unlink "./src/hbind/analyze_ligand.o";
unlink "./src/hbind/write_waters_pdb.o";
unlink "./src/hbind/number_ligand_atoms.o";
unlink "./src/hbind/initialize_unbump_matrices.o";
unlink "./src/hbind/intra_hbonds_flag.o";
unlink "./src/hbind/trace.o";
unlink "./src/hbind/sum_charges.o";
unlink "./src/hbind/read_hyd_defn.o";
unlink "./src/hbind/distance_matrices.o";
unlink "./src/hbind/bump_check.o";
unlink "./src/hbind/rotate.o";
unlink "./src/hbind/find_all_bumps.o";
unlink "./src/hbind/compute_all_rotation_angles.o";
unlink "./src/hbind/bitstrings.o";
unlink "./src/hbind/count_flexible_bonds.o";
unlink "./src/hbind/initialize.o";
unlink "./src/hbind/write_target_pdb.o";
unlink "./src/hbind/assign_hydrogens.o";
unlink "./src/hbind/eigen.o";		
unlink "./src/hbind/unbump_side_chains.o";
unlink "./src/hbind/compute_ligand_angles.o";
unlink "./src/hbind/calc_score_from_terms.o";
unlink "./src/hbind/assign_type.o";
unlink "./src/hbind/hashing.o";
unlink "./src/hbind/assign_fragments.o";
unlink "./src/hbind/score_complex.o";
unlink "./src/hbind/find_carbon_ring_centers.o";
unlink "./src/hbind/debug_funs.o";
unlink "./src/hbind/adj_list.o";
unlink "./src/hbind/hbond_check.o";
unlink "./src/hbind/dist_fun.o";
unlink "./src/hbind/compute_unbump_dependencies.o";
unlink "./src/hbind/read_mol2.o";
unlink "./src/hbind/find_cycles.o";
unlink "./src/hbind/trans_rotate.o";
unlink "./src/hbind/write_ligand_mol2.o";
unlink "./src/utils/hbind_itimer.o";
unlink "./src/utils/err_handle.o";
unlink "./src/utils/basics.o";
unlink "./bin/hbind";


$hbind_dir = '.';
$perl = $^X;

mkdir "$hbind_dir";
mkdir "$hbind_dir/bin";

# unpack the Perl scripts and change the path to the perl binary
open IN, "$hbind_dir/src/scripts/perl_scripts";
print "\n*** extracting Perl scripts ***\n";
while ( <IN> )
{
    if ( /^____FILE\_START/ )
    {
	chop;
	# grep the name of packed Perl script
	@line = split;
	$file = "$hbind_dir/bin/$line[1]";
	open OUT, ">$file";
	# read the line that includes the path to the Perl binary
	$_ = <IN>;
	print OUT "\#\!$perl\n";
	next;
    }
    if ( /^____FILE\_END/ )
    {
	# end of this script, close filehandle and change permissions
	close OUT;
	chmod 0755, "$file";
	next;
    }
    print OUT $_;
}
close IN;

# unpack the shell scripts
open IN, "$hbind_dir/src/scripts/shell_scripts";
print "\n*** extracting Shell scripts ***\n";
while ( <IN> )
{
    if ( /^____FILE\_START/ )
    {
	chop;
	# grep the name of packed Perl script
	@line = split;
	$file = "$hbind_dir/bin/$line[1]";
	open OUT, ">$file";
	next;
    }
    if ( /^____FILE\_END/ )
    {
	# end of this script, close filehandle and change permissions
	close OUT;
	chmod 0755, "$file";
	next;
    }
    print OUT $_;
}
close IN;

# compile the libHBIND_utils first since the other programs will need to
# link to it
chdir "src/utils";
print "\n*** compiling the static hbind utilties library ***\n\n";
system "make";

# compile HBIND
chdir "../hbind";
print "\n*** compiling HBIND ***\n\n";
system "make";
# compile compute_interaction_centers, check_connectivity, and
# generate_rasmol_script
chdir "../interactions";
print "\n*** compiling interactions auxiliaries ***\n\n";
system "make";
# compile average_template and unbiased_template
chdir "../template";
print "\n*** compiling template auxiliaries ***\n\n";
system "make";
print "\n*** HBIND installation completed ***\n";
print "\n\n*** Please see license.txt for licensing and use information ***\n\n";
