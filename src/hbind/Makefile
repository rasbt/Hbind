BINPATH = ../../bin

CC = gcc -v
EXENAME = hbind
COPTS = -O2 # Use for distribution

INCLUDES = -Iinc -I../utils/inc
VERSION = -DVERSION=\"1.0.0\"
CFLAGS = $(COPTS) $(INCLUDES) $(VERSION)
LIBS = ../utils/libhbindutils.a -lm

OBJS =	debug_funs.o \
	initialize.o  \
        distance_array.o \
	distance_matrices.o \
	read_mol2.o \
	assign_type.o \
	read_pdb.o \
	dist_fun.o \
	bitstrings.o  \
	read_flex_defn.o \
	read_rule.o \
	adj_list.o \
	assign_hydrogens.o \
	find_cycles.o \
	find_flexible_bonds.o \
	check_rule.o \
	score_complex.o \
	trans_rotate.o \
	read_waters.o \
	hbond_check.o \
	find_hyd_atoms.o \
	read_hyd_defn.o \
	intra_hbonds.o \
	intra_hbonds_flag.o \
	sum_charges.o \
	calc_score_from_terms.o \
	check_complementarity.o \
	count_flexible_bonds.o \
	trace.o \
	hashing.o  \
	read_pts_file.o \
	read_template.o \
	analyze_ligand.o \
	transform_molecule.o \
	assign_fragments.o \
	bump_check.o \
	check_connectivity.o \
	compute_all_rotation_angles.o \
	unbump_side_chains.o \
	compute_target_angles.o \
	compute_ligand_angles.o \
	unbump_anchor.o \
	unbump_translate.o \
	find_all_bumps.o \
	compute_unbump_dependencies.o \
	eigen.o \
	rotate.o \
	least_square_fit.o \
	find_carbon_ring_centers.o \
	initialize_unbump_matrices.o \
	insertion_sort.o \
	rotate_unbump_bonds.o \
	intra_bump_check.o \
	mean_field_minimization.o \
	number_ligand_atoms.o \
	print_interaction.o \
	quicksort.o \
	read_parameter_file.o \
	write_target_pdb.o \
	write_waters_pdb.o \
	write_log_file.o \
	write_ligand_mol2.o  \
	docking_features.o \
    match_triangles.o



$(BINPATH)/$(EXENAME): $(OBJS) main_hbind.o
	$(CC) $(COPTS) -o $(BINPATH)/$(EXENAME) $(OBJS) main_hbind.o \
	$(LIBS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o *.gcda *.gcno $(BINPATH)/$(EXENAME) $(BINPATH)/$(EXENAME)_hbind

touch:
	touch *.c
