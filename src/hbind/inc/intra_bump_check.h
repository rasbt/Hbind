#ifndef _INTRA_BUMP_CHECK_H
#define _INTRA_BUMP_CHECK_H

extern int  intra_ligand_bump ( molecule_pt ligand,
				float       max_overlap,
				int         *rotated_atoms,
				int         number_of_rotated_atoms );

extern int  check_for_intra_ligand_bump ( global_data_pt      global,
					  transform_matrix_pt matrix,
					  int                 flex_bond,
					  double              angle );

extern void  identify_second_grade_neighbors ( global_data_pt  global );

extern float  intra_molecular_bump_check ( atom_pt  atoms,
					   int      number_of_atoms );

extern double  bump_check_neighbors ( atom_pt  atoms,
				      int      *atoms_to_check,
				      int      number_of_atoms_to_check,
				      float    allowed_overlap );

#endif
