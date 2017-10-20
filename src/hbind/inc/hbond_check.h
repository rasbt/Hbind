#ifndef _HBOND_CHECK_H
#define _HBOND_CHECK_H
#include "types.h"

int  check_hydrogen_angle ( global_data_pt global,
				   int            target_index,
				   int            ligand_index,
				   float          *angle );

float  compute_target_hydrogen_angle ( atom_pt    atoms,
					      residue_pt residues,
					      int        index,
					      float      *acceptor_position,
					      float      *hydrogen_position,
					      global_data_pt  global);

int  compute_ligand_hydrogen_angle ( global_data_pt global,
					    int            ligand_donor_index,
					    float          *acceptor_position,
					    float          *angle );
int
target_preacc_angle(atom_pt acceptor, float *hydrogen_position,
                    atom_pt residue_atoms, int number_of_residue_atoms,
                    float *angle);


#endif
