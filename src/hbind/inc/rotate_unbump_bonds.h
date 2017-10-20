#ifndef _ROTATE_UNBUMP_BONDS_H
#define _ROTATE_UNBUMP_BONDS_H

extern int  rotate_unbump_bonds ( global_data_pt global,
				  int            *number_of_actual_rotations );

extern int  compute_atoms_to_rotate ( global_data_pt  global,
	     		              int             index,
			              int             *atom_indices );

#endif

