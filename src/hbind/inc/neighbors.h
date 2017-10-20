#ifndef _NEIGHBORS_H

#define _NEIGHBORS_H
/*
void identify_atom_neighbors(atom_pt atoms, int number_of_atoms);
*/

extern void  identify_water_neighbors ( atom_pt  waters,
					atom_pt  target_atoms,
					int      number_of_waters,
					int      number_of_target_atoms,
					global_data_pt global);
#endif

