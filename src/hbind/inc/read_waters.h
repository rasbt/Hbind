#ifndef _WATERS_H
#define _WATERS_H
#include "types.h"

int find_hbond_partner(atom_pt target, int num_target_atoms, atom_pt water);

int read_waters(atom_pt waters, int *num_waters, char *filename, atom_pt target,
                int num_target_atoms);

#endif
