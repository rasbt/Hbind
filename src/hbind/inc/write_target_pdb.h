#ifndef _WRITE_TARGET_PDB_H
#define _WRITE_TARGET_PDB_H
#include "types.h"

void write_target_pdb(residue_pt residues, int num_res, atom_pt atoms, 
                      int num_atoms, char *filename, int *rotations, 
                      float *positions);

#endif

