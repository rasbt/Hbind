#ifndef _PRINT_INTERACTION_H_
#define _PRINT_INTERACTION_H_
#include "types.h"
void print_interaction(atom_pt atom1, int atom1_type, int atom1_index, 
                       atom_pt atom2, int atom2_type, int atom2_index,
                       const char *interaction_type, FILE *fout);

#endif
