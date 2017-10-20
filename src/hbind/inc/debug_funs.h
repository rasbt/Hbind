
#ifndef DEBUG_FUNCTIONS_HEADER_FILE_INCLUDED
#define DEBUG_FUNCTIONS_HEADER_FILE_INCLUDED

void write_mol2(char *filename, float *positions, global_data_pt global);

void write_pdb(atom_pt atoms, int num_atoms, char *filename, float *positions);

const char* act2str(int act);

#endif
