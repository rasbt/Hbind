#ifndef _READ_PDB_H
#define _READ_PDB_H

void read_pdb(char *fname, atom_pt atoms, residue_pt residues,
              int what_to_read, int *number_of_atoms, int *number_of_residues);

#endif
