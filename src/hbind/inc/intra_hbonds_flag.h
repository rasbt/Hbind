#ifndef  _INTRA_HBONDS_REC_H
#define  _INTRA_HBONDS_REC_H

int intra_target_hbonds_flag(global_data_pt global, int index, 
                             short target_flag[MAX_PDB_ATOMS], FILE *fout);

extern void target_to_water_hbonds_flag ( global_data_pt  global,
					  int             index,
					  short           target_flag[MAX_PDB_ATOMS] );

extern void intra_ligand_hbonds_flag ( global_data_pt  global, 
				       int             index,
				       short           ligand_flag[MAX_PDB_ATOMS] );

extern void ligand_to_water_hbonds_flag ( global_data_pt  global,
					  int             index,
					  short           ligand_flag[MAX_PDB_ATOMS] );

int intra_target_polar_flag(residue_pt target_residues, atom_pt target_atoms,  
                            const int num_targ_atoms, short *target_flag, 
                            FILE *fout, float *num_of_intra_target_salt_bridges,
                            global_data_pt global);
#endif

