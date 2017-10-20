#ifndef _FIND_FLEXIBLE_BONDS_H
#define _FIND_FLEXIBLE_BONDS_H

extern int  count_flexible_ligand_bonds ( molecule_pt     molecule,
                                          short           ligand_flag[MAX_LIGAND_ATOMS] );

extern int  count_flexible_target_bonds ( global_data_pt  global,
                                          short           target_flag[MAX_PDB_ATOMS] );

#endif
