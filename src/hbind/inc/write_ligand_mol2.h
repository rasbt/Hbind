#ifndef _WRITE_LIGAND_MOL2_H
#define _WRITE_LIGAND_MOL2_H
#include "types.h"

/*! Write out the ligand corresponding to a docking
 *
 * If positions is NULL, use the positions in global.  Otherwise, use the 
 * positions vector for the atoms' positions.  
 */
void write_ligand_mol2(char *filename, float *positions, 
                       dock_feats_pt features, global_data_pt global);

#endif

