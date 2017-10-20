#ifndef _TRANSFORM_MOLECULE_H
#define _TRANSFORM_MOLECULE_H

int transform_molecule(global_data_pt global, int *template_matches,
                       int *matched_type, int *ligand_matches);
void reset_target_and_waters(global_data_pt global);

#endif

