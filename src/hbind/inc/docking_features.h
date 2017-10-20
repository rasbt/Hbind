#ifndef FEATURES_HEADER_FILE_INCLUDED
#define FEATURES_HEADER_FILE_INCLUDED

#include "types.h"

void save_features(dock_feats_pt features, score_data_pt score,
                   int number_of_bumps, float total_overlap);

void write_features_line(dock_feats_pt features, char *lig_fname, FILE *fout,
                         int docking);

void write_features_header(dock_feats_pt features, FILE *fp, char *lig_fname,
                           interaction_pt template_interactions, 
                           atom_pt ligand_atoms, atom_pt target_atoms);

void init_features(dock_feats_pt features);

#endif
