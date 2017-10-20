#ifndef  _MATCH_TRIANGLES_H
#define  _MATCH_TRIANGLES_H
#include "types.h"

void match_triangles(global_data_pt global);

int compute_triangles(global_data_pt global, int *array, int index, int length,
                      int offset);

void report_docking(dock_feats_pt features, moved_positions_pt moved_positions,
                    char* lig_fname, global_data_pt global);

int score(dock_feats_pt features, global_data_pt global, FILE *good_acts_file);

     
#endif
