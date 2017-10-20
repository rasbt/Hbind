#ifndef  _WRITE_WATERS_H
#define  _WRITE_WATERS_H
#include "types.h"

/*! Write the conserved waters 
 *
 * If either water_states OR water_positions is NULL the states and positions
 * are taken from the waters array.  Otherwise the positions and states
 * are taken from the water_positions and water_arrays respectively
 */
void write_waters_pdb(atom_pt waters, int num_waters, char *filename, 
                      int *water_states, float *water_positions);
#endif


