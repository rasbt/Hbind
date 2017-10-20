#ifndef  _COMPUTE_ALL_ROTATION_ANGLES_H
#define  _COMPUTE_ALL_ROTATION_ANGLES_H

extern double  compute_rotation_angle ( atom_pt            fixed_atom,
					atom_pt            rotated_atom,
					float              max_overlap,
				        transform_matrix_t matrix );

extern int  compute_all_rotation_angles ( global_data_pt  global );

#endif


