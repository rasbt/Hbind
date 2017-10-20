#ifndef _TRANS_ROTATE_H

#define _TRANS_ROTATE_H

extern void  compute_yzx_matrix ( double             angle[3],
				  double             trans[3],
				  transform_matrix_t matrix );

extern void  compute_transformation_matrix ( float              n[3],
					     float              ca[3],
					     float              cb[3],
					     transform_matrix_t matrix );

extern void  transform_point ( float              pos[3],
			       transform_matrix_t matrix );

extern void  transform_point_back ( float              pos[3],
				    transform_matrix_t matrix );

extern void  rotate_around_y_axis ( atom_pt atoms,
				    int     number,
				    int     level,
				    double  angle );

extern void  rotate_single_atom_around_y_axis ( float   *position,
						double  angle );

#endif
