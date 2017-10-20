#ifndef _UNBUMP_SIDE_CHAINS_H
#define _UNBUMP_SIDE_CHAINS_H

extern void  initialize_distance_table ( global_data_pt  global );

extern int  compute_matrix ( atom_pt            side_chain_atoms,
			     atom_pt            ligand_atom,
			     int                res_type,
			     int                bond,
			     transform_matrix_t matrix );

extern double  compute_unbump_angle ( double  x2,
				      double  d,
				      double  y1,
				      double  y2,
				      double  lxz,
				      double  l );

extern int  check_resolvability ( float              fixed_atom_pos[3],
				  float              rotated_atom_pos[3],
				  float              final_distance,
				  transform_matrix_t matrix );

extern int  unbump_side_chains ( global_data_pt  global );

#endif

