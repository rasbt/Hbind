#ifndef  _CHECK_COMPLEMENTARITY_H
#define  _CHECK_COMPLEMENTARITY_H
#include "types.h"

int is_repulsive_polar( global_data_pt global,
                        atom_t lig, int lig_atom_idx, /* ligand*/
                        atom_t tgt, int tgt_atom_idx, /* target*/
                        float lig_tgt_dist);

int is_repulsive_charged( global_data_pt global,
			  atom_t lig, int lig_atom_idx, /* ligand*/
			  atom_t tgt, int tgt_atom_idx, /* target*/
                        float lig_tgt_dist);
						


int is_hbond(atom_t a, atom_t b, int i, int j, float dist, float *angle, 
             global_data_pt  global);
int is_salt_bridge(int a_act, float a_charge, int b_act, float b_charge, 
                   float dist);
int is_metal_hbond(atom_t a, atom_t metal, float dist );
int is_metal_salt_bridge(atom_t lig, atom_t metal, float dist);
int is_hphobic_complementary( atom_t a , atom_t b, float dist);

static inline int
is_hbond_interaction(int A_act, int B_act)
{
  int act_sum = A_act + B_act;
  if(ACCEPTOR_DONOR == act_sum || DONOR_DONEPTOR == act_sum ||
     ACCEPTOR_DONEPTOR == act_sum || DONEPTOR_DONEPTOR == act_sum) return 1;
  return 0;
}

static inline int
is_polar_interaction(int A_act, int B_act)
{
  int act_sum = A_act + B_act;
  if(ACCEPTOR_DONOR == act_sum || DONOR_DONEPTOR == act_sum ||
     ACCEPTOR_DONEPTOR == act_sum || DONEPTOR_DONEPTOR == act_sum ||
     ACCEPTOR_DONOR_HYDROGEN == act_sum || DONEPTOR_DONOR_HYDROGEN == act_sum ||
     ACCEPTOR_METAL_1 == act_sum || ACCEPTOR_METAL_2 == act_sum ||
     DONEPTOR_METAL_1 == act_sum || DONEPTOR_METAL_2 == act_sum) return 1;
  return 0;
}

#endif
