#include <stdio.h>
#include <float.h>
#include <math.h>
#include <mymalloc.h>
#include <err_handle.h>
#include <dist_fun.h>
#include <distance_matrices.h>

static const float MAX_DIST = 8.0;
static const float MAX_SQRD_DIST = 64.0;
static const float SQRD_NBR_DIST = 64.0;


int non_covalently_bound2(atom_pt atom1, atom_pt atom2);

void
init_target_nbr_arrays(atom_pt target_atoms, const int num_atoms)
{
  int i, j; 
  float sq_dist;

  for(i = 0; i < num_atoms; ++i){
    for(j = 0; j < i; ++j){
      sq_dist = squared_dist(target_atoms[i].pos, target_atoms[j].pos);
      if(sq_dist <= SQRD_NBR_DIST && 
         non_covalently_bound2(&target_atoms[i], &target_atoms[j])){
        target_atoms[i].neighbors[target_atoms[i].num_nbrs] = &target_atoms[j];
        target_atoms[j].neighbors[target_atoms[j].num_nbrs] = &target_atoms[i];
        target_atoms[i].neighbor_dist[target_atoms[i].num_nbrs] = 
          target_atoms[j].neighbor_dist[target_atoms[j].num_nbrs] = 
          sqrt(sq_dist);
        ++(target_atoms[i].num_nbrs);
        ++(target_atoms[j].num_nbrs);
      }
    }
  }
}

int non_covalently_bound2(atom_pt atom1, atom_pt atom2)
{
  /* both atoms in the same residue, if level difference is 0, 1, 2 or 3
   * then they are covalently bound */
  if(atom1->residue_index == atom2->residue_index){
    if(diff_int(atom1->level, atom2->level) > 3 && atom1->res != HIS &&
       atom1->res != PRO && atom1->res != TYR && atom1->res != PHE &&
       atom1->res != TRP) return TRUE;
    else return FALSE;
  }
  /* atoms are main-chain atoms of neighbored residues, if they are both main 
   * chain atoms, it is assumed that they are covalently bound. this is ok 
   * since we do not change the main-chain conformation, thus there is no need 
   * to check for a bump here */
  if(diff_int ( atom1->residue_index, atom2->residue_index) == 1 &&
     atom1->level == ALPHA && atom2->level == ALPHA) return FALSE;
  /* disulfide bond */
  if(atom1->res == CYS && atom1->type == SG &&
     atom2->res == CYS && atom2->type == SG )return FALSE;
  return TRUE;
}

void
initialize_inter_dist_matrix(const float* targ_positions, 
                             const size_t num_targ_atoms,
                             const float* lig_positions, 
                             const size_t num_lig_atoms,
                             float *M, const float max_dist)
{
  size_t i, j;
  float *m = M;
  float sq_dist;
  const float* l_pos;
  const float* t_pos;
  const float max_sqrd_dist = max_dist * max_dist;
  const float* targ_pos_end = targ_positions + 3*num_targ_atoms;
  const float* lig_pos_end = lig_positions + 3*num_lig_atoms;

  for(l_pos = lig_positions; l_pos < lig_pos_end; l_pos += 3)
    for(t_pos = targ_positions; t_pos < targ_pos_end; t_pos += 3, ++m){
      sq_dist = squared_dist(l_pos, t_pos);
      *m = (sq_dist <= max_sqrd_dist ? sq_dist : FLT_MAX);
    }
}
