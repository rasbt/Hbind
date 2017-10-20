#include <float.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "types.h"
#include "basics.h"
#include "dist_fun.h"
#include "bitstrings.h"
#include "trans_rotate.h"
#include "intra_bump_check.h"
#include "insertion_sort.h"
#include <mymalloc.h>

/* #define TRACE */

int compute_atoms_to_rotate(global_data_pt global, int index, int *atom_indices)
{
  int res_index;
  int bond;
  int flex_bond;
  int offset; 
  int i;
  int count = 0;
  atom_pt target = global->target_atoms;
  molecule_pt ligand = global->ligand;
  atom_pt ligand_atoms = ligand->atoms;
  int *unbump_bonds = global->unbump_bonds;

  /* Check if this is a target bond */
  if(unbump_bonds[index] >= 0){
    res_index = unbump_bonds[index] / 10;
    bond = unbump_bonds[index] % 10;
    offset = global->target_residues[res_index].start_atom;
    for(i = 0; i < global->target_residues[res_index].number_of_atoms; i++){
      if(target[offset+i].level > bond) atom_indices[count++] = offset + i;
#ifdef TRACE
      if(target[offset+i].level == bond)
        printf("rotation bond: %s %s %s\n",
               global->target_residues[res_index].name,
               global->target_residues[res_index].num, target[offset+i].name);
#endif
    }
#ifdef TRACE
    printf("atoms to rotate:\n");
    for(i = 0; i < count; i++)
      printf("%4d %s %s %s (%5.3f %5.3f %5.3f)\n",
             atom_indices[i], target[atom_indices[i]].name,
             target[atom_indices[i]].residue, 
             target[atom_indices[i]].residue_num,
             target[atom_indices[i]].pos[X], target[atom_indices[i]].pos[Y],
             target[atom_indices[i]].pos[Z]);
#endif      
  /* otherwise, this is a ligand bond */
  }else{      
    flex_bond = ( (-1) * unbump_bonds[index] ) - 1;
#ifdef TRACE
    printf("unbump_bond = %d, flex_bond = %d, bond = %d, atom1 = %d, "
           "atom2 = %d (%d)\n", index, flex_bond,
           ligand->flexible_bonds[flex_bond],
           ligand->bonds[ligand->flexible_bonds[flex_bond]].atom1 + 1,
           ligand->bonds[ligand->flexible_bonds[flex_bond]].atom2 + 1,
           global->unbump_bonds[index]);
#endif
    if(ligand->bond_directions[flex_bond] == REVERSE){
      for(i = 0; i < ligand->number_of_atoms; i++)
        /* this atom is located 'right' from the rotation bond, i.e., at the 
         * side with the larger atom index, which is the side pointing away 
         * from the base */
        if(bitstring_get_bit(ligand_atoms[i].fragments, flex_bond))
          atom_indices[count++] = i;
    }else{
      for(i = 0; i < ligand->number_of_atoms; i++)
	/* this atom is located 'left' from the rotation bond */
        if(bitstring_get_bit(ligand_atoms[i].fragments, flex_bond) == 0)
          atom_indices[count++] = i;
    }
  }
  return count;
}

void identify_rotation_order(global_data_pt global,
                             sort_array_pt highest_probabilities)
{
  int i, j;

  for(j = 0; j < global->number_of_unbump_bonds; j++){
    highest_probabilities[j].value = 0.0;
    highest_probabilities[j].index1 = UNKNOWN;
    highest_probabilities[j].index2 = j;
#ifdef TRACE
    printf("0. j=%d, highest_probabilities[j].index1=%d, "
           "highest_probabilities[j].index2=%d\n", j, 
           highest_probabilities[j].index1, highest_probabilities[j].index2);
#endif
    for(i = 0; i < global->number_of_bumps; i++ ){
      if(global->unbump_data[i][j] == NULL) continue;
#ifdef TRACE
      printf("  iden:   global->unbump_data[%d][%d]->probability=%.8f\n", i, j,
             global->unbump_data[i][j]->probability);
#endif
      if(compare_double(global->unbump_data[i][j]->probability,
                        highest_probabilities[j].value) == 1){
        highest_probabilities[j].value = global->unbump_data[i][j]->probability;
        highest_probabilities[j].index1 = i;
      }
    }

#ifdef TRACE
    printf("1. j=%d, highest_probabilities[j].index1=%d, "
           "highest_probabilities[j].index2=%d\n", j, 
           highest_probabilities[j].index1, highest_probabilities[j].index2);
#endif
  }
  insertion_sort(highest_probabilities, global->number_of_unbump_bonds);

#ifdef TRACE
#if 0
  /* Not sure what was actually meant to be printed here since i and j are
   * not valid indices */
  printf("2. j=%d, highest_probabilities[j].index1=%d, "
         "highest_probabilities[j].index2=%d\n",i, 
         highest_probabilities[j].index1, highest_probabilities[j].index2);
#endif
#endif
}

int rotate_unbump_bonds(global_data_pt global, int *number_of_actual_rotations)
{
  transform_matrix_pt matrix;
  atom_pt             atom1,
                      atom2;
  sort_array_t        highest_probabilities[MAX_UNBUMP_BONDS];
  int                 atoms_to_rotate[MAX_NUMBER_OF_MOL2_ATOMS];
  double              angle,
                      total_intra_overlap;
  int                 bump,
                      bond,
                      residue,
                      number_of_atoms;
  int                 i, j, k;
  float tmp;
  unbump_data_pt **unbump_data = global->unbump_data;
  atom_pt target = global->target_atoms;
  atom_pt ligand = global->ligand->atoms;
  float orig_pos[3*MAX_NUMBER_OF_MOL2_ATOMS]; // just in case :(
  atom_pt *close_atoms;
  atom_pt *close_atom_p;
  atom_pt close_atom;
  size_t num_atoms;

  identify_rotation_order(global, highest_probabilities);

  /* now actually do the rotations */
  *number_of_actual_rotations = 0;
  for(i = 0; i < global->number_of_unbump_bonds; i++){
    if(highest_probabilities[i].index1 == UNKNOWN) continue;

    /* index1 is bump, index2 is bond */
    bump = highest_probabilities[i].index1;
    bond = highest_probabilities[i].index2;
    angle = unbump_data[bump][bond]->angle;
    matrix = unbump_data[bump][bond]->matrix;

    /* regular ligand atom */
    if(global->ligand_bumps[bump] >= 0)
      atom1 = &ligand[global->ligand_bumps[bump]];
    /* target atom (intramolecular bump) */
    else atom1 = &target[(-1)*(global->ligand_bumps[bump]+1)];
    /* regular target atom */
    if(global->target_bumps[bump] >= 0)
      atom2 = &target[global->target_bumps[bump]];
    /* water molecule */
    else atom2 = &global->waters[(-1)*(global->target_bumps[bump]+1)];

    /* it can happen that there was a dependency of between entries related 
     * to different molecules, e.g., after having done a ligand rotation, the
     * favorable target rotation to resolve the current bump is no longer 
     * necessary, this is not precomputed */
    if( fabs( dist_fun(atom1->pos, atom2->pos) - 
              unbump_data[bump][bond]->old_distance ) >= 0.001) continue;

    /* this is a target bond */
    number_of_atoms = compute_atoms_to_rotate(global, bond, atoms_to_rotate);
    if(global->unbump_bonds[bond] >= 0){
#ifdef TRACE
      printf("This is a target bond, %s %s (bond: %d)\n",
             global->target_residues[global->unbump_bonds[bond]/10].name,
             global->target_residues[global->unbump_bonds[bond]/10].num,
             global->unbump_bonds[bond]%10 );
#endif
      for(j = 0; j < number_of_atoms; j++){
        memcpy(orig_pos + 3*j, target[atoms_to_rotate[j]].pos, 
               3*sizeof(target[atoms_to_rotate[j]].pos[0]));
        transform_point(target[atoms_to_rotate[j]].pos, *matrix);
        rotate_single_atom_around_y_axis(target[atoms_to_rotate[j]].pos, angle);
        transform_point_back(target[atoms_to_rotate[j]].pos, *matrix);
#if 0
          printf("%s %s %s %5.3f %5.3f %5.3f\n",
                 global->target_atoms[atoms_to_rotate[j]].residue,
                 global->target_atoms[atoms_to_rotate[j]].residue_num,
                 global->target_atoms[atoms_to_rotate[j]].name,
                 global->target_atoms[atoms_to_rotate[j]].pos[X],
                 global->target_atoms[atoms_to_rotate[j]].pos[Y],
                 global->target_atoms[atoms_to_rotate[j]].pos[Z] );
#endif
      }
      total_intra_overlap = 
        bump_check_neighbors(target, atoms_to_rotate, number_of_atoms, 
                             global->intra_overlap);

      /* there is a total overlap of more than 1.0 Angstrom, do not accept 
       * this rotation */
      if(total_intra_overlap > 1.0){
#ifdef TRACE
       printf("Intra overlap of %5.3lf, rotation not accepted\n", 
              total_intra_overlap);
#endif
        /* Undo the rotation(s) */
        for(j = 0; j < number_of_atoms; j++)
          memcpy(target[atoms_to_rotate[j]].pos, orig_pos + 3*j, 
                 3*sizeof(target[atoms_to_rotate[j]].pos[0]));
      }else{
        (*number_of_actual_rotations)++;
        global->number_of_target_side_chain_rotations++;
        residue = target[atoms_to_rotate[0]].residue_index;
        global->target_rotations[residue] = YES;

        for(j = 0; j < global->ligand->number_of_atoms; j++)
          for(k = 0; k < number_of_atoms; k++){
            tmp = squared_dist(ligand[j].pos, target[atoms_to_rotate[k]].pos);
            global->target_ligand_sq_dists[j*global->number_of_target_atoms + 
                                           atoms_to_rotate[k]] = 
              (tmp > HYDRO_DIST_2 ? FLT_MAX : tmp);
          }

        if(total_intra_overlap > 0.0)
          global->target_intra_overlap[residue] = YES; 
        else global->target_intra_overlap[residue] = NO;
#ifdef TRACE
        printf ( "Target-bond rotation successful\n" );
#endif
      }
    /* this is a ligand bond */
    }else{
#ifdef TRACE
      printf("This is a ligand bond, rotating bond %d\n", 
             global->unbump_bonds[bond]);
#endif
      for(j = 0; j < number_of_atoms; j++ ){
        memcpy(orig_pos + 3*j, ligand[atoms_to_rotate[j]].pos, 
               3*sizeof(ligand[atoms_to_rotate[j]].pos[0]));
        transform_point ( ligand[atoms_to_rotate[j]].pos, *matrix );
        rotate_single_atom_around_y_axis(ligand[atoms_to_rotate[j]].pos, angle);
        transform_point_back(ligand[atoms_to_rotate[j]].pos, *matrix);
      }

      /* the rotation caused an intramolecular bump in the ligand, so 
       * rotate the side chain back and hope that the corresponding bump 
       * can be resolved in the next iteration */
      if(intra_ligand_bump(global->ligand, global->intra_overlap, 
                           atoms_to_rotate, number_of_atoms) == BUMP){
        for ( j = 0; j < number_of_atoms; j++ ) 
          memcpy(ligand[atoms_to_rotate[j]].pos, orig_pos + 3*j, 
                 3*sizeof(ligand[atoms_to_rotate[j]].pos[0]));
      /* otherwise update the distances for the rotated atoms in the lookup 
       * table */
      }else{ 
        for(j = 0; j < number_of_atoms; j++){
          for(k = 0; k < global->number_of_target_atoms; k++){
            tmp = squared_dist(ligand[atoms_to_rotate[j]].pos, target[k].pos);
            global->target_ligand_sq_dists[atoms_to_rotate[j] *
                                           global->number_of_target_atoms + k] =
              (tmp > HYDRO_DIST_2 ?  FLT_MAX : tmp);
          }
        }
			
        global->number_of_ligand_side_chain_rotations++;
        (*number_of_actual_rotations)++;
#ifdef TRACE
        printf ( "Ligand-bond rotation successfull\n" );
#endif
      } 
    } 
  } 
  
  return SUCCESS;
}
