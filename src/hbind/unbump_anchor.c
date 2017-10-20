/*
 *  This file contains routines to correct the position of ligands
 *  which failed the bump-check, but only bump no more than at
 *  MAX_BUMP_POINTS positions. The main routine in this file is
 *  'unbump_anchor()', which tries to move the ligand minimally away from
 *  the bumping protein atoms to get rid of all bumps. After this
 *  the bump-check has to be repeated for the new ligand position,
 *  since new bumps might have been caused by the transformation.
 */

#include <math.h>
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "unbump_translate.h"
#include "basics.h"
/*
 *  This function translates the ligand according to vector 'vector' 
 */
void  translate_ligand ( global_data_pt global,
			 vector_t       vector )
{
  atom_pt   atom;
  int       i, j;

  atom = global->ligand->atoms;
  for ( i = 0; i < global->ligand->number_of_atoms; i++ )
    {
      for ( j = 0; j < 3; j++ )
	atom->pos[j] += vector[j];
      atom++;
    }
}


/*
 *  The name of this function might be confusing: it computes the minimal
 *  translation along the global vector that is needed to get rid of
 *  all bumps, and so the maximal factor of all bumps has to be found.
 */
int compute_minimal_translation ( vector_pt target_vectors,
				  bump_pt   bump,
				  int       number_of_bumps,
				  vector_t  global_vector,
				  double    overlap,
				  double    *maximal_translation )
{
  int    i;
  double translation;
  
  *maximal_translation = 0.0;
  for ( i = 0; i < number_of_bumps; i++ )
    {
      if ( compute_translation_coefficient ( target_vectors[i],
					     global_vector,
					     bump[i].protein_atom_rad
					     + bump[i].ligand_atom_rad
					     - overlap
					     + 0.00001,
					     &translation ) == FAILURE )
	return FAILURE;
      if ( compare_double(translation, *maximal_translation) == 1 )
	  *maximal_translation = translation;
    }
  return SUCCESS;
}

int  unbump_anchor ( global_data_pt global,
		     bump_pt        bump,
		     int            number_of_bumps )
{
  vector_t  target_vectors[MAX_ANCHOR_BUMPS]; /* translation vector to ideal
					       position for each bump_point */
  vector_t  global_vector;        /* the global averaged translation vector */
  int       i, j;
  double    factor;

  for ( i = 0; i < number_of_bumps; i++ )
    for ( j = 0; j < 3; j++ )
      target_vectors[i][j] =
	bump[i].ligand_atom_pos[j] - bump[i].protein_atom_pos[j];

  compute_global_vector ( target_vectors,
			  number_of_bumps,
			  global_vector );

  if ( compute_minimal_translation ( target_vectors,
				     bump,
				     number_of_bumps,
				     global_vector,
				     (double) global->anchor_overlap,
				     &factor ) == FAILURE )
    return FAILURE;

  for ( i = 0; i < 3; i++ )
    global_vector[i] *= factor;

  if ( fabs ( global_vector[0] ) > global->anchor_translation
       || fabs ( global_vector[1] ) > global->anchor_translation
       || fabs ( global_vector[2] ) > global->anchor_translation )
    /* at least one component of the translation vector is larger than
       the threshold, i.e., one can no longer assume that the template
       points are matched, so this main-chain collisions are not
       resolvable */
    return FAILURE;
  translate_ligand ( global,
		     global_vector );
  return SUCCESS;
}  


