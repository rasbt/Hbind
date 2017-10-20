#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "types.h"
#include "basics.h"

double  update_mean_force ( global_data_pt  global )
{
  unbump_data_pt  **unbump_data;
  int             **dependents;
  int             *number_of_dependents;
  double          total_mean_force,
                  minimal_mean_force;
  int             dependent,
                  count;
  int             i, j, k, l;

  unbump_data = global->unbump_data;
  dependents = global->unbump_dependents;
  number_of_dependents = global->number_of_unbump_dependents;
  minimal_mean_force = 0.0;
  count = 0;
  for ( i = 0; i < global->number_of_bumps; i++ )
    for ( j = 0; j < global->number_of_unbump_bonds; j++ )
      if ( unbump_data[i][j] != NULL )
	{
	  unbump_data[i][j]->mean_force = 
	    (double) unbump_data[i][j]->force;

	  for ( k = 0; k < global->number_of_bumps; k++ )
	    /* check, if other bumps can be resolved by bond j,
	       this is benefical, i.e., a negative value has to be
	       added to the mean force of (i,j) */
	    if ( k != i )
	      if ( unbump_data[k][j] != NULL )
		{
		  if ( ( unbump_data[k][j]->angle > 0 
			 && unbump_data[i][j]->angle > 0 ) ||
		       ( unbump_data[k][j]->angle < 0 
			 && unbump_data[i][j]->angle < 0 ) )
		    /* both rotations are in the same direction, so by rotating
		       this bond we can resolve two collisions at once */
		      unbump_data[i][j]->mean_force -=
			  unbump_data[k][j]->probability 
			  * (double) unbump_data[k][j]->force;
		  else
		    /* the rotations are in opposite directions, i.e., by
		       applying one rotation, the other bump might get worse,
		       so add a penalty */
		    unbump_data[i][j]->mean_force +=
		      unbump_data[k][j]->probability 
		      * (double) unbump_data[k][j]->force;
		}
	  for ( k = 0; k < number_of_dependents[j]; k++ )
	    {
	      dependent = dependents[j][k];
	      /*printf ("%3d  %3d: %d\n",j,k,dependent);*/
	      for ( l = 0; l < global->number_of_bumps; l++ )
		if ( l != i
		     && unbump_data[l][j] == NULL
		     && unbump_data[l][dependent] != NULL )
		  /* we cannot use bond j to resolve bump l, but bond 
		     'dependent' is dependent of bond j, so by messing it 
		     up when rotating bond 'j', we loose the opportunity 
		     to resolve bump 'l' by rotating bond 'dependent', this 
		     is bad, so we add some positive force to the mean force */
		  unbump_data[i][j]->mean_force += 
		    unbump_data[l][dependent]->probability 
		    * (double) unbump_data[l][dependent]->force; 
	    }
	  if ( compare_double( unbump_data[i][j]->mean_force, minimal_mean_force ) == -1 )
	    minimal_mean_force = unbump_data[i][j]->mean_force;
	}

#ifdef TRACE
  printf ( "Mean-forces:\n" );
  for ( j = 0; j < global->number_of_unbump_bonds; j++ )
    {
      for ( i = 0; i < global->number_of_bumps; i++ )
	if ( unbump_data[i][j] != NULL )
	  printf ( "(%d,%d) force: %6.3lf probability: %6.3lf mean-force: %6.3lf\n", 
		   i, 
		   j,
		   unbump_data[i][j]->force,
		   unbump_data[i][j]->probability,
		   unbump_data[i][j]->mean_force );
	  /*
	  printf ( " %6.3lf", unbump_data[i][j]->mean_force );
	else
	  printf ( "  *.***" );
      printf ( "\n" );
	  */
    }	
#endif

  /* if the minimal mean force is negative, we have to shift all values before
     we can compute the average */
  minimal_mean_force *= -1.0;
  total_mean_force = 0.0;
  for ( j = 0; j < global->number_of_unbump_bonds; j++ )
    for ( i = 0; i < global->number_of_bumps; i++ )
      if ( unbump_data[i][j] != NULL )
	{
	  unbump_data[i][j]->mean_force += minimal_mean_force;
	  total_mean_force += unbump_data[i][j]->mean_force;
	  count++;
	}
  
  if(count > 0) return total_mean_force / count;
  else return 1.0;
}

void  update_probabilities ( global_data_pt  global,
			     double          average_mean_force )
{
  unbump_data_pt  **unbump_data;
  double          sum;
  int             i, j, k;

  unbump_data = global->unbump_data;
  for ( i = 0; i < global->number_of_bumps; i++ )
    for ( j = 0; j < global->number_of_unbump_bonds; j++ )
      if ( unbump_data[i][j] != NULL )
	{
	  sum = 0.0;
	  for ( k = 0; k < global->number_of_unbump_bonds; k++ )
	    if ( unbump_data[i][k] != NULL ) 
		{
	      sum += exp ( (-1) 
			   * unbump_data[i][k]->mean_force 
			   / average_mean_force );

#ifdef TRACE
	      printf("update: i=%d j=%d,k=%d, mean_force=%.16f,  average_mean_force=%.16f\n",i,j,k,unbump_data[i][k]->mean_force, average_mean_force);
#endif

		}
	  unbump_data[i][j]->probability =
	    exp ( (-1) 
		  * unbump_data[i][j]->mean_force 
		  / average_mean_force ) 
	    / sum;
	}

#ifdef TRACE
  for ( i = 0; i < 2; i++ )
    for ( j = 0; j < 2; j++ )
      if ( unbump_data[i][j] != NULL )
	{
	    printf("2. i=%d,j=%d  ",i,j);
	    printf("global->unbump_data[i][j]->probability=%.16f \n",global->unbump_data[i][j]->probability);
	}
#endif

}

void  mean_field_minimization ( global_data_pt  global )
{
  double          average_mean_force;
  int             i, j;
  int             k;

  unbump_data_pt  **unbump_data;  
  unbump_data = global->unbump_data;

#ifdef TRACE
  printf ( "Before Iteration 0:\n" );
  for(j = 0; j < global->number_of_unbump_bonds; j++ ){
    for(i = 0; i < global->number_of_bumps; i++ )
      if(unbump_data[i][j] != NULL )
        printf ( " %6.3lf", unbump_data[i][j]->probability );
      else printf ( "  *.***" );
    printf ( "\n" );
  }	
#endif
  
  for(k = 0; k < 10; k++){
    average_mean_force = update_mean_force(global);

#ifdef TRACE
    printf("average_mean_force = %5.3f\n", average_mean_force );
#endif

    update_probabilities(global, average_mean_force);

#ifdef TRACE
    printf ( "Iteration %d:\n", k );
    for(j = 0; j < global->number_of_unbump_bonds; j++){
      for(i = 0; i < global->number_of_bumps; i++ )
        if(unbump_data[i][j] != NULL ) 
          printf ( " %6.3lf", unbump_data[i][j]->probability ); 
        else printf ( "  *.***" );
      printf ( "\n" );
    }	
#endif
  }

#ifdef TRACE
  printf ( "Angles:\n" );
  for ( j = 0; j < global->number_of_unbump_bonds; j++ )
    {
      for ( i = 0; i < global->number_of_bumps; i++ )
	if ( unbump_data[i][j] != NULL )
	  printf ( " %6.3lf", unbump_data[i][j]->angle );
	else
	  printf ( "  *.***" );
      printf ( "\n" );
    }	  
  for ( i = 0; i < global->number_of_unbump_bonds; i++ )
    if ( global->unbump_bonds[i] >= 0 )
      {      
	printf ( "bond %2d ", i );
	printf ( "(%s,%d): ", 
		 global->target_residues[global->unbump_bonds[i] / 10].num,
		 global->unbump_bonds[i] % 10 );
	for ( j = 0; j < global->number_of_unbump_dependents[i]; j++ )
	  printf ( "%2d ", global->unbump_dependents[i][j] );
	printf ( "\n" );
      }      
    else
      {
	printf ( "bond %2d (%d): ", 
		 i, ( (-1) * global->unbump_bonds[i] ) - 1 );
	for ( j = 0; j < global->number_of_unbump_dependents[i]; j++ )
	  printf ( "%2d ", (-1 ) * (global->unbump_bonds[global->unbump_dependents[i][j]])-1);
	printf ( "\n" );
      }
  for ( i = 0; i < global->ligand->number_of_flexible_bonds; i++ )
    if ( global->ligand->bond_types[i] != RIGID )
      printf ( "%d (%d): %d - %d\n", 
	       i,
	       global->ligand->flexible_bonds[i],
	       global->ligand->bonds[global->ligand->flexible_bonds[i]].atom1 + 1,
	       global->ligand->bonds[global->ligand->flexible_bonds[i]].atom2 + 1 );
#endif
}
