#include <types.h>
#include "mymalloc.h"


void  sum_charges ( molecule_pt molecule )
{
  int      i, j, k, num_of_non_H_neighbors, neighbor, neighbor2;

  float    *temp_charge_sum = 0;
  int *number_of_neighbors = molecule->number_of_neighbors;
  atom_pt atoms = molecule->atoms;
  /* bond_pt bonds = molecule->bonds; */
  int **neighbors = molecule->neighbors;
  /* int **bond = molecule->bond_to_neighbor; */
  int C_neighbor_only = 0;
  
  temp_charge_sum = 
	(float *) mymalloc ( molecule->number_of_atoms * sizeof (float) );
  
  /* 1. Sum charges of all H neighbors to heavy atoms (non-H atoms) and record
   * the new charges in "charge_sum" entry.
   */
  for ( i = 0; i < molecule->number_of_atoms; i++ ){
    temp_charge_sum[i] = atoms[i].charge;
    if(atoms[i].type == H) continue;

#ifdef TRACE1
    printf("atoms[%d]  name = %s  charge = %6.4f\n", i, atoms[i].name, 
           atoms[i].charge);
#endif

    for ( j = 0; j < number_of_neighbors[i]; j++ ){
      neighbor = neighbors[i][j];
      if(atoms[neighbor].type == H){
        temp_charge_sum[i] += atoms[neighbor].charge;

#ifdef TRACE1	
        printf("- neighbor atoms[%d]  name = %s  charge = %6.4f", neighbor, 
               atoms[neighbor].name, atoms[neighbor].charge);
        printf("  temp_charge_sum %6.4f += %6.4f\n", temp_charge_sum[i], 
               atoms[neighbor].charge);
#endif			
      }
    }
  }
  
#ifdef TRACE
  printf ( "after suming all Hs charges: \n" );
  for( i = 0; i < molecule->number_of_atoms; i++ ){
    if(atoms[i].type != H) 
      printf("atoms[%d]  name = %s  charge = %6.4f  charge_sum = %6.4f\n", 
             i, atoms[i].name, atoms[i].charge, atoms[i].charge_sum);
  }
#endif
  
  /* 2. Sum charges of all non-H neighbors to atom N or O. If they only have C
   * neighbors and these C neighbors have other N or O neighbors, sum also 
   * the charges of their N or O neighbors.
   */
  for(i = 0; i < molecule->number_of_atoms; i++ ) {
    atoms[i].charge_sum = temp_charge_sum[i];
    if(atoms[i].type != N && atoms[i].type != O) continue;

#ifdef TRACE2
    printf("atoms[%d]  name = %s  charge_sum = %6.4f\n", i, atoms[i].name, 
           atoms[i].charge_sum);
#endif
    num_of_non_H_neighbors = 0;
    neighbor = -1;
    C_neighbor_only = 1;
    for ( j = 0; j < number_of_neighbors[i]; j++ ) {
      if(atoms[neighbors[i][j]].type == H) continue;

      neighbor = neighbors[i][j];
      atoms[i].charge_sum += temp_charge_sum[neighbor];
      num_of_non_H_neighbors++;
      if(atoms[neighbor].type != C) C_neighbor_only = 0;

#ifdef TRACE2		
      printf("-> neighbor atoms[%d]  name = %s  charge_sum = %6.4f", 
             neighbor, atoms[neighbor].name, atoms[neighbor].charge_sum);	
      printf("  charge_sum %6.4f += %6.4f\n", atoms[i].charge_sum, 
             temp_charge_sum[neighbor]);
#endif
    }

    /* If there is only 1 direct non_H neighbor, we need to sum up charges on
     * the 2nd layer neighbors 
     */
    if(num_of_non_H_neighbors == 1){
      for(j = 0; j < number_of_neighbors[neighbor]; j++){ 
        neighbor2 = neighbors[neighbor][j]; 
        if(atoms[neighbor2].type != H && neighbor2 != i){ 
          atoms[i].charge_sum += temp_charge_sum[neighbor2]; 
          num_of_non_H_neighbors ++;

#ifdef TRACE2		
          printf("  => neighbor atoms[%d]  name = %s  charge_sum = %6.4f", 
                 neighbor2, atoms[neighbor2].name, atoms[neighbor2].charge_sum);
          printf("  charge_sum %6.4f += %6.4f\n", atoms[i].charge_sum, 
                 temp_charge_sum[neighbor2]);
#endif			
        } 
      } 

    /* If all neighbors to atom i are carbons, we need to sum charges of any 
     * other N/O neighbor to these C neighbors 
     *          l----Carbon 
     * atom i---l---Carbon 
     *          l---Carbon------N/O 
     */
    }else if(C_neighbor_only == 1){ 
      for(j = 0; j < number_of_neighbors[i]; j++){ 
        neighbor = neighbors[i][j]; 
        if(atoms[neighbor].type != C) continue;

        for(k = 0; k < number_of_neighbors[neighbor]; k++){ 
          if(neighbors[neighbor][k] != i && 
             (atoms[neighbors[neighbor][k]].type == N || 
              atoms[neighbors[neighbor][k]].type == O )){
            atoms[i].charge_sum += temp_charge_sum[neighbors[neighbor][k]];
#ifdef TRACE2		
          printf("name = %s, temp_charge = %6.4f\n", 
                 atoms[neighbors[neighbor][k]].name,
                 temp_charge_sum[neighbors[neighbor][k]]);
          printf("  ** charge_sum %6.4f += %6.4f\n", atoms[i].charge_sum, 
                 temp_charge_sum[neighbors[neighbor][k]]);
#endif			
          }
        } 
      } 
    } 
  } 
   
#ifdef TRACE
  printf ( "after suming all neighbor charges: \n" );
  for ( i = 0; i < molecule->number_of_atoms; i++ ) {
     if ( atoms[i].type != H && atoms[i].type != C ) 
	printf ("atoms[%d]  name = %s  charge = %6.4f  charge_sum = %6.4f\n", 
		i, atoms[i].name, atoms[i].charge, atoms[i].charge_sum);
  }
#endif

  /* Now we want to remove those small charges and charges for non-heavy atoms */
  for(i = 0; i < molecule->number_of_atoms; i++ ) { 
    if ( atoms[i].type == H || atoms[i].type == C || atoms[i].type == S ) 
      atoms[i].charge_sum = 0.0; 
    else if(-MINIMAL_CHARGE < atoms[i].charge_sum && 
            atoms[i].charge_sum < MINIMAL_CHARGE) 
      atoms[i].charge_sum = 0.0;
  }

#ifdef TRACE
  printf ( "after removing small charges: \n" );
  for ( i = 0; i < molecule->number_of_atoms; i++ ) {
      if (atoms[i].charge_sum != 0.0 )
	printf ("atoms[%d]  name = %s  charge = %6.4f  charge_sum = %6.4f\n", 
		i, atoms[i].name, atoms[i].charge, atoms[i].charge_sum);
  }
#endif

  if(temp_charge_sum) free(temp_charge_sum);
}
