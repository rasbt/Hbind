#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "types.h"
#include "dist_fun.h"
#include "hash_class.h"
#include "mymalloc.h"
#include "err_handle.h"


#ifdef DEBUG
#include <assert.h>
#endif

typedef struct{
  triangle_parameters_t bounds; 
  int num_perimeter_buckets;
  int num_long_side_buckets;
  int num_short_side_buckets;
}hash_parameters_t;

void alloc_hash_table(global_data_pt global, hash_parameters_t* params);

int calc_table_parameters(global_data_pt global, hash_parameters_t* params);

void compute_combinations(global_data_pt global, hash_parameters_t* params, 
                          int *array, int index, int length);

int create_hash_table(global_data_pt global)
{
  int i, c, p, l, s;
  int triangles = 0, empty = 0, count = 0, max = 0;
  int *array;
  hash_parameters_t params;
  if(calc_table_parameters(global, &params) == FALSE) return FALSE;
  alloc_hash_table(global, &params);

  /* Build the index table*/
  global->template_triangles = 
    mymalloc(global->max_template_triangles * sizeof(triangle_t));
  global->number_of_template_triangles = 0;
  array = mymalloc(global->template_size * sizeof(int));
  for(i = 0; i < global->template_size; i++) array[i] = i;
  for(i = 2; i < global->template_size; i++){
    array[2] = array[i];
    compute_combinations(global, &params, array, 2, 3);
  }
  my_free(array);

  /* Sum up the categories*/
  for(c = 0; c < NUMBER_OF_TYPE_HASH_CLASSES; c++ )
    for(p = 0; p < params.num_perimeter_buckets; p++)
      for(l = 0; l < params.num_long_side_buckets; l++)
	for(s = 0; s < params.num_short_side_buckets; s++)
	  if(global->hash_table[c][p][l][s].number_of_entries > 0 ){
	    triangles += global->hash_table[c][p][l][s].number_of_entries;
	    count++;
	    if(global->hash_table[c][p][l][s].number_of_entries > max )
	      max = global->hash_table[c][p][l][s].number_of_entries;
	  }
	  else empty++;

  printf ( "multi-level hashing:\n" );
  printf ( "%d template triangles\n", global->number_of_template_triangles );
  printf ( "%d triangle pointers\n%d empty buckets\n%d non-empty buckets\n",
	   triangles, empty, count );

  if(triangles != 0)
    printf("%d triangles per bucket on average\n%d maximal number of "
	   "triangles\n", triangles / count, max);
  else{
    printf ("0 triangles identified. Program exiting\n");
    exit (1);
  }

  global->number_of_triangle_pointers = triangles;
  global->max_number_of_triangles = max;
  global->number_of_hash_buckets = count;
  return TRUE;
}

void compute_combinations(global_data_pt global, hash_parameters_t* params,
                          int *array, int index, int length)
{
  float          side[3];
  int            act[3];
  int orig, perimeter_class, longest_class, shortest_class, number_of_triangles;
  double perimeter;
  int i, p, l, s, c;
  int p_start, p_stop, l_start, l_stop, s_start, s_stop;
  int max = 0, min = 0, first = TRUE;
  triangle_parameters_t* bounds = &(params->bounds);
  int key_point;
  triangle_pt triangle;

  interaction_pt points = global->template_interactions;
  side[0] = dist_fun ( points[array[0]].pos, points[array[1]].pos );
  side[1] = dist_fun ( points[array[1]].pos, points[array[2]].pos );
  side[2] = dist_fun ( points[array[2]].pos, points[array[0]].pos );
  perimeter = side[0] + side[1] + side[2];
  for(i = 1; i < 3; i++){
    if(side[i] > side[max]) max = i;
    else if(side[i] < side[min]) min = i;
  }

  if(perimeter >= bounds->min_perimeter  && perimeter < bounds->max_perimeter)
    perimeter_class = 
      (int)(BUCKETS_PER_A_PERIMETER * (perimeter - bounds->min_perimeter));
  else perimeter_class = UNKNOWN;
  if(side[max] >= bounds->min_long_side && side[max] < bounds->max_long_side)
    longest_class = 
      (int)(BUCKETS_PER_A_LENGTH * (side[max] - bounds->min_long_side));
  else longest_class = UNKNOWN;
  if(side[min] >= bounds->min_short_side && side[min] < bounds->max_short_side)
    shortest_class = 
      (int)(BUCKETS_PER_A_LENGTH * (side[min] - bounds->min_short_side));
  else shortest_class = UNKNOWN;
  
  if(global->match_2_key_points)
    key_point = ( points[array[0]].key_point && points[array[1]].key2_point )
      || ( points[array[0]].key_point && points[array[2]].key2_point )
      || ( points[array[1]].key_point && points[array[0]].key2_point )
      || ( points[array[1]].key_point && points[array[2]].key2_point )
      || ( points[array[2]].key_point && points[array[0]].key2_point )
      || ( points[array[2]].key_point && points[array[1]].key2_point );
  else
    key_point = points[array[0]].key_point
      || points[array[1]].key_point
      || points[array[2]].key_point
      || points[array[0]].key2_point
      || points[array[1]].key2_point
      || points[array[2]].key2_point;

  if(perimeter_class != UNKNOWN && longest_class != UNKNOWN
     && shortest_class != UNKNOWN && key_point){
    for(i = 0; i < 3; i++ ) act[i] = points[array[i]].act;
    triangle = NULL;
    c = hash_class[act[0]+act[1]+act[2]];
    p_start = (0 > perimeter_class - 1 ? 0 : perimeter_class - 1);
    p_stop = (params->num_perimeter_buckets < perimeter_class + 2 ? 
      params->num_perimeter_buckets : perimeter_class + 2);
    l_start = (0 > longest_class - 1 ? 0 : longest_class - 1);
    l_stop = (params->num_long_side_buckets < longest_class + 2 ?
      params->num_long_side_buckets : longest_class + 2);
    s_start = (0 > shortest_class - 1 ? 0 : shortest_class - 1);
    s_stop = (params->num_short_side_buckets < shortest_class + 2 ? 
      params->num_short_side_buckets : shortest_class + 2);

    for(p = p_start; p < p_stop; p++)
      for(l = l_start; l < l_stop; l++)
	for(s = s_start; s < s_stop; s++){
          /* The first time we need to create the triangle in the array.
	     All the other times we just point to it.*/
	  if(first){
	    if(global->number_of_template_triangles 
	       >= global->max_template_triangles)
	      err_panic2("compute_combinations","too many template triangles");
	    triangle = 
              &global->template_triangles[global->number_of_template_triangles];
	    for(i = 0; i < 3; i++){
	      triangle->dist[i] = side[i];
	      triangle->index[i] = array[i];
	      triangle->act[i] = act[i];
	    }
	    global->number_of_template_triangles++;
	    first = FALSE;
	  }

	  number_of_triangles = 
	    global->hash_table[c][p][l][s].number_of_entries;
	  if(number_of_triangles == 0)
	    global->hash_table[c][p][l][s].entries =
	      (triangle_pt*) mymalloc(sizeof(triangle_pt));
	  else
	    global->hash_table[c][p][l][s].entries = (triangle_pt*) 
	      myrealloc(global->hash_table[c][p][l][s].entries,
                        (number_of_triangles + 1 ) * sizeof(triangle_pt) );
	  global->hash_table[c][p][l][s].entries[number_of_triangles] = 
	    triangle;
          global->hash_table[c][p][l][s].number_of_entries++;
	}
      }

  if(index > 0){
    orig = array[index-1];    /* store original entry at position index-1 */
    while(array[index-1] < array[index] - 1) {
      /* increase value at position index - 1 as long as it is
	 still smaller than the value at position index */
      array[index-1]++;
      /* recursively increase all elements in the array at
	 positions [0..index-2] */
      compute_combinations(global, params, array, index - 1, length);
    }
    /* restore original value at position index-1 */
    array[index-1] = orig;
  }
}

    

void alloc_hash_table(global_data_pt global, hash_parameters_t* params)
{
  int i, j, k, l;
  int nbuckets;
  hash_table_pt ***hash_table;
  hash_table_pt **perimeter_level_ptr;
  hash_table_pt *long_side_level_ptr;
  hash_table_pt short_side_level_ptr;
                                         
  global->hash_table = (hash_table_pt ***)
    mymalloc(NUMBER_OF_TYPE_HASH_CLASSES * sizeof(hash_table_pt **));
  hash_table = global->hash_table;

  /* Perimeter level*/
  nbuckets = params->num_perimeter_buckets * NUMBER_OF_TYPE_HASH_CLASSES ;
  *hash_table = (hash_table_pt**)mymalloc(nbuckets * sizeof(hash_table_pt*));
  perimeter_level_ptr = *hash_table + params->num_perimeter_buckets;
  for(i = 1; i < NUMBER_OF_TYPE_HASH_CLASSES; i++ ){
    hash_table[i] = perimeter_level_ptr;
    perimeter_level_ptr += params->num_perimeter_buckets;
  }

  /* Long side level*/
  nbuckets *= params->num_long_side_buckets;
  **hash_table = (hash_table_pt*) mymalloc(nbuckets * sizeof(hash_table_pt));
  long_side_level_ptr = **hash_table;
  for(i = 0; i < NUMBER_OF_TYPE_HASH_CLASSES; i++ )
    for(j = 0; j < params->num_perimeter_buckets; j++){
      hash_table[i][j] = long_side_level_ptr;
      long_side_level_ptr += params->num_long_side_buckets;
    }
 
  /* Short side level*/
  nbuckets *= params->num_short_side_buckets;
  ***hash_table = (hash_table_pt) mymalloc(nbuckets * sizeof(hash_table_t));
  short_side_level_ptr = ***hash_table;
  for(i = 0; i < NUMBER_OF_TYPE_HASH_CLASSES; i++ )
    for(j = 0; j < params->num_perimeter_buckets; j++)
      for(k = 0; k < params->num_long_side_buckets; k++){
        hash_table[i][j][k] = short_side_level_ptr;
        short_side_level_ptr += params->num_short_side_buckets;
        for(l = 0; l < params->num_short_side_buckets; l++){
          hash_table[i][j][k][l].entries = NULL;
          hash_table[i][j][k][l].number_of_entries = 0;
        }
      }
}

int calc_table_parameters(global_data_pt global, hash_parameters_t* params)
{
  double upper, lower;
  triangle_parameters_t* tri = &(global->triangle_parameters);
  
  upper = tri->max_perimeter;
  lower = tri->min_perimeter;
  if(upper < lower){
    err_error2("calc_table_parameters", "In the file defs.H the value"
              " for TRIANGLE_MIN_PERIMETER\nmay not be greater than the "
              "value for TRIANGLE_MAX_PERIMETER.\n");
    return FALSE;
  }
  params->num_perimeter_buckets = (int)
    ceil(BUCKETS_PER_A_PERIMETER * (upper - lower)); 
  params->bounds.max_perimeter = upper;
  params->bounds.min_perimeter = lower;

  upper = tri->max_long_side;
  lower = tri->min_long_side;
  if(upper < lower){
    err_error2("calc_table_parameters", "In the file defs.H the value"
              " for TRIANGLE_MIN_LONGEST_SIDE\nmay not be greater than the "
              "value for TRIANGLE_MAX_LONGEST_SIDE.\n");
    return FALSE;
  }
  params->num_long_side_buckets = (int)
    ceil(BUCKETS_PER_A_LENGTH * (upper - lower)); 
  params->bounds.max_long_side = upper;
  params->bounds.min_long_side = lower;

  upper = tri->max_short_side;
  lower = tri->min_short_side;
  if(upper < lower){
    err_error2("calc_table_parameters", "In the file defs.H the value"
              " for TRIANGLE_MIN_SHORTEST_SIDE\nmay not be greater than the "
              "value for TRIANGLE_MAX_SHORTEST_SIDE.\n");
    return FALSE;
  }
  params->num_short_side_buckets = (int)
    ceil(BUCKETS_PER_A_LENGTH * (upper - lower)); 
  params->bounds.max_short_side = upper;
  params->bounds.min_short_side = lower;

  return TRUE;
}
