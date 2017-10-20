#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"

void  print_template ( global_data_pt  global )
{
  interaction_pt  points;
  int             i;

  if ( global->template_interactions == NULL
       || global->template_size == 0 )
    err_panic2 ( "print_template", "no template_points");
  points = global->template_interactions;
  for ( i = 0; i < global->template_size; i++ )
    printf ( "%d: %s %8.3f %8.3f %8.3f\n",
	     i, 
	     points[i].act == DONOR       ? "DONOR    " 
	     : points[i].act == ACCEPTOR  ? "ACCEPTOR "  
	     : points[i].act == HYDROPHOB ? "HYDROPHOB" 
	     : "NOTHING  ",
	     points[i].pos[X],
	     points[i].pos[Y],
	     points[i].pos[Z] );
}

void  print_compounds ( global_data_pt  global )
{
  pts_compound_pt  compounds;
  int              i;
  
  compounds = global->compound_interactions->compounds;
  for ( i = 0; i < global->compound_interactions->number_of_compounds; i++ )
    printf ( "%-16s  %3d points   first: %d\n",
	     compounds[i].name,
	     compounds[i].number_of_points,
	     compounds[i].first_point );
  printf ( "%d points\n", global->compound_interactions->number_of_interaction_points );
  for ( i = 0; 
	i < global->compound_interactions->number_of_interaction_points;
	i++ )
    printf ( "%3d: act = %2d, index = %d\n",
	     i,
	     global->compound_interactions->act[i],
	     global->compound_interactions->atom_index[i] );
}

