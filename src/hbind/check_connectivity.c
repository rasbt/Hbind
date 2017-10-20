#include <stdio.h>
#include "types.h"
#include "defs.h"

void  label_atoms ( molecule_pt  molecule,
		    int          *label,
		    int          atom,
		    int          back )
{
  int  *neighbors;
  int  i;

  neighbors = molecule->neighbors[atom];
  if ( label[atom] == UNVISITED )
    {
      /* this is the first time we visit this node */
      label[atom] = VISITED;
      for ( i = 0; i < molecule->number_of_neighbors[atom]; i++ )
	/* visit all neighbors, but don't take direct way back */
	if ( neighbors[i] != back )
	  {
	    label_atoms ( molecule,
			  label,
			  neighbors[i],
			  atom );
	  }
    }
}
  
int  check_connectivity ( molecule_pt  molecule )
{
  int     labels[MAX_NUMBER_OF_MOL2_ATOMS];
  int     i;

  for ( i = 0; i < molecule->number_of_atoms; i++ )
    labels[i] = UNVISITED;
  label_atoms ( molecule, labels, 0, -1 );
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    if ( labels[i] == UNVISITED )
	return FAILURE;
  return SUCCESS;
}
