#include <stdio.h>
#include <string.h>
#include "types.h"
#include "err_handle.h"

int  check_atoms ( molecule_pt  molecule,
		   int          *label,
		   int          *start_counter,
		   int          bond,
		   int          atom,
		   int          back)
{
  int  *neighbors;
  int  result,
       cycles;
  int  i;
  char err_msg[FILENAME_MAX];

  result = 0;
  neighbors = molecule->neighbors[atom];
  if ( label[atom] == UNVISITED )
    {
      /* this is the first time we visit this node */
      label[atom] = VISITED;
      for ( i = 0; i < molecule->number_of_neighbors[atom]; i++ )
	/* visit all neighbors, but don't take direct way back */
	if ( neighbors[i] != back )
	  {
	    bond = molecule->bond_to_neighbor[atom][i];
	    cycles = check_atoms ( molecule,
				   label,
				   start_counter,
				   molecule->bond_to_neighbor[atom][i],
				   neighbors[i],
				   atom);
	    if ( cycles < 0 )
	      /* FATAL_FAILURE is -1 */
	      return FATAL_FAILURE;
	    if ( cycles > 0 )
	      /* the bond to this neighbor is in a cycle */
	      if ( molecule->bonds[bond].type < CYCLE_BOND )
		molecule->bonds[bond].type += CYCLE_BOND;
	    /* add the number of cycles this atom is in due to the 
	       the currten neighbor atom to the overall number of
	       cycles for the current atom */
	    result += cycles;
	  }
      if ( result > 0 )
	{
	  /* this atom is in a cycle */
	  if ( start_counter[atom] > 0 )
	    /* if the cycle has started at this atom, then don't return
	       CYCLE any longer */
	    {
	      result--;
	      start_counter[atom]--;
	    }
	  else
	    label[atom] = CYCLE;
	}
    }
  /* we have been here before, so this is the first atom in the cycle */
  else if ( label[atom] == VISITED ){
    start_counter[atom]++;
    result++;
  }else if ( label[atom] != CYCLE ){
    sprintf ( err_msg, "atom label: %d\n", label[atom] );
    fprintf ( stderr, err_msg);
    fprintf ( stdout, err_msg);
    err_print (err_msg);
    err_warning2 ( "check_atoms", "unknown label");
    return FATAL_FAILURE;
  }
  return result;
}

int  find_cycles ( molecule_pt  molecule)
{
  int     labels[MAX_NUMBER_OF_MOL2_ATOMS],
          start[MAX_NUMBER_OF_MOL2_ATOMS];
  int     i;

  memset(start, 0, molecule->number_of_atoms * sizeof(*start));
  for ( i = 0; i < molecule->number_of_atoms; i++ ) labels[i] = UNVISITED;
  if(check_atoms(molecule, labels, start, -1, 0, -1) == FATAL_FAILURE)
    return FAILURE;
  return SUCCESS;
}
