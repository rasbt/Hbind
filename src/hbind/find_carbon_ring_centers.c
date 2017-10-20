#include "defs.h"
#include "types.h"

/*
 *  This recursive routine does a depth-first search in the subgraph of 
 *  the molecule graph consisting only of carbons to identify cycles with
 *  up to MAX_CARBON_RING_SIZE atom. Once a ring has been found, the
 *  recursion immediately steps back and stores all node indices on the
 *  way back to the first node in the cycle in the array cycle.
 */
int  check_for_ring ( molecule_pt  molecule,
		      int          *mark,
		      int          *cycle,
		      int          index,
		      int          prohibited,
		      int          count )
{
  int  *neighbors;
  int  result,
       old_mark;
  int  i;

  if ( mark[index] < -1 )
    /* no carbon, so don't go any further */
    return NO_CYCLE;
  if ( count <= MAX_CARBON_RING_SIZE )
    /* go on only if we haven't visited MAX_CARBON_RING_SIZE atoms yet */
    {
      /* mark atom 'index' as atom number 'count' in the potential ring,
	 but save original mark */
      old_mark = mark[index];
      mark[index] = count;
      neighbors = molecule->neighbors[index];
      result = NO_CYCLE;
      for ( i = 0; i < molecule->number_of_neighbors[index]; i++ )
	/* check all neigbors */
	if ( mark[neighbors[i]] > -2 
	     && mark[neighbors[i]] <= 1 
	     && neighbors[i] != prohibited )
	  /* consider only carbons, if you bump into an atom of the current
	     ring other than the first one, stop the search, and don't go 
	     directly back to where you came from */
	  {
	    if ( mark[neighbors[i]] == 1 )
	      /* we have found the first node in our cycle */
	      result = CYCLE;
	    else
	      /* check the next neighbor */
	      result = check_for_ring ( molecule,
					mark,
					cycle,
					neighbors[i],
					index,
					count + 1 );
	    if ( result == CYCLE )
	      {
		/* store index of this ring atom in array 'cycle' */
		cycle[count-1] = index;
		/* this is to mark that this atom is in a cycle that has been
		   already found, i.e. no new search can start at this atom
		   to prevent marking the same cycle twice or more */
		mark[index] = -1;
		return CYCLE;
	      }
	  }
      /* this atom is not part of a cycle, so give back it's original mark */
      mark[index] = old_mark;
    }
  return NO_CYCLE;
}

/*  This routine searches all carbon rings with up to MAX_CARBON_RING_SIZE
 *  atoms. It is not very efficently programmed, but since it is part of
 *  the global preprocessing when initializing the screening database, it 
 *  has to be only done once per database. Markings for all atoms are
 *  stored in an array 'mark', the entries have the following meanings:
 *   -2   no carbon, ignore this atom during the search
 *   -1   carbon already included in a cycle
 *   0    an unvisited carbon or carbon not included in a cycle yet
 *  The reason for distinguishing between the last two cases is that
 *  a search will never start with an atom that is already included in 
 *  another ring to prevent identifying the same ring more than once.
 */
void  find_carbon_ring_centers ( molecule_pt  molecule )
{
  atom_pt  atoms;
  float    avg[3];
  int      mark[MAX_NUMBER_OF_MOL2_ATOMS],
           cycle[MAX_CARBON_RING_SIZE],
           *neighbors;
  int      count, 
           number_of_carbon_rings;
  int      i, j;

  atoms = molecule->atoms;
  number_of_carbon_rings = 0;
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    /* mark all carbons */
    if ( atoms[i].type == C )
      mark[i] = 0;
    else
      mark[i] = -2;
  for ( i = 0; i < MAX_CARBON_RING_SIZE; i++ )
    cycle[i] = -1;
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      if ( mark[i] == 0 )
	/* atom i is a carbon that is not included in any of the rings
	   found yet */
	{
	  /* mark it as the first atom in the potential cycle */
	  mark[i] = 1;
	  neighbors = molecule->neighbors[i];
	  for ( j = 0; 
		j < molecule->number_of_neighbors[i]; 
		j++ )
	    if ( mark[neighbors[j]] != -2 )
	      /* neigbor is carbon */
	      if ( check_for_ring ( molecule,
				    mark,
				    cycle,
				    neighbors[j],
				    i,
				    2 ) == CYCLE )
		/* a ring starting at atom 'i' has been found */
		{
		  cycle[0] = i;
		  /* do not check the other neighbors, we already found
		     one ring */
		  break;
		}
	  /* set mark back to non-ring carbon */
	  mark[i] = 0;
	}
      if ( cycle[0] >= 0 )
	/* a cycle has been found */
	{
	  /* ok, mark it as a ring carbon */
	  mark[i] = -1;
	  /* compute the center of this ring as the average of the
	     positions of it's atoms */
	  avg[X] = avg[Y] = avg[Z] = 0.0;
	  count = 0;
	  for ( j = 0; j < MAX_CARBON_RING_SIZE; j++ )
	    if ( cycle[j] >= 0 )
	      {
		count++;
		avg[X] += atoms[cycle[j]].pos[X];
		avg[Y] += atoms[cycle[j]].pos[Y];
		avg[Z] += atoms[cycle[j]].pos[Z];
		/*
		  printf ( "%2d: %3s %7.3f   %7.3f   %7.3f\n",
			 cycle[j],
			 atoms[cycle[j]].name,
			 atoms[cycle[j]].pos[X],
			 atoms[cycle[j]].pos[Y],
			 atoms[cycle[j]].pos[Z] );
			 */
		cycle[j] = -1;
	      }
	  molecule->carbon_ring_centers[number_of_carbon_rings][X] =
	    avg[X] / count;
	  molecule->carbon_ring_centers[number_of_carbon_rings][Y] =
	    avg[Y] / count;
	  molecule->carbon_ring_centers[number_of_carbon_rings][Z] =
	    avg[Z] / count;
	  number_of_carbon_rings++;
	}
    }
  molecule->number_of_carbon_rings = number_of_carbon_rings;
}
