#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "bitstrings.h"

#define TRAC

void  assign_fragment ( molecule_pt molecule, 
			int         *label,
			int         root_atom,
			int         prohibited_atom,
			int         fragment_index )
{
  int  *neighbors;
  int  i;

  neighbors = molecule->neighbors[root_atom];
  for ( i = 0; i < molecule->number_of_neighbors[root_atom]; i++ )
    if ( neighbors[i] != prohibited_atom 
	 && label[neighbors[i]] != VISITED
	 && molecule->bonds[molecule->bond_to_neighbor[root_atom][i]].type
	 != FLEX )
      {
	molecule->atoms[neighbors[i]].fragment = fragment_index;
	label[neighbors[i]] = VISITED;
	assign_fragment ( molecule,
			  label,
			  neighbors[i],
			  root_atom,
			  fragment_index );
      }
}

/*
 *  Recursive routine to do a depth-first search in the molecule graph
 *  to set a bit in the fragments string of all atoms in the fragment
 *  on one side of the flexible bond with index 'fragment_index'. Root
 *  of the search is 'root_atom', search all paths but not that including
 *  atom with index 'prohibited_atom'.
 */
void  mark_atoms ( molecule_pt molecule, 
		   int         *label,
		   int         root_atom,
		   int         prohibited_atom,
		   int         fragment_index )
{
  int  *neighbors;
  int  i;

  neighbors = molecule->neighbors[root_atom];
  for ( i = 0; i < molecule->number_of_neighbors[root_atom]; i++ )
    if ( neighbors[i] != prohibited_atom && label[neighbors[i]] != VISITED )
      {
	if ( molecule->atoms[neighbors[i]].number 
	     > molecule->atoms[root_atom].number )
	  molecule->atoms[neighbors[i]].fragment = fragment_index + 1;
	label[neighbors[i]] = VISITED;
	mark_atoms ( molecule,
		     label,
		     neighbors[i],
		     root_atom,
		     fragment_index );
	bitstring_set_bit ( molecule->atoms[neighbors[i]].fragments,
			    fragment_index );
      }
}

/*
 *  Routine that marks for all atoms in the corresponding bit string to
 *  which of the two fragments that exist on both sides of a flexible
 *  bond they belong. For that, all atoms behind that adjacent atoms with the
 *  larger index are marked by setting the corresponding bit in the
 *  fragment string.
 */
void  assign_fragments ( molecule_pt  molecule )
{
  bond_pt  bond;
  int      label[MAX_NUMBER_OF_MOL2_ATOMS];
  int      *counter;
  int      fragment_index;
  int      i, j;
  
  /* mark each atom as belonging to fragment no. 0 */
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      molecule->atoms[i].fragment = 0;
      bitstring_clear_all ( molecule->atoms[i].fragments );
    }
      
  for ( i = 0; i < molecule->number_of_flexible_bonds; i++ )
    {
      /* label each atoms as unvisited, this is necessary, because the
	 graph includes cycles, and we only have to visit each atom once */
      for ( j = 0; j < molecule->number_of_atoms; j++ )
	label[j] = UNVISITED;
      bond = &molecule->bonds[molecule->flexible_bonds[i]];
      /* we have already visited the atoms connected to the bond */
      label[bond->atom1] = label[bond->atom2] = VISITED;
      /* the atom with the largest index is the root atom, all atoms in the
	 molecule on this side of the bond are marked */
      if ( molecule->atoms[bond->atom1].number 
	   < molecule->atoms[bond->atom2].number )
      /*
      if ( bond->atom1 < bond->atom2 )
	   */
	{
	  /* mark that atom with index 'bond->atom2' is located "right"
	     from the bond, since its index is larger */
	  bitstring_set_bit ( molecule->atoms[bond->atom2].fragments,
			      i );
	  /* this atom belongs to fragment i+1 */
	  molecule->atoms[bond->atom2].fragment = i + 1;
	  mark_atoms ( molecule,
		       label,
		       bond->atom2,
		       bond->atom1,
		       i );
	}
      else
	{
	  /* mark that atom with index 'bond->atom1' is located "right"
	     from the bond */
	  bitstring_set_bit ( molecule->atoms[bond->atom1].fragments,
			      i );
	  /* this atom belongs to fragment i+1 */
	  molecule->atoms[bond->atom1].fragment = i + 1;
	  mark_atoms ( molecule,
		       label,
		       bond->atom1,
		       bond->atom2,
		       i );
	}
    }
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      molecule->atoms[i].fragment = UNKNOWN;
      label[i] = UNVISITED;
    }
  fragment_index = -1;

  if ( molecule->number_of_flexible_bonds > 0 )
    /* there are molecules without any rotatable bonds, in that case we
       still have to assign fragment numbers */
    for ( i = 0; i < molecule->number_of_flexible_bonds; i++ )
      {
	bond = &molecule->bonds[molecule->flexible_bonds[i]];
	if ( molecule->atoms[bond->atom1].fragment == UNKNOWN )
	  {	
	    label[bond->atom1] = VISITED;
	    fragment_index++;
	    molecule->atoms[bond->atom1].fragment = fragment_index;
	    assign_fragment ( molecule,
			      label,
			      bond->atom1,
			      bond->atom2,
			      fragment_index );
	  }
	if ( molecule->atoms[bond->atom2].fragment == UNKNOWN )
	  {
	    label[bond->atom2] = VISITED;
	    fragment_index++;
	    molecule->atoms[bond->atom2].fragment = fragment_index;
	    assign_fragment ( molecule,
			      label,
			      bond->atom2,
			      bond->atom1,
			      fragment_index );
	  }
      }
  else
    for ( i = 0; i < molecule->number_of_atoms; i++ )
      molecule->atoms[i].fragment = 0;
#ifdef TRACE
  for ( i = 0; i < molecule->number_of_atoms; i++ )
    {
      printf ( "%2d: %d ", i + 1, molecule->atoms[i].fragment );
      bitstring_print_part ( molecule->atoms[i].fragments,
			     molecule->number_of_flexible_bonds );
    }
#endif
  for ( i = 0; i < molecule->number_of_bonds; i++ )
    {
      bond = &molecule->bonds[i];
      bond->fragment1 = molecule->atoms[bond->atom1].fragment;
      bond->fragment2 = molecule->atoms[bond->atom2].fragment;
    }
  for ( i = 0; i < molecule->number_of_fragments; i++ )
    molecule->number_of_fragment_neighbors[i] = 0;
  /* now construct the fragment graph, which consists of all fragments and
     the flexible bonds that connnect them as adjacency matrices */
  for ( i = 0; i < molecule->number_of_fragments; i++ )
    for ( j = 0; j < molecule->number_of_flexible_bonds; j++ )
      {
	bond = &molecule->bonds[molecule->flexible_bonds[j]];
	counter = &molecule->number_of_fragment_neighbors[i];
	if ( molecule->atoms[bond->atom1].fragment == i )
	  {
	    molecule->fragment_neighbors[i][*counter] 
	      = molecule->atoms[bond->atom2].fragment;
	    molecule->bond_to_fragment_neighbor[i][*counter] = j;
	    (*counter)++;
	  }
	else
	  if ( molecule->atoms[bond->atom2].fragment == i )
	    {
	      molecule->fragment_neighbors[i][*counter] 
		= molecule->atoms[bond->atom1].fragment;
	      molecule->bond_to_fragment_neighbor[i][*counter]
		= j;
	      (*counter)++;
	  }
      }	
}

