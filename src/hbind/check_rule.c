#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "check_rule.h"

int  check_prohibited_rule ( molecule_pt  molecule,
			     int          base_atom,
			     int          bond_atom,
			     rule_list_pt rule )
{
  int      *neighbors;
  atom_pt  atom;
  int      result,
           number_of_neighbors,
           number_of_matches;
  int      i;

  neighbors = molecule->neighbors[base_atom];
  number_of_neighbors = molecule->number_of_neighbors[base_atom];
  number_of_matches = 0;
  for ( i = 0; i < number_of_neighbors; i++ )
    if ( neighbors[i] != bond_atom )
      /* don't consider atom on the other end of the bond we are checking */
      {
	atom = &molecule->atoms[neighbors[i]];
	result = check_atom ( rule->type,
			      rule->orbit,
			      atom->type, 
			      atom->orbit );
	if ( result == SUCCESS )
	  {
	    /*	    printf("prohibited match!\n");*/
	    number_of_matches++;
	    /* atom matches, so check, if we have counted enough of this type,
	     if so, since this is a prohibited rule, return FAILURE */
	    if ( number_of_matches == rule->number )
	      {
		/*		printf("FAILURE!!!!\n");*/
		return FAILURE;
	      }
	    /* now check, if this is a nested rule */
	    if ( rule->prohibited != NULL )
		result = check_prohibited_rule ( molecule,
						 neighbors[i],
						 base_atom,
						 rule->prohibited );
	    if ( result == FAILURE )
	      return FAILURE;
	    if ( rule->required != NULL )
	      result = check_required_rule ( molecule,
					     neighbors[i],
					     base_atom,
					     rule->required );
	    if ( result == FAILURE )
	      return FAILURE;
	  }
      }
  /*  printf("SUCCESS!!!\n");*/
  /* return SUCCESS, since rule is fulfilled (or satisfied?) */
  return SUCCESS;
}

int  check_required_rule ( molecule_pt  molecule,
			   int          base_atom,
			   int          bond_atom,
			   rule_list_pt rule )
{
  int      *neighbors;
  atom_pt  atom;
  int      result,
           number_of_neighbors,
           number_of_matches;
  int      i;

  neighbors = molecule->neighbors[base_atom];
  number_of_neighbors = molecule->number_of_neighbors[base_atom];
  number_of_matches = 0;
  for ( i = 0; i < number_of_neighbors; i++ )
    if ( neighbors[i] != bond_atom )
      /* don't consider atom on the other end of the bond we are checking */
      {
	atom = &molecule->atoms[neighbors[i]];
	result = check_atom ( rule->type,
			      rule->orbit,
			      atom->type, 
			      atom->orbit );
	if ( result == SUCCESS )
	  {
	    number_of_matches++;
	    if ( rule->prohibited != NULL )
	      result = check_prohibited_rule ( molecule,
					       neighbors[i],
					       base_atom,
					       rule->prohibited );
	    if ( result == FAILURE )
	      return FAILURE;
	    if ( rule->required != NULL )
	      result = check_required_rule ( molecule,
					     neighbors[i],
					     base_atom,
					     rule->required );
	    if ( result == FAILURE )
	      return FAILURE;
	  }	    
    }
  /* check, if we have found enough atoms of the required type, if so, the
     rule is fulfilled, so return SUCCESS */
  if ( number_of_matches >= rule->number )
    return SUCCESS;
  return FAILURE;
}

