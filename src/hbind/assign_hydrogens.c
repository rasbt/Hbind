/*
 *
 *    assign_hydrogens.c   Volker Schnecke     Mon Feb 16 21:17:11 EST 1998
 *
 *    function: assign_hydrogens()
 *
 *    assigns hydrogens to those nitrogens and oxygens where they seem to 
 *    be missing in the structure (adds only atoms but no positions)
 *
 */
#include <stdio.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"

/*
 *  This routine checks all nitrogens and oxygens in 'molecule' by counting 
 *  the number of bonds (considering double bonds as two bonds, triple as 
 *  three) and adding hydrogens, if less than three (two for oxygen) bonds 
 *  are found.
 */
int  assign_hydrogens ( molecule_pt  molecule)
{
  int i;
  char errmsg[FILENAME_MAX];
  int total_bonds = 0;
  int number_of_atoms = molecule->number_of_atoms;

  molecule->number_of_added_hydrogens = 0;
  for ( i = 0; i < number_of_atoms; i++ )
    if ( molecule->atoms[i].type == N || molecule->atoms[i].type == O )
      /* only check nitrogens and oxygens */
      {
	if ( molecule->atoms[i].type == N )
	  {
	    if ( (molecule->atoms[i].orbit == SP3) || (molecule->atoms[i].orbit == SP2) || (molecule->atoms[i].orbit == SP1) || (molecule->atoms[i].orbit == AR) || (molecule->atoms[i].orbit == AM))
	      {
		/* printf("orbit %d %d\n", i, molecule->atoms[i].orbit);*/
		/* nitrogen should have three bonds, otherwise it is assumed
		   that hydrogens are not included in the structure */
		total_bonds = 3;
	      }
	    else
	      {
		total_bonds = 4;
	      }
	  }
	else 
	  {
	    if ( ( molecule->atoms[i].orbit == SP3 ) || ( molecule->atoms[i].orbit == SP2 ) )
	      {
		/* oxygens should have two bonds */
		/*printf("orbit %d %d\n", i, molecule->atoms[i].orbit);*/
		total_bonds = 2;
	      }
	  }
	/*	if ( molecule->number_of_neighbors[i] != total_bonds )  toneroma 09FEB07 added bond_order check and used it instead, since it gives a better indication of valency */
	if ((  molecule->bond_order[i] != total_bonds ) && (( molecule->atoms[i].type != O ) && ( molecule->atoms[i].orbit != CO2 )) && (( molecule->atoms[i].type != N ) && ( molecule->atoms[i].orbit != PL3 )))
	  {
	    /*	    printf ( "%d %d %s %d %d, #bonds %s BO: %d\n", i, molecule->atoms[i].type, molecule->atoms[i].type_str,  molecule->number_of_neighbors[i], total_bonds, molecule->bonds[i].type_str, molecule->bond_order[i] ); */
	    sprintf (errmsg, "mol2 file incorrectly protonated. Protonate and re-analyze ligand %s\n", molecule->name );
	    err_error2 ( "assign_hydrogens", errmsg);
	    return FAILURE;
	  }
	else if ((  molecule->bond_order[i] > 2 ) && (( molecule->atoms[i].type == O ) && ( molecule->atoms[i].orbit == CO2 )))
	  {
	    /*	    printf ( "%d %d %s %d %d, #bonds %s BO: %d\n", i, molecule->atoms[i].type, molecule->atoms[i].type_str,  molecule->number_of_neighbors[i], total_bonds, molecule->bonds[i].type_str, molecule->bond_order[i] ); */
	    sprintf (errmsg, "mol2 file incorrectly protonated. Protonate and re-analyze ligand %s\n", molecule->name );
	    err_error2 ( "assign_hydrogens", errmsg);
	    return FAILURE;
	  }
	else if ((  molecule->number_of_neighbors[i] != 3 ) && (( molecule->atoms[i].type == N ) && ( molecule->atoms[i].orbit == PL3 )))
	  {
	    /*	    printf ( "%d %d %s %d %d, #bonds %s BO: %d\n", i, molecule->atoms[i].type, molecule->atoms[i].type_str,  molecule->number_of_neighbors[i], total_bonds, molecule->bonds[i].type_str, molecule->bond_order[i] ); */
	    sprintf (errmsg, "mol2 file incorrectly protonated. Protonate and re-analyze ligand %s\n", molecule->name );
	    err_error2 ( "assign_hydrogens", errmsg);
	    return FAILURE;
 	  }
      }
  return SUCCESS;
}
