#include "types.h"

int  assign_number ( molecule_pt  ligand, 
		     int          atom, 
		     int          prohibited,
		     int          number )
{
  int  *neighbors;
  int  i;

  ligand->atoms[atom].number = number;
  neighbors = ligand->neighbors[atom];
  for ( i = 0; i < ligand->number_of_neighbors[atom]; i++ )
    if ( neighbors[i] != prohibited
	 && ligand->atoms[neighbors[i]].number == -1 )
      number = assign_number ( ligand,
			       ligand->neighbors[atom][i],
			       atom,
			       number + 1 );
  return number;
}

void  number_ligand_atoms ( molecule_pt  ligand )
{
  int  i;

  for ( i = 0; i < ligand->number_of_atoms; i++ )
    ligand->atoms[i].number = -1;
  assign_number ( ligand, 0, -1, 0 );
}
