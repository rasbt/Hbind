#include <types.h>

void  construct_adjacency_list ( molecule_pt molecule )
{
  bond_pt  bonds;
  int      **neighbors;
  int      **bond;
  int      *number_of_neighbors;
  int      *bond_order;
  int      i, atom1, atom2;
  
  number_of_neighbors = molecule->number_of_neighbors;
  bond_order = molecule->bond_order;
  for(i = 0; i < molecule->number_of_atoms; i++){
    number_of_neighbors[i] = 0;
    bond_order[i] = 0;
  }
  bonds = molecule->bonds;
  neighbors = molecule->neighbors;
  number_of_neighbors = molecule->number_of_neighbors;
  bond_order = molecule->bond_order;
  bond = molecule->bond_to_neighbor;
  for(i = 0; i < molecule->number_of_bonds; i++){
    atom1 = bonds[i].atom1;
    atom2 = bonds[i].atom2;
    neighbors[atom1][number_of_neighbors[atom1]] = atom2;
    bond[atom1][number_of_neighbors[atom1]++] = i;
    neighbors[atom2][number_of_neighbors[atom2]] = atom1;
    bond[atom2][number_of_neighbors[atom2]++] = i;
      
    if ( molecule->bonds[i].type == SINGLE ) {
      bond_order[atom1] +=  1;
      bond_order[atom2] +=  1;      
    } else if ( molecule->bonds[i].type == DOUBLE ) {
      bond_order[atom1] +=  2;
      bond_order[atom2] +=  2;      
    } else if ( molecule->bonds[i].type == TRIPLE ) {
      bond_order[atom1] +=  3;
      bond_order[atom2] +=  3;      
    } else if ( molecule->bonds[i].type == QUADRUPLE ) {
      bond_order[atom1] +=  4;
      bond_order[atom2] +=  4;      
    } else if ( molecule->bonds[i].type == DELOCALIZED ) {
      bond_order[atom1] +=  1;
      bond_order[atom2] +=  1;      
    } else if ( molecule->bonds[i].type == AMIDE ) {
      bond_order[atom1] +=  1;
      bond_order[atom2] +=  1;      
    } else if ( molecule->bonds[i].type == AROMATIC ) {
      bond_order[atom1] +=  1.5;
      bond_order[atom2] +=  1.5;      
    }
  }
}
