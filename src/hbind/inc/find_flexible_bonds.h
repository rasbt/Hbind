#ifndef _FIND_FLEXIBLE_BONDS_H
#define _FIND_FLEXIBLE_BONDS_H
#include "types.h"

int find_flexible_bonds(molecule_pt molecule, flex_bond_defn_pt flex_bond_rules,
                        int number_of_rules);

#endif
