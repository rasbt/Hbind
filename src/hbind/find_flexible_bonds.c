#include <stdio.h>
#include <string.h>
#include "types.h"
#include "check_rule.h"
#include "err_handle.h"

/*  This function checks, if the bond connecting 'atom1' and 'atom2'
 *  might be a flexible bond based on the rules for the adjacent
 *  atoms to a flexible bond. If the atoms match, the number returned
 *  depends on the order: if atom1 matches atom[0] in the rule and
 *  atom2 matches atom[1], then rule_number+1 is returned, if atom1
 *  matches atom[1] and atom2 matches atom[0] then -1*(rule_number+1)
 *  is returned, if the bond is not rotatable, then 0 is returned.
 */
int  check_bond_atoms ( flex_bond_defn_pt  rules,
			int                number_of_rules,
			atom_pt            atom1,
			atom_pt            atom2 )
{
  int  rule_type1,
       rule_orbit1,
       rule_type2,
       rule_orbit2;
  int  i;
  
  i = number_of_rules;
  /*  for ( i = 0; i < number_of_rules; i++ )
      {*/
      rule_type1 = rules[i].atom[0].type;
      rule_orbit1 = rules[i].atom[0].orbit;
      rule_type2 = rules[i].atom[1].type;
      rule_orbit2 = rules[i].atom[1].orbit;
      if ( ( rule_type1 == atom1->type
	     || rule_type1 == ANY )
	   && ( rule_orbit1 == atom1->orbit
		|| rule_orbit1 == ANY )
	   && ( rule_type2 == atom2->type
		|| rule_type2 == ANY )
	   && ( rule_orbit2 == atom2->orbit
		|| rule_orbit2 == ANY ) )
	return i+1;
      if ( ( rule_type2 == atom1->type
	     || rule_type2 == ANY )
	   && ( rule_orbit2 == atom1->orbit
		|| rule_orbit2 == ANY )
	   && ( rule_type1 == atom2->type
		|| rule_type1 == ANY )
	   && ( rule_orbit1 == atom2->orbit
		|| rule_orbit1 == ANY ) )
	return (-1) * (i+1);

  return 0;
}

int find_flexible_bonds(molecule_pt molecule, flex_bond_defn_pt flex_bond_rules,
                        int number_of_rules)
{
  rule_pt  bond1,
           bond2;
  int      rule,
           atom1,
           atom2,
           result1,
           result2,
           bond_flex;
  int      i, j, k;
  int      rule_used;
  
  
  molecule->number_of_flexible_bonds = 0;
  for( i = 0; i < molecule->number_of_bonds; i++ ){
    if( molecule->bonds[i].type != SINGLE) continue;

    rule_used = 99;
    /* check if the atoms match one of the rules */
    for ( k = 0; k < number_of_rules; k++ ){

      rule = check_bond_atoms(flex_bond_rules, k,
                              &molecule->atoms[molecule->bonds[i].atom1],
                              &molecule->atoms[molecule->bonds[i].atom2] );
      if ( rule == 0) continue;
      if ( rule > 0 ){ 
        atom1 = molecule->bonds[i].atom1;
        atom2 = molecule->bonds[i].atom2;
      }else{
        atom1 = molecule->bonds[i].atom2;
        atom2 = molecule->bonds[i].atom1;
        rule *= -1;
      }
	      
      result1 = SUCCESS;
      result2 = SUCCESS;
      bond1 = &flex_bond_rules[rule-1].atom[0];
      bond2 = &flex_bond_rules[rule-1].atom[1];
      bond_flex = flex_bond_rules[rule-1].flex;
      for(j = 0; j < MAX_RULES_PER_BONDED_ATOM && result1 == SUCCESS; j++ )
        if ( bond1->prohibited[j] != NULL )
          result1 = check_prohibited_rule(molecule, atom1, atom2, 
                                          bond1->prohibited[j] );
      for(j = 0; j < MAX_RULES_PER_BONDED_ATOM && result1 == SUCCESS; j++ )
        if(bond1->required[j] != NULL )
          result1 = check_required_rule(molecule, atom1, atom2, 
                                        bond1->required[j] );
      for ( j = 0; j < MAX_RULES_PER_BONDED_ATOM && result2 == SUCCESS; j++ )
        if ( bond2->prohibited[j] != NULL )
          result2 = check_prohibited_rule(molecule, atom2, atom1,
                                          bond2->prohibited[j] );
      for ( j = 0; j < MAX_RULES_PER_BONDED_ATOM && result2 == SUCCESS; j++ )
        if ( bond2->required[j] != NULL ) 
          result2 = check_required_rule(molecule, atom2, atom1, 
                                        bond2->required[j] );

      if ( result1 == SUCCESS && result2 == SUCCESS ){
        rule_used = rule;
        if (bond_flex == FLEX ) molecule->bonds[i].type = FLEX;
        if (bond_flex == RIGID ) molecule->bonds[i].type = RIGID;
      }
    }
    if(molecule->bonds[i].type == FLEX)
      molecule->flexible_bonds[molecule->number_of_flexible_bonds++] = i;

    
    if(molecule->number_of_flexible_bonds > MAX_NUMBER_OF_FLEXIBLE_BONDS ){
      err_warning2("find_flexible_bonds",
                   "more than MAX_NUMBER_OF_FLEXIBLE_BONDS flexible bonds");
      return FAILURE;
    }
  }
  molecule->number_of_fragments = molecule->number_of_flexible_bonds + 1;
  return SUCCESS;
}
