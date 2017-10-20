#ifndef _CHECK_RULE_H
#define _CHECK_RULE_H

static inline int 
check_atom(int rule_atom, int rule_orbit, int atom, int orbit)
{
  if((rule_atom == ANY || rule_atom == atom ) && 
     (rule_orbit == ANY || rule_orbit == orbit)) return SUCCESS;
  return FAILURE;
}
#if 0
extern int  check_atom ( int  rule_atom,
			 int  rule_orbit,
			 int  atom,
			 int  orbit );
#endif

extern int  check_prohibited_rule ( molecule_pt  molecule,
				    int          base_atom,
				    int          bond_atom,
				    rule_list_pt rule );

extern int  check_required_rule ( molecule_pt  molecule,
				  int          base_atom,
				  int          bond_atom,
				  rule_list_pt rule );

#endif
