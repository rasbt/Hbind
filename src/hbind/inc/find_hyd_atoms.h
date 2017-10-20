#ifndef _FIND_HYD_ATOMS_H
#define _FIND_HYD_ATOMS_H

extern int  check_hyd_atom ( molecule_pt  molecule,
			     hyd_defn_pt  rules,
			     int          index );

extern void  find_hyd_atoms ( molecule_pt molecule,
			      hyd_defn_pt rules );

#endif
