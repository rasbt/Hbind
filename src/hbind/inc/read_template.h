#ifndef _READ_TEMPLATE_H
#define _READ_TEMPLATE_H
#include "types.h"

int  find_atom_index ( global_data_pt global,
			      char           *residue,
			      char           *residue_num,
			      char           *atom_name );

int  read_template ( char           *filename,
			    global_data_pt global );

#endif
