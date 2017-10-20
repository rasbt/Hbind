#include <string.h>
#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "match_triangles.h"
#include "read_hyd_defn.h"
#include "check_connectivity.h"
#include "assign_hydrogens.h"
#include "read_mol2.h"
#include "err_handle.h"
#include "find_cycles.h"
#include "adj_list.h"
#include "find_flexible_bonds.h"
#include "find_carbon_ring_centers.h"
#include "find_hyd_atoms.h"
#include "sum_charges.h"

int  get_compounds_from_args ( global_data_pt  global,
			       char            **argv,
			       int             argc )
{
  pts_compound_pt compounds;
  atom_pt         atom;
  char            file[FILENAME_MAX];
  char            *help;
  int             counter,
                  compound_counter;
  int             i, j;
  char            err_msg[FILENAME_MAX];
  FILE            *lig_file;

  printf ( "\ngetting compounds to screen from command-line arguments\n\n" );
  counter = 0;
  compounds = global->compound_interactions->compounds;
  for( i = 4; i < argc; i++ ){
    compound_counter = 0;      
    compounds[i-4].first_point = counter;

    if((lig_file = open_mol2(argv[i])) == NULL) return FATAL_FAILURE;
    if(read_mol2(lig_file, argv[i], global, NULL) != SUCCESS){
      fclose(lig_file);
      return FATAL_FAILURE;
    }
    fclose(lig_file);

    sum_charges ( global->ligand );     
    strcpy ( file, argv[i] );
    help = file + strlen ( file );
    while ( help != file && *help != '/' ) help--;
    if ( help != file ) {
      *help = '\0';
      help++;
      strcpy ( global->compound_dir, file );
    }else sprintf ( global->compound_dir, "." );

      help[strlen(help)-5] = '\0';
      construct_adjacency_list ( global->ligand );
      if ( check_connectivity ( global->ligand ) == FAILURE ) {
	  sprintf ( err_msg, "ERROR: unconnected atoms in file %s\n", argv[i] );
	  fprintf ( stderr, err_msg);
	  err_print(err_msg);
	  err_warning2 ( "screen_single_compounds", "skipping ligand");
	}
      find_carbon_ring_centers ( global->ligand );
      assign_hydrogens ( global->ligand);
      find_cycles ( global->ligand);
      find_flexible_bonds ( global->ligand, global->flex_bond_rules,
			    global->number_of_flex_bond_rules);
      find_hyd_atoms ( global->ligand, &global->hyd_atom_rules );
      for ( j = 0; j < global->ligand->number_of_atoms; j++ )
	if ( global->ligand->atoms[j].act != NOTHING ) {
	    atom = &global->ligand->atoms[j];
	    global->compound_interactions->act[counter] = atom->act;
	    global->compound_interactions->atom_index[counter] = j;
	    global->compound_interactions->pos[3*counter+X] = atom->pos[X];
	    global->compound_interactions->pos[3*counter+Y] = atom->pos[Y];
	    global->compound_interactions->pos[3*counter+Z] = atom->pos[Z];
	    counter++;
	    compound_counter++;
	  }
      for ( j = 0; j < global->ligand->number_of_carbon_rings; j++ ) {
	  global->compound_interactions->act[counter] = HYDROPHOB;
	  global->compound_interactions->atom_index[counter] 
	    = global->ligand->carbon_ring_atom[j];
	  global->compound_interactions->pos[3*counter+X] 
	    = global->ligand->carbon_ring_centers[j][X];
	  global->compound_interactions->pos[3*counter+Y] 
	    = global->ligand->carbon_ring_centers[j][Y];
	  global->compound_interactions->pos[3*counter+Z] 
	    = global->ligand->carbon_ring_centers[j][Z];
	  counter++;
	  compound_counter++;
	}
      compounds[i-4].number_of_points = compound_counter;
      strcpy ( compounds[i-4].name, help );
      printf ( "Name of compound %d: %s\n", i - 4, help );
    }
  return argc - 4;
}

int  get_compounds_from_file ( global_data_pt  global,
			       FILE            *fp )
{
  pts_compound_pt compounds;
  atom_pt         atom;
  char            file[FILENAME_MAX],
                  line[256];
  char            *help;
  int             counter,
                  compound_counter,
                  number_of_compounds;
  int             j;
  char            err_msg[FILENAME_MAX];
  FILE            *lig_file;

  printf ( "\ngetting compounds to screen from file\n\n" );
  /* reset the stream, since we already read the first line of the file to
     figure out, if it is a database definition file */
  rewind ( fp );
  counter = 0;
  compound_counter = 0;      
  number_of_compounds = 0;
  compounds = global->compound_interactions->compounds;
  while(fgets ( line, sizeof ( line ), fp ) )
    if(line[0] != ' ' && line[0] != '\0' && line[0] != '\n' ){ 
      sscanf ( line, "%s", file );
      compounds[number_of_compounds].first_point = counter;

      if((lig_file = open_mol2(file)) == NULL) return FATAL_FAILURE;
      if(read_mol2(lig_file, file, global, NULL) != SUCCESS){
        fclose(lig_file);
        return FATAL_FAILURE;
      }
      fclose(lig_file);

      sum_charges ( global->ligand );     
	help = file + strlen ( file );
	while ( help != file && *help != '/' )
	  help--;
	if ( help != file )
	  {
	    *help = '\0';
	    help++;
	    strcpy ( global->compound_dir, file );
	  }
	else
	  sprintf ( global->compound_dir, "." );
	help[strlen(help)-5] = '\0';
	construct_adjacency_list ( global->ligand );
	if ( check_connectivity ( global->ligand ) == FAILURE )
	  {
	    sprintf ( err_msg, "ERROR: unconnected atoms in file %s\n", file );
	    fprintf ( stderr, err_msg);
	    err_print(err_msg);
	    err_warning2 ( "screen_single_compounds", "skipping ligand");
	  }
	find_carbon_ring_centers ( global->ligand );
	assign_hydrogens ( global->ligand);
	find_cycles ( global->ligand);
	find_flexible_bonds ( global->ligand, global->flex_bond_rules,
			      global->number_of_flex_bond_rules);
	find_hyd_atoms ( global->ligand, 
			 &global->hyd_atom_rules );
	for ( j = 0; j < global->ligand->number_of_atoms; j++ )
	  if ( global->ligand->atoms[j].act != NOTHING )
	    {
	      atom = &global->ligand->atoms[j];
	      global->compound_interactions->act[counter] = atom->act;
	      global->compound_interactions->atom_index[counter] = j;
	      global->compound_interactions->pos[3*counter+X] = atom->pos[X];
	      global->compound_interactions->pos[3*counter+Y] = atom->pos[Y];
	      global->compound_interactions->pos[3*counter+Z] = atom->pos[Z];
	      counter++;
	      compound_counter++;
	    }
	for ( j = 0; j < global->ligand->number_of_carbon_rings; j++ )
	  {
	    global->compound_interactions->act[counter] = HYDROPHOB;
	    global->compound_interactions->atom_index[counter] 
	      = global->ligand->carbon_ring_atom[j];
	    global->compound_interactions->pos[3*counter+X] 
	      = global->ligand->carbon_ring_centers[j][X];
	    global->compound_interactions->pos[3*counter+Y] 
	      = global->ligand->carbon_ring_centers[j][Y];
	    global->compound_interactions->pos[3*counter+Z] 
	      = global->ligand->carbon_ring_centers[j][Z];
	    counter++;
	    compound_counter++;
	  }
	compounds[number_of_compounds].number_of_points = compound_counter;
	strcpy ( compounds[number_of_compounds].name, help );
	number_of_compounds++;
      }
  return number_of_compounds;
}


void  screen_single_compounds ( global_data_pt  global,
				char            **argv,
				int             argc,
				FILE            *fp )
{
  pts_compound_pt compounds;
  char            file[FILENAME_MAX];
  int             array[MAX_NUMBER_OF_MOL2_ATOMS];
  int             result;
  int             i, j;

  err_panic2("screen_single_compounds", 
            "This method of running is obselete and does not work correctly");
  /* we only want to screen single compounds */
  sprintf ( file, "%s/params/hbond.defn", global->hbind_dir );
  read_hyd_defn ( file, &global->hyd_atom_rules);

  if ( fp == NULL )
    /* we haven't opened a file that lists the compounds we have to screen,
       so they are passed as arguments when HBIND was called */
    global->compound_interactions->number_of_compounds 
      = get_compounds_from_args ( global, argv, argc );
  else
    /* the filepointer is not NULL, i.e., we have to read the compounds
       out of the file */
    global->compound_interactions->number_of_compounds 
      = get_compounds_from_file ( global, fp );
  compounds = global->compound_interactions->compounds;
  for ( i = 0; i < global->compound_interactions->number_of_compounds; i++ )
    {
      global->last = TRIANGLE_MATCH;
      global->number_of_screened_compounds++;
      global->compound_index = i;
      global->binding_modes_counter = 0;
      strcpy ( global->ligand_file_name, compounds[i].name );
      printf ( "%s\n", global->ligand_file_name );
      fflush ( stdout );
      for ( j = 0; j < compounds[i].number_of_points; j++ ) array[j] = j;
      for ( j = 2, result = SUCCESS; 
	    j < compounds[i].number_of_points && result == SUCCESS; j++ ){
	  array[2] = array[j];
	  result = compute_triangles ( global, array, 2, 3,
				       compounds[i].first_point );
	}
      if ( global->last == PASSED )
	{
	  printf ( "%s: %d %s\n",
		   global->ligand_file_name,
		   global->binding_modes_counter == 0 ? 
		   1 : global->binding_modes_counter,
		   global->binding_modes_counter == 0 ? 
		   "orientation" : "orientations" );
	  global->number_of_potential_ligands++;
	}
      else
	/* store for the statistics how far this compound made it */
	global->filter_counter[global->last]++;
      fflush ( stdout );
    }
}
