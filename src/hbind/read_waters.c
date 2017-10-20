#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "read_waters.h" 
#include "err_handle.h"
#include "dist_fun.h" 

int find_hbond_partner(atom_pt target, int num_target_atoms, atom_pt water)
{
  float    min_dist,
           distance;
  int      partner;
  int      i;

  partner = UNKNOWN;
  
  min_dist = 99.9;
  for ( i = 0; i < num_target_atoms; i++ ){
    if( target[i].act == NOTHING ) continue;

    distance = dist_fun ( target[i].pos, water->pos );
    if ( distance < min_dist ){
      partner = i;
      min_dist = distance;
    }
  }

  if (partner == UNKNOWN){
    err_warning2 ( "find_hbond_partner", "unable to find hbond partner");
    return FAILURE;
  }
  else return partner;
}


int read_waters(atom_pt waters, int *num_waters, char *filename, atom_pt target,
                int num_target_atoms)
{
  FILE           *in;
  char           line[MAX_TEMPLATE_LINELENGTH];
  int            i;
  int errsv;
  atom_pt water = waters;
  
  *num_waters = 0;
  in = fopen ( filename, "r" );
  if ( in == NULL ) {
    errsv = errno;
    sprintf(line, "Unable to open the waters file: %s", strerror(errsv));
    err_warning2("read_waters", line);
    return FAILURE;
  }

  while ( fgets ( line, MAX_TEMPLATE_LINELENGTH, in ) != NULL ){ 
    /* comment line, empty line, regular ATOM entry or no water, skip */
    if ( line[0] != 'H' || (strncmp ( line + 17, "HOH", 3 ) != 0 && 
                            strncmp ( line + 17, "WAT", 3 ) != 0 && 
                            strncmp ( line + 17, "H2O", 3 ) != 0 ) ) continue;

    water->number = atoi ( line + 22 );
    water->pos[X] = (float) atof ( line + 30 );
    water->pos[Y] = (float) atof ( line + 38 );
    water->pos[Z] = (float) atof ( line + 46 );
    water->prediction = (float) atof ( line + 54 );
      
    water->type = WATER;
    water->orbit = UNKNOWN;
    water->rad = WATER_RAD;
    water->state = CONSERVED;

    /* find for the closest polar target atom, the water is considered
     * to be connected to this atom, i.e., it will be displaced when
     * this atom moves */
    water->target_atom_index =  
      find_hbond_partner(target, num_target_atoms, water);
    (*num_waters)++;
    water++;
    if(*num_waters >= MAX_BINDING_SITE_WATERS )
      err_panic2("read_waters", "more than MAX_BIDING_SITE_WATERS waters");
  }
  fclose ( in );

  water = waters;
  printf ( "CHECK: water file read\n" );
  for (i = 0; i < 3; i++)
    printf("   Wat %d: %8.3f %8.3f %8.3f %4.2f\n", i, water[i].pos[X],
	   water[i].pos[Y], water[i].pos[Z], water[i].prediction);

  return SUCCESS;
}
       
