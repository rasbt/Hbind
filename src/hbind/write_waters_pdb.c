#include <stdio.h>
#include <errno.h>
#include <string.h>
#include "err_handle.h"
#include "defs.h"
#include "types.h"

void write_waters_pdb(atom_pt waters, int num_waters, char *filename, 
                      int *water_states, float *water_positions)
{
  atom_pt  water;
  FILE     *fp;
  int      i;
  char     err_msg[2*FILENAME_MAX];
  int errsv;

  if(num_waters <= 0) return;

  fp = fopen ( filename, "w" );
  if ( fp == NULL ) {
    errsv = errno;
    sprintf(err_msg, "Could not open the file: %s\n\t%s", filename, 
            strerror(errsv));
    fprintf(stderr, err_msg);
    err_print(err_msg);
    err_panic2 ( "write_waters_pdb", "file open failed");
  }

  water = waters;
  if(water_states && water_positions){
    for(i = 0; i < num_waters; i++, waters++){
      if(water_states[i] != CONSERVED) continue;

      fprintf(fp, "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  1.00\n",
              water->number, water->number, water_positions[3*i],
              water_positions[3*i + 1], water_positions[3*i + 2]);
    }
  }else{
    for(i = 0; i < num_waters; i++, waters++){
      if(water->state != CONSERVED) continue;

      fprintf(fp, "HETATM %4d  O   HOH %5d     %7.3f %7.3f %7.3f  1.00  1.00\n",
              water->number, water->number, water->pos[0], water->pos[1], 
              water->pos[2]);
    }
  }
  fclose ( fp );
}
