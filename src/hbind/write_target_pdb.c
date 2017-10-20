#include <stdio.h>
#include <string.h>
#include <errno.h>
#include "err_handle.h"
#include "defs.h"
#include "types.h"
#include "dist_fun.h"

void write_target_pdb(residue_pt residues, int num_res, atom_pt atoms, 
                      int num_atoms, char *filename, int *rotations, 
                      float *positions)
{
  FILE     *fp;
  int      count;
  int      i, j;
  char     new_res_num[20];
  char     check_new_res_num[20];
  char     err_msg[2*FILENAME_MAX];
  int errsv;
  atom_pt atom;
  float *atom_pos = 0;

  fp = fopen ( filename, "w" );
  if ( fp == NULL ) {
    errsv = errno;
    sprintf(err_msg, "Could not open the file: %s\n\t%s", filename, 
            strerror(errsv));
    fprintf(stderr, err_msg);
    err_print(err_msg);
    err_panic2("write_target_pdb", "file open failed");
  }

  count = 1;
  for(i = 0; i < num_res; i++){
    if(rotations[i] != YES) continue;

    atom = atoms + residues[i].start_atom;
    if(positions) atom_pos = positions + 3*residues[i].start_atom;
    else atom_pos = atom->pos;
    for(j = 0; j < residues[i].number_of_atoms; j++){
        fprintf(fp, "ATOM  %5d  %-3s%1c%-3s %1c%5s    %7.3f %7.3f %7.3f"
                "  1.00  1.00\n", count, atom->name, atom->alt_location,
                atom->residue, atom->chain_id, atom->residue_num,
                atom_pos[0], atom_pos[1], atom_pos[2]);
      count++;
      atom++;
      atom_pos += 3;
    }
  }

  fclose ( fp );
}
