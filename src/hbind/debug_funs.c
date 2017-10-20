#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "types.h"
#include "err_handle.h"

const char*
act2str(int act)
{

  if(act == NOTHING) return "NOTHING";
  else if(act == ACCEPTOR) return "ACCEPTOR";
  else if(act == DONOR) return "DONOR";
  else if(act == HYDROPHOB) return "HYDROPHOB";
  else if(act == DONEPTOR) return "DONEPTOR";
  else if(act == DONOR_HYDROGEN) return "DONOR_HYDROGEN";
  else if(act == METAL_1) return "METAL_1";
  else if(act == METAL_2) return "METAL_2";

#if 0
/* summed interaction types*/
ACCEPTOR_ACCEPTOR
DONOR_DONOR
HYDROPHOB_HYDROPHOB
ACCEPTOR_DONOR
ACCEPTOR_DONEPTOR
DONOR_DONEPTOR
DONEPTOR_DONEPTOR
ACCEPTOR_DONOR_HYDROGEN
DONEPTOR_DONOR_HYDROGEN
ACCEPTOR_METAL_1
ACCEPTOR_METAL_2
DONEPTOR_METAL_1
DONEPTOR_METAL_2
#endif

}

void 
write_mol2(char *filename, float *positions, global_data_pt global)
{
  FILE *fp;
  int i;
  char err_msg[2*FILENAME_MAX];
  int errsv;
  float *position;

  molecule_pt ligand = global->ligand;
  atom_pt atoms = global->ligand->atoms;
  bond_pt bonds = global->ligand->bonds;
  /* do not print the added hydrogens that have been necessary to
     check the donor/acceptor and flexible bond rules */
  int num_atoms = ligand->number_of_atoms - ligand->number_of_added_hydrogens;
  /* do not print bonds to added hydrogens */
  int num_bonds = ligand->number_of_bonds - ligand->number_of_added_hydrogens;

  fp = fopen ( filename, "w" );
  if ( fp == NULL ) {
    errsv = errno;
    sprintf(err_msg, "Could not open the file: %s\n\t%s", filename,
            strerror(errsv));
    fprintf(stderr, err_msg);
    err_print(err_msg);
    err_panic2("write_mol2", "file open failed");
  }

  fprintf(fp, "@<TRIPOS>MOLECULE\n" );
  fprintf(fp, "%s\n", ligand->name );
  fprintf(fp, "%5d %5d %5d     0     0\n", num_atoms, num_bonds,
          ligand->number_of_substructures );
  fprintf ( fp, "SMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n" );
  for(i = 0; i < num_atoms; i++ ){
    if(positions) position = positions + 3*i;
    else position = atoms[i].pos;
    fprintf(fp, "  %5d %-8s %9.4f %9.4f %9.4f %-9s%2d %s       %7.4f \n",
            atoms[i].atom_number, atoms[i].name, position[0], position[1],
            position[2], atoms[i].type_str, atoms[i].subst_id,
            atoms[i].subst_name, atoms[i].charge );
  }

  fprintf ( fp, "@<TRIPOS>BOND\n" );
  for(i = 0; i < num_bonds; i++)
    fprintf(fp, "%5d %5d %5d %s\n", bonds[i].number,
            atoms[bonds[i].atom1].atom_number,
            atoms[bonds[i].atom2].atom_number, bonds[i].type_str);

  fprintf ( fp, "@<TRIPOS>SUBSTRUCTURE\n" );
  for(i = 0; i < ligand->number_of_substructures; i++)
    fprintf( fp, "%s", ligand->substructure[i]);
  fclose ( fp );
}

void
write_pdb(atom_pt atoms, int num_atoms, char *filename, float *positions)
{
  FILE     *fp;
  int       j;
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
    err_panic2("write_pdb", "file open failed");
  }

  atom = atoms;
  atom_pos = positions;
  for(j = 0; j < num_atoms; j++){
    fprintf(fp, "ATOM  %5d  %-3s%1c%-3s %1c%5s    %7.3f %7.3f %7.3f"
            "  1.00  1.00\n", atom->atom_number, atom->name, 
            atom->alt_location, atom->residue, atom->chain_id, 
            atom->residue_num, atom_pos[0], atom_pos[1], atom_pos[2]);
    atom++;
    atom_pos += 3;
  }

  fclose ( fp );
}
