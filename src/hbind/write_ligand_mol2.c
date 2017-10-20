#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "docking_features.h"
#include "err_handle.h"
#include "calc_score_from_terms.h"
#include "sf_weights.h"

#define AFFINITY_SCORE 62  
#define ORIENTATION_SCORE 63

void write_ligand_mol2(char *filename, float *positions, 
                       dock_feats_pt features, global_data_pt global)
{
  FILE *fp;
  int i;
  char err_msg[2*FILENAME_MAX];
  int errsv;
  float *position;

  molecule_pt ligand = global->ligand;
  atom_pt atoms = global->ligand->atoms;
  bond_pt bonds = global->ligand->bonds;
  atom_pt target_atoms = global->target_atoms;

  int num_atoms = ligand->number_of_atoms - ligand->number_of_added_hydrogens;

  int num_bonds = ligand->number_of_bonds - ligand->number_of_added_hydrogens;

  fp = fopen ( filename, "w" );
  if ( fp == NULL ) {
    errsv = errno;
    sprintf(err_msg, "Could not open the file: %s\n\t%s", filename,
            strerror(errsv));
    fprintf(stderr, err_msg);
    err_print(err_msg);
    err_panic2 ( "write_ligand_mol2", "file open failed");
  }
  

  features->binding_modes_counter = global->binding_modes_counter;
  
  write_features_header(features, fp, global->ligand_file_name, 
                        global->template_interactions, atoms, target_atoms);
 
  fprintf(fp, "@<TRIPOS>MOLECULE\n" );

  fprintf(fp, "%s%s\n", features->ligand_name_noconf, features->ligand_conf );
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
