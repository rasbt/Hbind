#include <stdio.h> 
#include <stdlib.h>
#include "print_interaction.h"
#include "err_handle.h"
#include "dist_fun.h"

const char INTERACTION_SEPARATOR = '|';

int print_atom(FILE *fout, atom_pt a, int type);

void print_interaction(atom_pt atom1, int atom1_type, int atom1_index, 
                       atom_pt atom2, int atom2_type, int atom2_index,
                       const char *interaction_type, FILE *fout)
{
  fprintf(fout, "%s%c", interaction_type, INTERACTION_SEPARATOR);

  if(atom1_type == atom2_type && atom2_index < atom1_index){
    if(!print_atom(fout, atom2, atom2_type))
      err_panic2("print_interaction", "Second atom has unknown type");
    fprintf(fout, "%c", INTERACTION_SEPARATOR);
    if(!print_atom(fout, atom1, atom1_type))
      err_panic2("print_interaction", "First atom has unknown type");
  }else{
    if(!print_atom(fout, atom1, atom1_type))
      err_panic2("print_interaction", "First atom has unknown type");
    fprintf(fout, "%c", INTERACTION_SEPARATOR);
    if(!print_atom(fout, atom2, atom2_type))
      err_panic2("print_interaction", "Second atom has unknown type");
  }
  fprintf(fout, "%c", INTERACTION_SEPARATOR);
  fprintf(fout, "%f\n", dist_fun(atom1->pos, atom2->pos));
}

int print_atom(FILE *fout, atom_pt a, int type)
{
  int rv = 1;
  char *end;
  if(type == TARGET) fprintf(fout, "T,");
  else if(type == LIGAND) fprintf(fout, "L,");
  else if(type == WATER) fprintf(fout, "W,");
  else rv = 0;

  if(rv){
    fprintf(fout, "%c,%s,%-s,%s,%d,", a->chain_id, a->residue, 
            a->residue_num, a->name, a->atom_number);
    if(a->act == DONOR) fprintf(fout, "D,");
    else if(a->act == ACCEPTOR) fprintf(fout, "A,");
    else if(a->act == DONEPTOR) fprintf(fout, "N,");
    else fprintf(fout, ",");
  }

  return rv;
}
