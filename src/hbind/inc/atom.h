#ifndef ATOM_HEADER_FILE_INCLUDED
#define ATOM_HEADER_FILE_INCLUDED

#include <defs.h>

/*! data structure describing a bit string */
typedef struct {
  unsigned char *bits;            /*!< bits, stored in a byte vector */
  unsigned int len;               /*!< number of bits in this bit string */
}*bitstring_pt, bitstring_t;

/*! atom entry in mol2 or pdb file */
typedef struct atom_struct{
  /* components for all atom types */
  char          name[5];                            /* atom name in the file */
  float *pos;
  float         rad;                                 /* van der Waals radius */
  float         charge;                      /* partial charge for this atom */
  float         charge_sum;
  int           atom_number; /* the number assigned to this atom in the file */
  int           type;                                    /* type of the atom */
  int           act;                   /* DONOR, ACCEPTOR, DONEPTOR, NOTHING */
#if 0
  int           *neighbors;  /* indices of neighbors for intramolecular bump */
#endif
  struct atom_struct **neighbors;
  float         *neighbor_dist;
  size_t        num_nbrs;
  size_t        num_added_nbrs;
  int           hydro;                               /* hydrophilicity value */

  /* components only relevant for mol2 atoms */
  char          type_str[7];                /* string defining the atom type */
  int           orbit;            /* orbit or whatever follows the . in type */
  int           number;         /*!< ligand ==> # assigned for fragmentation */
  bitstring_pt  fragments;                  /* fragments the atom is part of */
  int           fragment;  /* number of the fragment this atom is located in */
  int           subst_id;                                 /* substructure id */
  char          subst_name[MAX_SUBST_NAME_LENGTH];      /* substructure name */

  /* components only relevant for PDB atoms */
  int           level;                        /* level of atom in side-chain */
  int           residue_index; /* index of the residue in the residue vector */
  int           res;                          /* the residue type identifier */
  char          residue[4];                                  /* residue type */
  char          alt_location;               /* alternate location identifier */
  char          chain_id;                                /* chain identifier */
  char          insertion_code;                            /* insertion code */
  char          residue_num[6];                            /* residue number */

  /* components only relevant for waters */
  int                target_atom_index;  /* atom it is connected to (H-bond) */
  int                state;                        /* CONSERVED or DISPLACED */
  float              prediction;           /* the likelyhood of conservation */
}*atom_pt, atom_t;


#endif
