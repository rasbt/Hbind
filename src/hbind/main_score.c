#define _MAIN_
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <popt.h>
#include "defs.h"
#include "types.h"
#include <basics.h>
#include "mymalloc.h"
#include "read_pdb.h"
#include "err_handle.h"
#include "initialize.h"
#include "read_mol2.h"
#include "read_flex_defn.h"
#include "read_hyd_defn.h"
#include "adj_list.h"
#include "read_waters.h"
#include "intra_hbonds.h"
#include "find_hyd_atoms.h"
#include "unbump_side_chains.h"
#include "sum_charges.h"
#include "assign_hydrogens.h"
#include "find_cycles.h"
#include "find_flexible_bonds.h"
#include "score_complex.h"
#include "analyze_ligand.h"
#include "docking_features.h"
#include "match_triangles.h"
#include "distance_matrices.h"
#include <unistd.h>
#include <transform_molecule.h>

#define DISPLAY_HELP_ONLY -1
#define FATAL_ERROR 0
#define BUILD_INTERACTIONS_TABLE 1
#define PRINT_INTERACTIONS 2
#define PRINT_SALTBRIDGES 3
#define PRINT_SUMMARY 4

typedef struct{
  char *prot_fname;
  char *sc_fname;
  char *lig_fname;
  char *waters_fname;
  char *lig_list_fname;
  char *water_list_fname;
  int build_interact_tbl;
  int print_interactions;
  int print_saltbridges;
  int print_summary;
}*cmdline_opts_pt, cmdline_opts_t;

typedef struct node{
  dock_feats_t features;
  char lig_name[FILENAME_MAX];
  struct node *next;
}*features_node_pt, features_node_t;

void build_interact_tbl(features_node_pt features_head,
                        features_node_pt features_last, atom_pt atoms,
                        int num_atoms, residue_pt residues, int num_residues,
                        const char* prot_fname, atom_pt ligand_atoms, dock_feats_pt features, print_saltbridges);

int parse_cmdline(const int argc, const char **argv, cmdline_opts_pt opts);

int load_ligand(char *filename, global_data_pt global);

void set_global_junk(global_data_pt global);

int hphob_sidechain(int res);
int hphob_prot_atom(int res, int type);

int read_moved_target(char *pdb_fname, atom_pt target_atoms);

char cwd[100000];

int main(const int argc, const char **argv)
{
  int rv = 0;
  int USE_TARGET_FILES;
  dock_feats_pt features;
  global_data_pt global;
  cmdline_opts_t cmdline_opts;

  FILE *lig_list_fp;
  char lig_fname[FILENAME_MAX];
  char pdb_fname[FILENAME_MAX];
  int errsv;
  char *pos;
  features_node_pt cnode;
  const features_node_pt features_head =
    (features_node_pt) mymalloc(sizeof(features_node_t));
  features_node_pt features_last = features_head;
  features_head->next = 0;

#ifdef TRACE
  for(rv = 0; rv < argc; ++rv){
    fprintf(stderr, "argv[%d] %s\n", rv, argv[rv]);
  }
  rv = 0;
#endif
  printf ( "\nHBIND Version: %s\n", VERSION);

  rv = parse_cmdline(argc, argv, &cmdline_opts);
  if(rv == FATAL_ERROR) return -1;
  if(rv == DISPLAY_HELP_ONLY) return 0;

  global = initialize_global_data_structure();
  read_pdb(cmdline_opts.prot_fname, global->target_atoms,
           global->target_residues, ALSO_HETERO,
           &global->number_of_target_atoms, &global->number_of_target_residues);
  set_global_junk(global);

  /* Score 1 protein and 1 ligand */
  rv = load_ligand(cmdline_opts.lig_fname, global);
  
  getcwd(cwd, sizeof(cwd));

  if (cmdline_opts.lig_fname[0] == '.') 
      memmove(cmdline_opts.lig_fname, cmdline_opts.lig_fname+1, strlen(cmdline_opts.lig_fname));

  if (cmdline_opts.lig_fname[0] == '/') 
      memmove(cmdline_opts.lig_fname, cmdline_opts.lig_fname+1, strlen(cmdline_opts.lig_fname));

  printf ( "Ligand file: %s/%s", cwd, cmdline_opts.lig_fname);
  
  if(rv != SUCCESS){
    fprintf(stderr, "Reading of %s: failed\n", cmdline_opts.lig_fname);
    exit(-1);
  }

  initialize_inter_dist_matrix(global->target_atom_positions,
                               global->number_of_target_atoms,
                               global->ligand->atom_positions,
                               global->ligand->number_of_atoms,
                               global->target_ligand_sq_dists,
                               DONT_CARE_BUMP_DISTANCE);
  features = &global->current_orientation;
  init_features(features);
  init_features(&global->best_orientation);
  if(cmdline_opts.print_interactions) score_complex(global, features, stdout);
  else score_complex(global, features, 0);


  features_last->next =
    (features_node_pt) mymalloc(sizeof(features_node_t));
  cnode = features_last;
  features_last = features_last->next;
  features_last->next = 0;
  init_features(&cnode->features);
  init_features(&global->best_orientation);
  if(cmdline_opts.print_interactions)
    score_complex(global, &cnode->features, stdout);
  else score_complex(global, &cnode->features, 0);

  if(cmdline_opts.print_summary){
    write_features_line(features, cmdline_opts.lig_fname, stdout, JUST_SCORE);
  }

  /* Build the target matchprints file if desired */
  if(cmdline_opts.build_interact_tbl){
    fprintf(stdout, "\n\n");
    build_interact_tbl(features_head, features_last, global->target_atoms,
                       global->number_of_target_atoms,
                       global->target_residues,
                       global->number_of_target_residues,
                       cmdline_opts.prot_fname, global->ligand->atoms, &cnode->features, 
                       cmdline_opts.print_saltbridges);
  }

/* Score the protein versus each ligand listed in the ligand file and
 * move the side chains if given a sidechain file */

return 0;
}

void build_interact_tbl(features_node_pt features_head,
                        features_node_pt features_last, atom_pt atoms,
                        int num_atoms, residue_pt residues, int num_residues,
                        const char* prot_fname, atom_pt ligand_atoms, 
                        dock_feats_pt features, print_saltbridges)
{
  int num_rows = 0;   /*<! Number of ligands */
  int num_cols = 0;/*<! Number of features (hbonding atoms + hphob sidechains)*/
  int i = 0;
  int j = 0;
  int *cells;     /*<! Think spreadsheet cells */
  int *atom_idz;  /*<! Mapping from atom[i] to a column in cells */
  int *row = 0;
  int col = 0;
  FILE *fout = stdout;  /*<! In the future, we can write to file if desired */
  int *col_sums;
  char **labels;
  char target_name[FILENAME_MAX];
  char *pos;
  int len;


  int lig_act;
  int prot_act;
  int template_ref;
  int lig_idx;
  int ligh_idx;
  int targ_idx;




  fprintf(fout, "++++++++++++++++++++++++++ HBind Interaction Table +++++++++++++++++++++++++\n");
  fprintf(fout, "#            | Ligand Atom -- Protein  Atom | Bond   D-H-A  Ligand-Protein\n");
  fprintf(fout, "#            |  #  type    -- RES   #  type | Dist.  Angle  Interaction\n");
  for(i = 0; i < features->number_of_hbonds; i++ ){
    lig_idx = features->ligand_hbond_idz[i];
    targ_idx = features->target_hbond_idz[i];

    if(targ_idx >= 0){
      if ((features->hbond_angles[i] == 0.0) && features->hbond_dists[i] <= 3.5){
          fprintf(fout, "| metal   %3d %3d  %-5s   -- %-3s %3s %-3s "
           " %6.3f  N/A    ", i+1,
           ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
           atoms[targ_idx].residue,
           atoms[targ_idx].residue_num, atoms[targ_idx].name,
           features->hbond_dists[i]);
           }
      else{
          fprintf(fout, "| hbond   %3d %3d  %-5s   -- %-3s %3s %-3s "
           " %6.3f  %3.1f  ", i+1,
           ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
           atoms[targ_idx].residue,
           atoms[targ_idx].residue_num, atoms[targ_idx].name,
           features->hbond_dists[i], features->hbond_angles[i]);

      }
    }
    lig_act = ligand_atoms[lig_idx].act;
    prot_act = atoms[targ_idx].act;

    if (lig_act == DONOR){
        if (prot_act == ACCEPTOR || prot_act == DONEPTOR)
            fprintf(fout, "Donor - Acceptor\n");
        else if (prot_act == METAL_1 || prot_act == METAL_2)
            fprintf(fout, "Metal interaction\n");
    }
    else if (lig_act == ACCEPTOR){
        if (prot_act == DONOR || prot_act == DONEPTOR)
            fprintf(fout, "Acceptor - Donor\n");
        else if (prot_act == METAL_1 || prot_act == METAL_2)
            fprintf(fout, "Metal interaction\n");
    }
    else if (lig_act == DONEPTOR){
        if (prot_act == ACCEPTOR)
            fprintf(fout, "Donor - Acceptor\n");
        else if (prot_act == DONOR)
            fprintf(fout, "Acceptor - Donor\n");
        else if (prot_act == METAL_1 || prot_act == METAL_2)
            fprintf(fout, "Metal interaction\n");
        else if (prot_act == DONEPTOR)
            fprintf(fout, "Doneptor - Doneptor\n", lig_idx, prot_act);
    }
    else if ((lig_act == METAL_1 || lig_act == METAL_2) && (prot_act == METAL_1 || prot_act == METAL_2))
        fprintf(fout, "Metal interaction\n");
    else
        fprintf(fout, "N/A\n");

  }

  if(print_saltbridges){

    for(i = 0; i < features->number_of_salt_bridges; i++){
      lig_idx = features->ligand_salt_bridge_idz[i];
      targ_idx = features->target_salt_bridge_idz[i];

      lig_act = ligand_atoms[lig_idx].act;
      prot_act = atoms[targ_idx].act;

      if (prot_act == METAL_1 || prot_act == METAL_2 || lig_act == METAL_1 || lig_act == METAL_2) {
          fprintf(fout, "| ld_metal");
        }
      else {
              fprintf(fout, "| saltb   ");
        }


      fprintf(fout, "%3d %3d  %-5s   -- %-3s %3s %-3s  %6.3f  N/A    ", i+1,
      ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
      atoms[targ_idx].residue, atoms[targ_idx].residue_num,
      atoms[targ_idx].name, features->salt_bridge_dists[i]);

      if (lig_act == DONOR){
          if (prot_act == ACCEPTOR || prot_act == DONEPTOR)
              fprintf(fout, "Donor - Acceptor\n");
          else if (prot_act == METAL_1 || prot_act == METAL_2)
              fprintf(fout, "Metal interaction\n");
      }
      else if (lig_act == ACCEPTOR){
          if (prot_act == DONOR || prot_act == DONEPTOR)
              fprintf(fout, "Acceptor - Donor\n");
          else if (prot_act == METAL_1 || prot_act == METAL_2)
              fprintf(fout, "Metal interaction\n");
      }
      else if (lig_act == DONEPTOR){
          if (prot_act == ACCEPTOR)
              fprintf(fout, "Donor - Acceptor\n");
          else if (prot_act == DONOR)
              fprintf(fout, "Acceptor - Donor\n");
          else if (prot_act == METAL_1 || prot_act == METAL_2)
              fprintf(fout, "Metal interaction\n");
          else if (prot_act == DONEPTOR)
              fprintf(fout, "Doneptor - Doneptor\n");
      }
      else if ((lig_act == METAL_1 || lig_act == METAL_2) && (prot_act == METAL_1 || prot_act == METAL_2))
          fprintf(fout, "Metal interaction\n");
      else
          fprintf(fout, "N/A\n");

    }
  }


/* END: NEW STUFF FOR H-BOND Interaction Details, Sebastian Raschka Nov 23, 2016 */

  if(atom_idz) free(atom_idz);
  if(cells) free(cells);
  atom_idz = 0;
  cells = 0;
}

void
print_interactions()



{

}

int parse_cmdline(const int argc, const char **argv, cmdline_opts_pt opts)
{
  int single_docking = 0;
  int file_list = 0;
  char header[2048];
  opts->prot_fname = 0;
  opts->sc_fname = 0;
  opts->lig_fname = 0;
  opts->waters_fname = 0;
  opts->lig_list_fname = 0;
  opts->water_list_fname = 0;
  opts->build_interact_tbl = 1;
  opts->print_interactions = 0;
  opts->print_saltbridges = 0;

  /* Run the "old" school way where protein comes first and then ligand */
  if((argc == 3 || argc == 4) && argv[1][0] != '-'){
    opts->prot_fname = (char *) mymalloc((strlen(argv[1]) + 1) * sizeof(char));
    strcpy(opts->prot_fname, argv[1]);
    opts->lig_fname = (char *) mymalloc((strlen(argv[2]) + 1) * sizeof(char));
    strcpy(opts->lig_fname, argv[2]);
    if(argc == 4){
      opts->waters_fname = (char *) mymalloc((strlen(argv[3]) + 1) *
                                             sizeof(char));
      strcpy(opts->lig_fname, argv[3]);
    }
    return SUCCESS;
  }

  snprintf(header, 2048, "Example: hbind -p <target>.pdb -l <ligand>.mol2\n",
           argv[0], argv[0]);


  struct poptOption mainOptionsTable[] = {
    { "protein", 'p', POPT_ARG_STRING, &opts->prot_fname, 0,
      "Path to protein PDB file", 0},
    { "ligand", 'l', POPT_ARG_STRING, &opts->lig_fname, 0,
      "Path to ligand mol2 file (in docked conformation)", 0},
    { "saltbridges", 's', POPT_ARG_NONE, &opts->print_saltbridges, 0,
      "Include saltbridges in the output", 0},
    { "summary", 's', POPT_ARG_NONE, &opts->print_summary, 0,
      "Include a summary table in the output", 0},
    POPT_AUTOHELP
    POPT_TABLEEND
  };

  poptContext optCon = poptGetContext(argv[0], argc, argv, mainOptionsTable, 0);
  poptSetOtherOptionHelp(optCon, header);

  if(argc < 2) {
    poptPrintUsage(optCon, stderr, 0);
    return DISPLAY_HELP_ONLY;
  }

  int rc;
  /* Process the options */
  for(rc = 0; (rc = poptGetNextOpt(optCon)) >= 0; ){
    switch(rc){
    case BUILD_INTERACTIONS_TABLE:
      opts->build_interact_tbl = BUILD_INTERACTIONS_TABLE;
      break;
    case PRINT_INTERACTIONS:
      opts->print_interactions = PRINT_INTERACTIONS;
      break;
    case PRINT_SALTBRIDGES:
      opts->print_saltbridges = PRINT_SALTBRIDGES;
      break;
    default:
      fprintf(stderr, "Error in processing command line arguments\n");
      return FATAL_ERROR;
    }
  }


  /* Always print the interaction table,
  not the individual interactions due to information overkill,
  Sebastian Raschka 2016 */
  opts->build_interact_tbl = 1;
  opts->print_interactions = 0;


  /* An error occurred during option processing */
  if (rc < -1) {
    fprintf(stderr, "%s: %s\n", poptBadOption(optCon, POPT_BADOPTION_NOALIAS),
            poptStrerror(rc));
    return FATAL_ERROR;
  }

  if(!opts->prot_fname){
    fprintf(stderr, "A protein filename is required\n");
    return FATAL_ERROR;
  }

  if(!opts->lig_fname){
    fprintf(stderr, "A ligand filename is required\n");
    return FATAL_ERROR;
  }

  if(opts->prot_fname && opts->lig_fname) single_docking = 1;
  if(opts->prot_fname && opts->lig_list_fname) file_list = 1;


  poptFreeContext(optCon);
  return SUCCESS;
}

int load_ligand(char *filename, global_data_pt global)
{
  FILE *MOL2;
  int rv = FATAL_FAILURE;

  if((MOL2 = open_mol2(filename)) == NULL) return FATAL_FAILURE;
  if((rv = read_mol2(MOL2, filename, global, NULL)) != SUCCESS){
    fclose(MOL2);
    return rv;
  }
  fclose(MOL2);

  /* -- analyze_ligand() depends on interactions loaded from a pts file --
   * Do it the "hard" way */
  construct_adjacency_list ( global->ligand );
  sum_charges(global->ligand);
  find_hyd_atoms(global->ligand, &global->hyd_atom_rules);
  find_cycles(global->ligand);
  find_flexible_bonds(global->ligand, global->flex_bond_rules,
                        global->number_of_flex_bond_rules);

  if(global->ligand_flag != 0) free(global->ligand_flag);
  global->ligand_flag =
    (short *) mymalloc(3 * global->ligand->number_of_atoms * sizeof(short) );

  return SUCCESS;
}

void set_global_junk(global_data_pt global)
{
  int i;
  char file[FILENAME_MAX + 1];

  strncpy(global->hbind_dir, ".", FILENAME_MAX );
  if(*global->hbind_dir == '\0')
    err_panic("main_score", "HBIND_DIR environment variable not set");

  sprintf(file, "%s/params/flex.defn", global->hbind_dir );
  global->number_of_flex_bond_rules =
    read_flex_defn(file, global->flex_bond_rules);
  sprintf(file, "%s/params/hbond.defn", global->hbind_dir );
  read_hyd_defn(file, &global->hyd_atom_rules);

  /* make backup of target atom positions, since the positions in
   * 'global->target_atoms' are modified when target side chains
   * are rotated during bump resolvement */
  for ( i = 0; i < global->number_of_target_atoms; i++ )
    global->orig_target_atom_act[i] = global->target_atoms[i].act;
  memcpy(global->orig_target_atom_positions, global->target_atom_positions,
         3*global->number_of_target_atoms *
         sizeof(*global->target_atom_positions));
  memset(global->target_rotations, 0, global->number_of_target_residues *
         sizeof(*global->target_rotations));

  /* this is the array for the lookup table of inter-atomic distances,
     before doing the very first bump-check after transforming a ligand,
     this array is filled with the distances between all pairs of ligand
     and target atoms, so that we will avoid most of the calls of
     'dist_fun()' during the modeling of the induced complementarity - Volker*/
  global->target_ligand_distances =
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS *
                        global->number_of_target_atoms * sizeof (float) );
  global->target_ligand_sq_dists =
    (float *) mymalloc (MAX_NUMBER_OF_MOL2_ATOMS *
                        global->number_of_target_atoms * sizeof (float) );
  distance_array(&global->target_dists_array, global->target_atoms,
                 global->number_of_target_atoms, 4.0, 5.0);
  // not used -- need to allocate memory if we decide to use them
  //init_target_nbr_arrays(global->target_atoms, global->number_of_target_atoms);

  /* Allocate & initialize memory for flag arrays used in the scoring function
   */
  global->target_flag =
    (short *) mymalloc(global->number_of_target_atoms * sizeof(short));

#ifdef NON_METALBONDED_REPULSIVE
  /* Build a list of metal atom indices */
  global->number_of_metals = 0;
  global->metal_atom_indices =
    (int*) mymalloc (MAX_TARGET_METAL_ATOMS * sizeof(int));

  for(i = 0; i < global->number_of_target_atoms; i++)
    if(global->target_atoms[i].act == METAL_1 ||
       global->target_atoms[i].act == METAL_2)
      global->metal_atom_indices[global->number_of_metals++] = i;

  if(MAX_TARGET_METAL_ATOMS <= global->number_of_metals)
    err_panic2("main_score", "Number of metals exceed MAX_TARGET_METAL_ATOMS");
#endif

  global->number_of_waters = 0;
}

int hphob_sidechain(int res)
{
  switch(res){
  case ALA:
    return 1;
    break;
  case ARG:
    return 1;
    break;
  case ASP:
    return 1;
    break;
  case ASN:
    return 1;
    break;
  case CYS:
    return 1;
    break;
  case GLN:
    return 1;
    break;
  case GLU:
    return 1;
    break;
  case HIS:
    return 1;
    break;
  case ILE:
    return 1;
    break;
  case LEU:
    return 1;
    break;
  case LYS:
    return 1;
    break;
  case MET:
    return 1;
    break;
  case PCA:
    return 1;
    break;
  case PHE:
    return 1;
    break;
  case PRO:
    return 1;
    break;
  case THR:
    return 1;
    break;
  case TRP:
    return 1;
    break;
  case TYR:
    return 1;
    break;
  case VAL:
    return 1;
    break;
  default:
    return 0;
    break;
  };
}

int hphob_prot_atom(int res, int type)
{
  if(res == ALA && type ==   CB) return 1;
  else if(res == ARG && type ==   CB) return 1;
  else if(res == ARG && type ==   CG) return 1;
  else if(res == ASN && type ==   CB) return 1;
  else if(res == ASP && type ==   CB) return 1;
  else if(res == CYS && type ==   CB) return 1;
  else if(res == CYS && type ==   SG) return 1;
  else if(res == GLN && type ==   CB) return 1;
  else if(res == GLN && type ==   CG) return 1;
  else if(res == GLU && type ==   CB) return 1;
  else if(res == GLU && type ==   CG) return 1;
  else if(res == HIS && type ==   CB) return 1;
  else if(res == ILE && type ==   CB) return 1;
  else if(res == ILE && type ==  CG1) return 1;
  else if(res == ILE && type ==  CG2) return 1;
  else if(res == ILE && type ==  CD1) return 1;
  else if(res == LEU && type ==   CB) return 1;
  else if(res == LEU && type ==   CG) return 1;
  else if(res == LEU && type ==  CD1) return 1;
  else if(res == LEU && type ==  CD2) return 1;
  else if(res == LYS && type ==   CB) return 1;
  else if(res == LYS && type ==   CG) return 1;
  else if(res == LYS && type ==   CD) return 1;
  else if(res == MET && type ==   CB) return 1;
  else if(res == MET && type ==   CG) return 1;
  else if(res == MET && type ==   SD) return 1;
  else if(res == MET && type ==   CE) return 1;
  /* Need to check which atoms are valid in PCA */
  else if(res == PCA && type ==   CB) return 1;
  else if(res == PCA && type ==   CG) return 1;
  else if(res == PHE && type ==   CB) return 1;
  else if(res == PHE && type ==   CG) return 1;
  else if(res == PHE && type ==  CD1) return 1;
  else if(res == PHE && type ==  CD2) return 1;
  else if(res == PHE && type ==  CE1) return 1;
  else if(res == PHE && type ==  CE2) return 1;
  else if(res == PHE && type ==   CZ) return 1;
  else if(res == PRO && type ==   CB) return 1;
  else if(res == PRO && type ==   CG) return 1;
  else if(res == THR && type ==  CG2) return 1;
  else if(res == TRP && type ==   CB) return 1;
  else if(res == TRP && type ==   CG) return 1;
  else if(res == TRP && type ==  CD2) return 1;
  else if(res == TRP && type ==  CE2) return 1;
  else if(res == TRP && type ==  CZ2) return 1;
  else if(res == TRP && type ==  CZ3) return 1;
  else if(res == TRP && type ==  CH2) return 1;
  else if(res == TYR && type ==   CB) return 1;
  else if(res == TYR && type ==   CG) return 1;
  else if(res == TYR && type ==  CD1) return 1;
  else if(res == TYR && type ==  CD2) return 1;
  else if(res == TYR && type ==  CE1) return 1;
  else if(res == TYR && type ==  CE2) return 1;
  else if(res == VAL && type ==   CB) return 1;
  else if(res == VAL && type ==  CG1) return 1;
  else if(res == VAL && type ==  CG2) return 1;
  return 0;
}

/* Assumes that we have a target_pdb as written by write_target_pdb()
 * This means the rotated sidechains are in order and the atoms are in order
 * so that we do not need to do a linear scan per atom or per residue
 */
int
read_moved_target(char *pdb_fname, atom_pt target_atoms)
{
  int i;
  FILE *pdb_file;
  char line[82];
  char pos_str[9];
  /* Note we haven't allocated memory for tmp_atom.pos -- don't use it */
  atom_t tmp_atom;
  atom_pt cur_atom = target_atoms;
  double tmp;

  pdb_file = hbind_fopen(pdb_fname, "r");
  if(pdb_file == NULL){
    printf("Unable to open the moved sidechains file: (%s)\n", pdb_fname);
    return FAILURE;
  }

  while(fgets(line, 82, pdb_file)){
    /* Need to have at least 54 characters or we miss part of the positions */
    if(strlen(line) < 54) continue;
    /* Only handle ATOM lines at the present */
    if(strncmp(line, "ATOM  ", 6) != 0) continue;

    sscanf(line + 12, "%4s", tmp_atom.name);
    /* atom number is bogus -- don't use it */
    strncpy(tmp_atom.residue, line + 17, 3);
    tmp_atom.residue[3] = 0;
    tmp_atom.alt_location = *(line + 16);
    tmp_atom.chain_id = *(line + 21);
    /* insertion code is considered as part of the residue number */
    strncpy(tmp_atom.residue_num, line + 22, 5);
    tmp_atom.residue_num[5] = 0;

    /* Look for the atom in the target that corresponds to the current line */
    for( ; cur_atom; cur_atom++)
      if(tmp_atom.chain_id == cur_atom->chain_id &&
         strcmp(tmp_atom.residue_num, cur_atom->residue_num) == 0 &&
         tmp_atom.alt_location == cur_atom->alt_location &&
         strcmp(tmp_atom.residue, cur_atom->residue) == 0 &&
         strcmp(tmp_atom.name, cur_atom->name) == 0) break;

    /* Hit the end of the target atoms before end of file */
    if(!cur_atom){
      printf("Unable to find the protein atom corresponding to the moved "
             "sidechain atom\n%s\n", line);
      return FAILURE;
    }

    /* Copy the position from line to the target atoms */
    for(i = 0; i < 3; ++i){
      strncpy(pos_str, line + (30 + 8*i), 8);
      pos_str[8] = 0;
      if(!hbind_strtod(pos_str, &tmp)){
        printf("Offending line\n%s\n", line);
        return FAILURE;
      }
      cur_atom->pos[i] = (float) tmp;
    }
  }
  fclose(pdb_file);
  return SUCCESS;

}
