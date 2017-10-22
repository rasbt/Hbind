#include <stdio.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <string.h>
#include "dist_fun.h"
#include "err_handle.h"
#include "hbond_check.h"
#include "intra_hbonds.h"
#include "mymalloc.h"
#include <intra_hbonds_flag.h>
#include "count_flexible_bonds.h"
#include "check_complementarity.h"
#include "calc_score_from_terms.h"
#include "sf_weights.h"
#include "print_interaction.h"
#include "docking_features.h"

#define between( dist, lower, upper ) ( ((dist) >= (lower - 0.0005)) && ((dist) < (upper + 0.0005)) )

float diff ( float a, float b )
{
  if ( a > b ) return a - b;
  else if ( a < b ) return b - a;
  return 0.0;
}

void print_atom_flag(char *residue_num, char *name, int act, int flag)
{
  printf("%4s %6s ", name, residue_num);
  switch(act){
  case ACCEPTOR:
    printf("  ACCEPTOR ");
    break;
  case DONOR:
    printf("  DONOR    ");
    break;
  case DONEPTOR:
    printf("  DONEPTOR ");
    break;
  case HYDROPHOB:
    printf("  HYDROPHOB");
    break;
  default:
    printf("  *****    ");
    break;
  }

  switch(flag){
  case INITIAL:
    printf("  INITIAL\n");
    break;
  case INTERFACIAL_ATOM:
    printf("  INTERFACIAL_ATOM\n");
    break;
  case METAL_DIRECT_HBOND:
    printf("  METAL_DIRECT_HBOND\n");
    break;
  case SALT_BRIDGE:
    printf("  SALT_BRIDGE\n");
    break;
  case DIRECT_HBOND:
    printf("  DIRECT_HBOND\n");
    break;
  case UNSAT_CHARGE:
    printf("  UNSAT_CHARGE\n");
    break;
  case UNSAT_POLAR:
    printf("  UNSAT_POLAR\n");
    break;
  case INTRA_LIGAND_HBOND:
    printf("  INTRA_LIGAND_HBOND\n");
    break;
  case INTRA_LIGAND_SALT_BRIDGE:
    printf("  INTRA_LIGAND_SALT_BRIDGE\n");
    break;
  case INTRA_TARGET_HBOND:
    printf("  INTRA_TARGET_HBOND\n");
    break;
  case INTRA_TARGET_SALT_BRIDGE:
    printf("  INTRA_TARGET_SALT_BRIDGE\n");
    break;
  default:
    printf("  *****\n");
    break;
  }
}

/*
 *  This function evaluates a ligand in it's neighborhood. Various interactions
 *  are evaluated - most of these are protein-ligand interactions. The scoring function
 *  expression has various terms which are calculated in this function.
 *  Term weights derived from an extensive research study of large data sets of protein-ligand
 *  dockings and complexes. These weights are availables as part of a 2D array in file
 *  inc/sf_weights.h
 */

void score_complex(global_data_pt global, dock_feats_pt features,
                   FILE *good_acts_file)
{
  float dist;
  float sq_dist;
  float hbond_angle;
  float affinity_scores[MAX_NUM_OF_SCORE_FUNCTIONS];
  float sum_pos[3];
  int i, j, k; /* Loop indices */
  char prev_target_atom_residue[4];
  char prev_target_atom_name[5];
  char prev_target_atom_alt_location;

  int number_of_metal_ligand_bonds = 0;
  int number_of_hbonds = 0;
  int number_of_water_mediated_hbonds = 0;
  int number_of_salt_bridges = 0;
  int number_of_unsat_ligand_polar = 0;
  int number_of_unsat_ligand_charge = 0;
  int number_of_unsat_target_polar = 0;
  int number_of_unsat_target_charge = 0;
  int number_of_interfacial_ligand_atoms = 0;
  int number_of_ligand_nonH_atoms = 0;
  int number_of_exposed_hydro_ligand_atoms = 0;
  int number_of_flexible_interfacial_ligand_bonds = 0;
  int number_of_flexible_interfacial_target_bonds = 0;
  int number_of_hphob_hphob_contact = 0;
  int number_of_hphob_hphil_contact = 0;
  int number_of_hphob_hphob_contact_of_one_lig_atm = 0;
  int number_of_target_ligand_phob_phob_interactions = 0;
  int number_of_target_ligand_phil_phil_interactions = 0;
  int number_of_target_ligand_phob_phil_interactions = 0;
  int number_of_target_ligand_interactions = 0;
  int number_of_contacts_of_one_lig_atm = 0;
  int old_hbind_num_of_hphob_hphob_neighbor = 0;
  int found_hbond = 0;
  int found_metal_hbond = 0;
  int found_water_hbond = 0;
  int found_salt_bridge = 0;
  int found_interface = 0;
  int number_of_ligand_neighbors = 0;
  int number_of_ligand_carbons = 0;
  int number_of_exposed_ligand_carbons = 0;
  float total_sum_hphob_hphob = 0.0;
  float total_avg_hphob_hphob = 0.0;
  float total_sum_hphil_hphil = 0.0;
  float total_diff_hphob_hphob = 0.0;
  float total_diff_hphob_hphil = 0.0;
  float total_diff_hphil_hphil = 0.0;
  float target_hydro = 0.0;
  float ligand_hydro = 0.0;
  float sum_hphob_hphob = 0.0;
  float sum_hphil_hphil = 0.0;
  float diff_hphob_hphob = 0.0;
  float diff_hphob_hphil = 0.0;
  float diff_hphil_hphil = 0.0;
  float old_hbind_sum_hphob_hphob = 0.0;
  float old_hbind_total_sum_hphob_hphob = 0.0;
  float avg_hydro = 0.0;
  float diff_hydro = 0.0;
  float diff_interaction_hphob = 0.0;
  float local_hphob_environ = 0.0;
#ifdef SCORING_WATERS
  float **ligand_water_dist = NULL;
  float **target_water_dist = NULL;
  short **water_water_hbond = NULL;       /* records if hbond between waters */
  atom_pt water = global->waters;
#endif
  short *ligand_flag = NULL;     /* status of ligand atom */
  short *target_flag = NULL;     /* status of protein atom */
  short *neighbor_count = NULL;  /* the total number of neighbors of each ligand atom */

  /* get pointers to objects */
  atom_pt target = global->target_atoms;
  atom_pt ligand = global->ligand->atoms;
  float *raw_terms = features->score_component_terms;

  /* assign default score value */
  for(i = 0; i < MAX_NUM_OF_SCORE_FUNCTIONS; i++)
    affinity_scores[i] = global->score_cutoff;

  ligand_flag = global->ligand_flag;
  target_flag = global->target_flag;
  memset(ligand_flag, 0, 3 * global->ligand->number_of_atoms * sizeof(short));
  memset(target_flag, 0, global->number_of_target_atoms * sizeof(short));
  neighbor_count = global->ligand_flag + global->ligand->number_of_atoms;

  /****************************************************************************
   * Go through all ligand and target atoms, and calculate the parameters for
   * the scoring function.  The following distance is defined in file
   * inc/defs.h
   * 	HYDRO_DIST       4.5   
   * 	SALT_BRIDGE      2.5 ~ 4.5
   *	HBOND            2.4 ~ 3.5
   ***************************************************************************/
  for(i = 0; i < global->ligand->number_of_atoms; i++){
    if(ligand[i].type == H) continue;

    /*---------------------------------------------------------
     	We will count the hydrophobic complementarity terms
    	when ligand atom and it's average protein neighbors are:
    	1. hphob - avg hphob
    	2. hphob - avg hphil or hphil - avg hphob
    	3. hphil - avg hphil
  	--------------------------------------------------------*/
    sum_hphob_hphob = 0.0;
    sum_hphil_hphil = 0.0;
    diff_hphob_hphob = 0.0;
    diff_hphob_hphil = 0.0;
    diff_hphil_hphil  = 0.0;
    number_of_hphob_hphob_contact_of_one_lig_atm = 0;
    old_hbind_sum_hphob_hphob = 0.0;
    old_hbind_num_of_hphob_hphob_neighbor = 0;
    local_hphob_environ = 0.0;
    number_of_contacts_of_one_lig_atm = 0;
    diff_interaction_hphob = 0.0;
    prev_target_atom_residue[0] = '\0';
    prev_target_atom_name[0] = '\0';
    prev_target_atom_alt_location = 0;
    number_of_ligand_neighbors = 0;

    /* At first, we try to shift the value about 235,
     * (original range around 0 ~ 635 is changed to -235 ~ 400 )
     * any negative value is considered to be hydrophobic;
     * any positive value is considered to be hydrophilic */
    ligand_hydro = (float) ligand[i].hydro - HYDROPHOB_VALUE_SHIFT;
    number_of_ligand_nonH_atoms ++;
    raw_terms[NUMBER_OF_LIGAND_NONH_ATOMS]++;
    if(ligand[i].type == C) number_of_ligand_carbons++;

    /*-------------------------------------------------------
     * go through all the target atoms
     *---------------------------------------------------------*/
    for(j = 0; j < global->number_of_target_atoms; j++){
      if(target[j].type == H) continue;

      dist = FLT_MAX;
      hbond_angle = FLT_MAX;
      sq_dist =
        global->target_ligand_sq_dists[i*global->number_of_target_atoms + j];

      /* count the protein neighbors of this ligand atom */
      if(ligand_hydro < 0.0 && ligand[i].type == C &&
         MIN_BURY_DIST_2 <= sq_dist && sq_dist <= MAX_BURY_DIST_2)
        neighbor_count[i]++; /*NOTE1*/

      /* count the protein neighbors of this ligand atom
       * for all heavy atoms */
      target_hydro = (float) target[j].hydro - HYDROPHOB_VALUE_SHIFT;

      /* Reset flags that track bond detections.*/
      found_interface = 0;
      found_metal_hbond = 0;
      found_salt_bridge = 0;
      found_hbond = 0;
      found_water_hbond = 0;

      /* search for direct-interfacial atoms and mark them, HYDRO_DIST=4.5
       * changed value to 4.5005*/
      if(sq_dist < HYDRO_DIST_2){
        if(ligand_flag[i] == INITIAL) ligand_flag[i] = INTERFACIAL_ATOM;
        if(target_flag[j] == INITIAL) target_flag[j] = INTERFACIAL_ATOM;

        found_interface = 1;
        old_hbind_sum_hphob_hphob += target[j].hydro;
        old_hbind_num_of_hphob_hphob_neighbor ++;
        number_of_ligand_neighbors++;
      }

      /* we don't want duplicated interaction between the same pair of atoms.
       * Some side chains in protein may have more than one orientations */
      if((strcmp ( prev_target_atom_residue, target[j].residue) == 0
         && strcmp ( prev_target_atom_name, target[j].name) == 0
         && prev_target_atom_alt_location != target[j].alt_location))
        continue;

      strcpy ( prev_target_atom_residue, target[j].residue );
      strcpy ( prev_target_atom_name, target[j].name );
      prev_target_atom_alt_location = target[j].alt_location;

      /*========== CHECK FOR POLAR/CHARGE ATOM INTERACTIONS ==============*/
      if((target[j].act != NOTHING || target[j].charge != 0 ) &&
         (ligand[i].act != NOTHING || ligand[i].charge_sum != 0 ) ){

        /*----------------------------------------------------------------
         * step I: categorise for direct-interactions in following order:
         * 1. metal-hbond
         * 2. salt-bridge
         * 3. hbond
         *----------------------------------------------------------------*/

        /* I.1  mark direct ligand-target metal-hbond
         * metals are treated as donors, defined in read_pdb.c  */

         if(is_hbond(ligand[i], target[j], i, j, sq_dist, &hbond_angle,
                  global)){
           ligand_flag[i] = DIRECT_HBOND;
           target_flag[j] = DIRECT_HBOND;
           found_hbond = 1;
 	        /* NOTE: Hbonds are printed below */
         }

        else if(is_metal_hbond(ligand[i], target[j], sq_dist)){
          ligand_flag[i] = METAL_DIRECT_HBOND;
          target_flag[j] = METAL_DIRECT_HBOND;
          found_metal_hbond = 1;

          if(number_of_hbonds >= MAX_TOTAL_HBONDS)
            err_warning2("check_hydro_compl",
                         "more than MAX_TOTAL_METAL_HBOND metal-hbond");
          else{
            features->ligand_hbond_idz[number_of_hbonds] = i;
            features->target_hbond_idz[number_of_hbonds] = j;
            features->hbond_angles[number_of_hbonds] = 0;
            dist = sqrt(sq_dist);
            features->hbond_dists[number_of_hbonds] = dist;
            number_of_hbonds++;
            raw_terms[NUMBER_OF_HBONDS]++;

	    /*            number_of_metal_ligand_bonds ++;
			  raw_terms[NUMBER_OF_METAL_HBONDS]++;*/
            if(good_acts_file)
              print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                                "Metal", good_acts_file);
#ifdef PRINT_INTERACTIONS
            print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                              "Metal", stdout);
#endif
          }
        }

        /* Count charged ligand-atom metal ion pairs within salt-bridge distance
         * but beyond metal-bond distance as regular salt bridges*/
        else if(is_metal_salt_bridge(ligand[i], target[j], sq_dist)){
          ligand_flag[i] = SALT_BRIDGE;
          target_flag[j] = SALT_BRIDGE;
          found_salt_bridge = 1;

          if ( number_of_salt_bridges >= MAX_TOTAL_SALT_BRIDGES )
            err_warning2("check_hydro_compl",
		        "more than MAX_TOTAL_SALT_BRIDGES salt bridges");
          else{
            features->ligand_salt_bridge_idz[number_of_salt_bridges] = i;
            features->target_salt_bridge_idz[number_of_salt_bridges] = j;
            dist = sqrt(sq_dist);
            features->salt_bridge_dists[number_of_salt_bridges] = dist;

            number_of_salt_bridges++;
            raw_terms[NUMBER_OF_SALT_BRIDGES]++;
            if(good_acts_file)
              print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                                "MetalSB", good_acts_file);
#ifdef PRINT_INTERACTIONS
            print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                              "MetalSB", stdout);
#endif
          }

        /* I.2  mark salt-bridge */
        }else if(is_salt_bridge(target[j].act, target[j].charge, ligand[i].act,
                                ligand[i].charge_sum, sq_dist)){
         /* check if this is a salt bridge, in the CSD and for the
          * target all charged atoms were assigned a charge of -1.0 or +1.0.
          * Charge for protein atoms are assigned in read_pdb.c and charge
          * for ligand atoms are read in from the input mol2 file.
          * so we assume that every charge smaller than -0 or larger
          * than +0 specified an atom that can form a salt bridge */
          ligand_flag[i] = SALT_BRIDGE;
          target_flag[j] = SALT_BRIDGE;
          found_salt_bridge = 1;

          if ( number_of_salt_bridges >= MAX_TOTAL_SALT_BRIDGES )
            err_warning2("check_hydro_compl",
                         "more than MAX_TOTAL_SALT_BRIDGES salt bridges");

          else{
            features->ligand_salt_bridge_idz[number_of_salt_bridges] = i;
            features->target_salt_bridge_idz[number_of_salt_bridges] = j;
            dist = sqrt(sq_dist);
            features->salt_bridge_dists[number_of_salt_bridges] = dist;
            number_of_salt_bridges++;
            raw_terms[NUMBER_OF_SALT_BRIDGES]++;
            if(good_acts_file)
              print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                                "SB", good_acts_file);
#ifdef PRINT_INTERACTIONS
            print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                              "SB", stdout);
#endif
          }
        /* I.3  mark direct ligand-target normal hbond  */
        }
        else if(is_hbond(ligand[i], target[j], i, j, sq_dist, &hbond_angle,
                 global)){
          ligand_flag[i] = DIRECT_HBOND;
          target_flag[j] = DIRECT_HBOND;
          found_hbond = 1;
	        /* NOTE: Hbonds are printed below */
        }


        /***********************************************************************
         * step IV.  count hbonds newly found. Remember: current two atoms might
         * have been marked in the previous loops,  only "found_hbond" can track
         * the new hbond
         **********************************************************************/
        if(found_hbond == 1){
          if( number_of_hbonds >= MAX_TOTAL_HBONDS)
            err_warning2("score_complex",
                         "more than MAX_TOTAL_HBONDS h-bonds");
          else{
            features->ligand_hbond_idz[number_of_hbonds] = i;
            features->target_hbond_idz[number_of_hbonds] = j;
            features->hbond_angles[number_of_hbonds] = hbond_angle;
            dist = sqrt(sq_dist);
            features->hbond_dists[number_of_hbonds] = dist;
            number_of_hbonds++;
            raw_terms[NUMBER_OF_HBONDS]++;
          }
          if(good_acts_file)
            print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                              "HB", good_acts_file);
#ifdef PRINT_INTERACTIONS
          print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i, "HB",
		            stdout);
#endif
        }
      }/*========== END OF CHECK FOR POLAR/CHARGE ATOMS INTERACTIONS =========*/

      if(found_water_hbond == 1){
        err_warning2("score_complex",
                     "Handling of water hbonds is not fully implemented");
        number_of_water_mediated_hbonds++;
      }
      /*================= CHECK FOR HYDROPHOBIC COMPLEMENTARITIES ============*/
      /********************************************************************
       *               step V.   count hydrophobicity values
       *********************************************************************/
      if(found_interface == 1 || found_water_hbond == 1){
        dist = (dist < FLT_MAX ? dist : sqrt(sq_dist));
        local_hphob_environ += target_hydro;
        number_of_contacts_of_one_lig_atm++;
        raw_terms[TOTAL_TARGET_HPHOB_CONTACT] += target_hydro;
        raw_terms[NORM_TARGET_HPHOB_CONTACT] += (target_hydro/(5.0-dist));

        /* both i and j are hphob  */
        if(ligand_hydro < 0.0 && target_hydro < 0.0){
          number_of_target_ligand_phob_phob_interactions++;
          sum_hphob_hphob += target_hydro;
          diff_hphob_hphob += fabs( ligand_hydro - target_hydro );
          raw_terms[TOTAL_TARGET_HPHOB_HPHOB] += target_hydro;
          raw_terms[NORM_TARGET_HPHOB_HPHOB] += (target_hydro/(5.0-dist));
          number_of_hphob_hphob_contact++;
          raw_terms[NUMBER_OF_HPHOB_HPHOB_CONTACT]++;
          number_of_hphob_hphob_contact_of_one_lig_atm++;
          features->target_hphob_contacts[j]++;
          if(good_acts_file)
            print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                              "HPh", good_acts_file);
#ifdef PRINT_INTERACTIONS
          print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                            "HPh", stdout);
#endif

	      /* both i and j are hphil  */
        }else if ( ligand_hydro >= 0.0 && target_hydro >= 0.0 ){
          number_of_target_ligand_phil_phil_interactions++;
          sum_hphil_hphil += target_hydro;
          diff_hphil_hphil += fabs( ligand_hydro - target_hydro );
          /* Either i - hphob, j - hphil  or  i - hphil, j - hphob */
        }else if((ligand_hydro < 0.0 && target_hydro >= 0.0)
	               || (ligand_hydro >= 0.0 &&  target_hydro < 0.0)){
          diff_hphob_hphil += fabs( ligand_hydro - target_hydro );
          number_of_hphob_hphil_contact++;
          raw_terms[NUMBER_OF_HPHOB_HPHIL_CONTACT]++;
          number_of_target_ligand_phob_phil_interactions++;
          if(ligand_hydro < 0.0)
            raw_terms[INCREASE_HPHOB_ENVIRON] += fabs(400.0 - target_hydro);
          if(target_hydro < 0.0)
            raw_terms[INCREASE_HPHOB_ENVIRON] += fabs(400.0 - ligand_hydro);
#ifdef PRINT_INTERACTIONS
          print_interaction(&target[j], TARGET, j, &ligand[i], LIGAND, i,
                            "PhB_PhL", stdout);
#endif
        }
      } /******* end of step V.  *******/
    /*============ END OF CHECK FOR HYDROPHOBIC COMPLEMENTARITIES ============*/
    } /* end of for loop over target atoms */

    /*====== CALCULATE WATER CONTRIBUTION TO FINAL HYDROPHOBICITY VALUE ======*/
    if(ligand_flag[i] != INITIAL){
      number_of_interfacial_ligand_atoms++;
      raw_terms[NUMBER_OF_INTERFACIAL_LIGAND_ATOMS]++;

#ifdef SCORING_WATERS
      /**********************************************************************
       * The new scoring function has been developed, trained and tested in *
       * absence of water molecules. Though, code below implements a preli- *
       * inary water handling, it may make the results worse or better. To  *
       * include water handling,  please type '#define  SCORING_WATERS'     *
       * before the #ifdef statement below                                  *
       **********************************************************************/

      /* count water hydrophobic values */
      if(global->number_of_waters > 0 ){
        /* i - hphob, water - hphil  */
        if(ligand_hydro < 0.0 ){
          for(k = 0; k < global->number_of_waters; k++){
            if(ligand_water_dist[i][k] >= HYDRO_DIST + 0.0005 ||
               water[k].state != CONSERVED ) continue;  /* NOTE1*/

            diff_hphob_hphil += fabs( ligand_hydro - WATER_HYDRO );
            number_of_hphob_hphil_contact ++;
            number_of_target_ligand_phob_phil_interactions++;
            raw_terms[NUMBER_OF_HPHOB_HPHIL_CONTACT]++;
            old_hbind_sum_hphob_hphob += WATER_HYDRO;
            old_hbind_num_of_hphob_hphob_neighbor ++;
#ifdef PRINT_INTERACTIONS
            print_interaction(&water[k], WATER, k, &ligand[i], LIGAND, i,
                              "phB_phL", stdout);
#endif /* endif PRINT_INTERACTIONS */
          }
        /* i - hphil, water - hphil  */
        }else{
         for(k = 0; k < global->number_of_waters; k++){
           /* NOTE1*/
           if(ligand_water_dist[i][k] >= HYDRO_DIST + 0.0005) continue;
           number_of_target_ligand_phil_phil_interactions++;
           sum_hphil_hphil += WATER_HYDRO - HYDROPHOB_VALUE_SHIFT;
           diff_hphil_hphil += fabs( ligand_hydro - WATER_HYDRO);
         }
        }
      }
#endif /* endif SCORING_WATERS section*/

      raw_terms[TOTAL_HPHOB_COMP] += (635.0 - fabs(ligand_hydro - (local_hphob_environ/number_of_contacts_of_one_lig_atm)))/635.0;
      if(ligand_flag[i] != INITIAL)
        raw_terms[TOTAL_LIGAND_HYDRO] += (400.0 - ligand_hydro);
      if(ligand_hydro < 0.0 )
        raw_terms[TOTAL_BURIED_HPHOB] += ligand_hydro;
      if(ligand_hydro < 0.0 && sum_hphob_hphob != 0.0 )
        sum_hphob_hphob += ligand_hydro;
      else if( ligand_hydro >= 0.0 && sum_hphil_hphil != 0.0 )
        sum_hphil_hphil += ligand_hydro;

      total_sum_hphob_hphob += sum_hphob_hphob;
      raw_terms[TOTAL_SUM_HPHOB_HPHOB] += sum_hphob_hphob;
      total_sum_hphil_hphil += sum_hphil_hphil;
      raw_terms[TOTAL_SUM_HPHIL_HPHIL] += sum_hphil_hphil;
      total_diff_hphob_hphob += diff_hphob_hphob;
      raw_terms[TOTAL_DIFF_HPHOB_HPHOB] += diff_hphob_hphob ;
      total_diff_hphob_hphil += diff_hphob_hphil ;
      raw_terms[TOTAL_DIFF_HPHOB_HPHIL] += diff_hphob_hphil ;
      total_diff_hphil_hphil += diff_hphil_hphil;
      raw_terms[TOTAL_DIFF_HPHIL_HPHIL] += diff_hphil_hphil;
      if(sum_hphob_hphob != 0.0 &&
        number_of_hphob_hphob_contact_of_one_lig_atm != 0){
        total_avg_hphob_hphob += sum_hphob_hphob /
          (float) number_of_hphob_hphob_contact_of_one_lig_atm;
        raw_terms[TOTAL_AVG_HPHOB_HPHOB] += sum_hphob_hphob /
          (float) number_of_hphob_hphob_contact_of_one_lig_atm;
        raw_terms[TOTAL_HPHOB_HPHOB_COMP] += (635.0 - fabs(ligand_hydro - ((sum_hphob_hphob - ligand_hydro) / (float) number_of_hphob_hphob_contact_of_one_lig_atm)))/635.0;  /* toneroma 05JUN06 - new hphob term (works better)*/
      }

      /* calculate the old HBIND hphob-hphob term */
      if ( old_hbind_num_of_hphob_hphob_neighbor > 0 ){
        /* hydrophobic values will be negative and in addition to
        * this expanded, so that hydrophobic-hydrophobic complementary
        * gives a better contribution to the overall complementarity */

        if ( ligand_hydro < 0.0 ) ligand_hydro *= (-1.0);
        else ligand_hydro = 0.0;
        avg_hydro = old_hbind_sum_hphob_hphob /
          (float) old_hbind_num_of_hphob_hphob_neighbor - HYDROPHOB_VALUE_SHIFT;
        if ( avg_hydro < 0.0 ) avg_hydro *= (-1.0);
        else avg_hydro = 0.0;
        diff_hydro = fabs ( ligand_hydro - avg_hydro );
        if ( diff_hydro < 32.0 ) diff_hydro = 32.0;

        /* take the average value of the hydrophilicities for
        * computing the overall contribution of this ligand atom */
        old_hbind_total_sum_hphob_hphob +=
          ( ligand_hydro + avg_hydro ) / 2.0 / diff_hydro;
        raw_terms[OLD_HBIND_TOTAL_SUM_HPHOB_HPHOB] +=
          ( ligand_hydro + avg_hydro ) / 2.0 / diff_hydro;
      }
    }

    if( ligand[i].type == C && number_of_ligand_neighbors == 0)
      number_of_exposed_ligand_carbons++; /* This ligand carbon is exposed.*/
  } /* end of for loop of ligand atoms */

  /*****************************************************************************
   * count the number of interfacial rotatable bonds in ligand & protein
   * (not including the terminal bonds)
   ****************************************************************************/
  number_of_flexible_interfacial_ligand_bonds =
    count_flexible_ligand_bonds ( global->ligand, ligand_flag );
  raw_terms[NUMBER_OF_FLEXIBLE_INTERFACIAL_LIGAND_BONDS ] =
    number_of_flexible_interfacial_ligand_bonds;
  number_of_flexible_interfacial_target_bonds =
    count_flexible_target_bonds ( global, target_flag );

  /*****************************************************************************
   * Mark interactions for intra-ligand, ligand-water, intra-target and
   * target-water Only those interfacial atoms will be checked.
   ****************************************************************************/
  for ( i = 0; i < global->ligand->number_of_atoms; i++ ){
    if( ligand_flag[i] == INTERFACIAL_ATOM && ligand[i].act != NOTHING
	     && ligand[i].type != H){
      /* mark intra-ligand salt-bridge */
      for( k = 0; k < global->ligand->number_of_atoms; k++ )
        if( i != k && ligand[k].act != NOTHING && ligand[i].type != H ){
          sq_dist = squared_dist(ligand[i].pos, ligand[k].pos);
          if(is_salt_bridge(ligand[i].act, ligand[i].charge_sum, ligand[k].act,
                            ligand[k].charge_sum, sq_dist )){
            ligand_flag[i] = INTRA_LIGAND_SALT_BRIDGE;
            if(ligand_flag[k] == INTERFACIAL_ATOM)
              ligand_flag[k] = INTRA_LIGAND_SALT_BRIDGE;
            break;
          }
        }
      /* if non-salt bridge, mark ligand-water hbond or intra-ligand hbond */
      if ( ligand_flag[i] == INTERFACIAL_ATOM ){
        intra_ligand_hbonds_flag ( global, i, ligand_flag );
	/* mark ligand-water hbond */
        if ( ligand_flag[i] == INTERFACIAL_ATOM )
	  ligand_to_water_hbonds_flag ( global, i, ligand_flag );
      }
    }
  }

  intra_target_polar_flag(global->target_residues, global->target_atoms,
                          global->number_of_target_atoms, target_flag, stdout,
                          &raw_terms[NUMBER_OF_INTRA_TARGET_SALT_BRIDGES],
                          global);

  /**************************************************************************
   * go through all the ligand atoms again to count unsatisfied and repulsive
   ***************************************************************************/
  for(i = 0; i < global->ligand->number_of_atoms; i++ ){
    /*SA+MIZ 2/17/05*/
    if(ligand_flag[i] != INTERFACIAL_ATOM || ligand[i].act == NOTHING
       || ligand[i].type == H) continue;

    /****************************************************************
     * So far we know the ligand atom i does not have any partner (either
     * interation nor repulsive partner). So it's marked as unsatisfied
     *
     * Check for unsatisfied charge (no repulsive) for ligand atom i
     ******************************************************************/
    if( ligand[i].charge_sum != 0 ){
      ligand_flag[i] = UNSAT_CHARGE;
      number_of_unsat_ligand_charge++;
      raw_terms[NUMBER_OF_UNSAT_LIGAND_CHARGE]++;
#ifdef PRINT_INTERACTIONS
      printf("L Unsat Charge: %9.4f %9.4f %9.4f %d\n", ligand[i].pos[0],
             ligand[i].pos[1], ligand[i].pos[2], ligand[i].act);
#endif
    /********************************************************
     * Check for unsatisfied polar (A/D/N) for ligand atom i
     ************************************************************/
    }else if ( ligand[i].act != NOTHING ){
      ligand_flag[i] = UNSAT_POLAR;
      number_of_unsat_ligand_polar ++;
      raw_terms[NUMBER_OF_UNSAT_LIGAND_POLAR]++;

#ifdef PRINT_INTERACTIONS
      printf("L Unsat Polar: %9.4f %9.4f %9.4f %d\n", ligand[i].pos[0], ligand[i].pos[1],
             ligand[i].pos[2], ligand[i].act); /* 08MAY06 toneroma*/
#endif
    }
  }

#ifdef TRACE
  for(i = 0; i < global->ligand->number_of_atoms; i++ )
    if(ligand_flag[i] != INITIAL)
      print_atom_flag("", ligand[i].name, ligand[i].act, ligand_flag[i]);
#endif

  /*-------------------------- Iterate over TARGET atoms ---------------------*/
  for(j = 0; j < global->number_of_target_atoms; j++){
    if(target_flag[j] != INITIAL)
      raw_terms[TOTAL_TARGET_HYDRO] += (635 - target[j].hydro);
    if(target_flag[j] != INTERFACIAL_ATOM || target[j].act == NOTHING) continue;

    /****************************************************************
     * Check for unsatisfied charge (no repulsive) for target atom J
     *****************************************************************/
    if(target[j].charge != 0){
      target_flag[j] = UNSAT_CHARGE;
      number_of_unsat_target_charge++;
      raw_terms[NUMBER_OF_UNSAT_TARGET_CHARGE]++;
#ifdef PRINT_INTERACTIONS
      /* 08MAY06 toneroma*/
      printf("T Unsat Charge: %7.3f %7.3f %7.3f %d\n", target[j].pos[0],
             target[j].pos[1], target[j].pos[2], target[j].act);
#endif

    /********************************************************
     * Check for unsatisfied polar (A/D/N) for target atom j
     **********************************************************/
    }else if ( target[j].act != NOTHING ){
      target_flag[j] = UNSAT_POLAR;
      number_of_unsat_target_polar++;
      raw_terms[NUMBER_OF_UNSAT_TARGET_POLAR]++;
#ifdef PRINT_INTERACTIONS

      printf("T: %7.3f %7.3f %7.3f %d\n", target[j].pos[0],
             target[j].pos[1], target[j].pos[2], target[j].act);
#endif
    }
  }  

#ifdef TRACE
  for(j = 0; j < global->number_of_target_atoms; j++)
    if(target_flag[j] != INITIAL)
      print_atom_flag(target[j].residue_num, target[j].name, target[j].act,
                      target_flag[j]);
#endif

  /*****************************************************************************
   * go through all the ligand atoms again to count hydrophobic atoms with no
   * neighbor
   ****************************************************************************/
  for( i = 0; i < global->ligand->number_of_atoms; i++){
    /* check hydrophobic carbons without protein neighbor, neighbors are only
     * counted for hydrophobic carbons in ligand in the previous loops */
    /* we need to check if this ligand atoms is buried in ligand.  If the
     * total number of ligand neighbors is over MIN_NUMBER_BURY_NEIGHBOR, it
     * is a buried atom; otherwise, it is an exposed one and waters can come
     * close to it */
    if(ligand[i].type != C || ligand[i].hydro >= HYDROPHOB_VALUE_SHIFT ||
       neighbor_count[i] > 0) continue;

    memset(sum_pos, 0, 3*sizeof(*sum_pos));
    for(j = 0; j < global->ligand->number_of_atoms; j++ ){
      if( ligand[i].type == H ) continue;

      sq_dist = squared_dist ( ligand[i].pos, ligand[j].pos );
      /* if  ( dist > MIN_BURY_DIST && dist < MAX_BURY_DIST ) */
      if(MIN_BURY_DIST_2 <= sq_dist && sq_dist <= MAX_BURY_DIST_2){
        neighbor_count[i]++;
        for(k = 0; k < 3; ++k) sum_pos[k] += ligand[j].pos[k];
      }
    }

    /*printf("carbon %d: %s has %d neighbor\n", i, ligand[i].name, neighbor_count[i] );*/
    for(k = 0; k < 3; ++k) sum_pos[k] /= (float) neighbor_count[i];
    sq_dist = squared_dist ( ligand[i].pos, sum_pos );
    /*printf("carbon %d: %s has %d neighbor, dist to center is %f\n", i, ligand[i].name, neighbor_count[i], dist );*/
    /* NOTE1*/
    if((neighbor_count[i] < MIN_NUMBER_BURY_NEIGHBOR || sq_dist >= 1.0 )
       &&( ligand[i].type != H )){
      number_of_exposed_hydro_ligand_atoms ++;
      raw_terms[NUMBER_OF_EXPOSED_HYDRO_LIGAND_ATOMS]++;
#ifdef PRINT_INTERACTIONS
      printf(" Neighbors : %d, center %6.3f", neighbor_count[i], sqrt(sq_dist));
      printf("L: %9.4f %9.4f %9.4f %d\n", ligand[i].pos[0], ligand[i].pos[1], ligand[i].pos[2], ligand[i].act); 
#endif

    }
  }

  /****************************************************************************
   * Compute terms and "scoring function" values
   ****************************************************************************/
  raw_terms[TOTAL_SUM_HPHOB_HPHOB] /= 1000;
  raw_terms[TOTAL_DIFF_HPHOB_HPHOB] /= 1000;
  raw_terms[TOTAL_DIFF_HPHOB_HPHIL] /= 1000;
  raw_terms[TOTAL_SUM_HPHIL_HPHIL] /= 1000;
  raw_terms[TOTAL_DIFF_HPHIL_HPHIL] /= 1000;
  raw_terms[TOTAL_TARGET_HPHOB_HPHOB] /= 1000;
  raw_terms[NORM_TARGET_HPHOB_HPHOB] /= 1000;
  raw_terms[TOTAL_TARGET_HPHOB_CONTACT] /= 1000;
  raw_terms[NORM_TARGET_HPHOB_CONTACT] /= 1000;
  raw_terms[TOTAL_TARGET_HYDRO] /= 1000;
  raw_terms[TOTAL_LIGAND_HYDRO] /= 1000;
  raw_terms[INCREASE_HPHOB_ENVIRON] /= 1000;
  raw_terms[NUMBER_OF_UNSAT_POLAR] =
    raw_terms[NUMBER_OF_UNSAT_TARGET_POLAR] +
    raw_terms[NUMBER_OF_UNSAT_LIGAND_POLAR];
  raw_terms[NUMBER_OF_UNSAT_CHARGE] =
    raw_terms[NUMBER_OF_UNSAT_TARGET_CHARGE] +
    raw_terms[NUMBER_OF_UNSAT_LIGAND_CHARGE] ;
  raw_terms[RATIO_OF_INTERFACIAL_LIGAND_HEAVY_ATOMS] =
    (float) number_of_interfacial_ligand_atoms /
    (float) number_of_ligand_nonH_atoms;
  raw_terms[NUMBER_OF_LIGAND_FLEXIBLE_BONDS] =
    global->ligand->number_of_flexible_bonds;
  raw_terms[NUMBER_OF_ALL_INTERFACIAL_FLEXIBLE_BONDS] =
    number_of_flexible_interfacial_ligand_bonds +
    number_of_flexible_interfacial_target_bonds;
  raw_terms[TOTAL_AVG_HPHOB_HPHOB] /= 1000;
  raw_terms[NORM_HPHOB_COMP] = (float) number_of_target_ligand_interactions /
    (float) number_of_ligand_nonH_atoms;
  raw_terms[NORM_HPHOB_COMP_2] = (float) number_of_target_ligand_interactions /
    (float) number_of_interfacial_ligand_atoms;
  raw_terms[NORM_HPHOB_COMP_3] = (float) number_of_target_ligand_interactions /
   (float) number_of_ligand_nonH_atoms /
   (float) number_of_interfacial_ligand_atoms;
  raw_terms[CONTACT_HPHOB_HPHOB] =
    number_of_target_ligand_phob_phob_interactions;
  raw_terms[CONTACT_HPHIL_HPHIL] =
    number_of_target_ligand_phil_phil_interactions;
  raw_terms[CONTACT_HPHOB_HPHIL] =
    number_of_target_ligand_phob_phil_interactions;
  raw_terms[HBOND_RATIO] = (float)( number_of_hbonds + number_of_hbonds + number_of_salt_bridges + number_of_salt_bridges + number_of_metal_ligand_bonds + number_of_metal_ligand_bonds)/((float)( number_of_hbonds + number_of_hbonds + number_of_salt_bridges + number_of_salt_bridges + number_of_metal_ligand_bonds + number_of_metal_ligand_bonds + number_of_unsat_ligand_polar + number_of_unsat_ligand_charge + number_of_unsat_target_polar + number_of_unsat_target_charge) + 0.00001);
  raw_terms[NORM_POLAR] = (number_of_hbonds + number_of_salt_bridges +
                           number_of_metal_ligand_bonds) *
                          number_of_ligand_nonH_atoms;
  raw_terms[NORM_UNSAT] = (number_of_unsat_ligand_polar +
                           number_of_unsat_ligand_charge +
                           number_of_unsat_target_polar +
                           number_of_unsat_target_charge) *
                          number_of_ligand_nonH_atoms;
  raw_terms[INTERMEDIATE_OVERLAP] = global->intermediate_overlap;
  raw_terms[TOTAL_OVERLAP] = features->total_overlap;
  score_from_terms(raw_terms, affinity_scores, MAX_NUM_OF_SCORE_FUNCTIONS);

#define AFFINITY_SCORE 62
#define ORIENTATION_SCORE 63
  features->orient_score = affinity_scores[63]; /* scoring func 26*/
  features->affi_score = affinity_scores[62]; /* scoring func 62*/

/*------- OUTPUT RAW TERMS -------------*/
#ifdef DISPLAY_RAW_TERMS

  for(i = 0; i < 35; i++ ) printf("%5.3f\t", raw_terms[i]);
  printf(" \n");
#endif


  features->number_of_hbonds = number_of_hbonds ;
  features->number_of_salt_bridges = number_of_salt_bridges;
  features->number_of_metal_ligand_bonds = number_of_metal_ligand_bonds;
  features->number_of_interfacial_unsatisfied_polar_atoms =
    raw_terms[NUMBER_OF_UNSAT_POLAR];
  features->number_of_interfacial_unsatisfied_charged_atoms =
    raw_terms[NUMBER_OF_UNSAT_CHARGE];
  features->buried_carbons = 1.0 - (float) number_of_exposed_ligand_carbons /
    (float) number_of_ligand_carbons;
}
