#include <string.h>
#include "types.h"
#include "hbind_itimer.h"
#include "sf_weights.h"
#include "calc_score_from_terms.h"
#include "mymalloc.h"

#define AFFINITY_SCORE 62  		/* affinity trained funtion # 66 */

void compute_scores(dock_feats_pt features);

void write_features_line(dock_feats_pt features, char *lig_fname, FILE *fout,
                         int docking, global_data_pt global)
{



  compute_scores(features);

  fprintf(fout,
          "+++++++++++++++++ Summary +++++++++++++++++++\n"
          /* "| Buried Protein Hydrophobic Term     : %5.3f\n" */
          /* "| Hydrophobic Complementarity Term    : %5.3f\n" */
          /* "| Polar Component Term                : %5.3f\n" */
          "| Protein-Ligand Hydrophobic Contacts : %5.0f\n"
          "| Protein-Ligand H-bonds              : %5d\n"
          "| Protein-Ligand Salt-bridges         : %5d\n"
          "| Metal-Ligand Bonds                  : %5d\n"
          /*"| Interfacial Unsatisfied Polar Atoms : %5d\n"*/
          /*"| Interfacial Unsatisfied Polar Atoms : %5d\n"*/
          /*"| Buried Carbons (%%)                  : %5.0f\n"*/
          "| SLIDE OrientScore                   : %5.3f\n"
          "| SLIDE AffiScore (heavy atoms)       : %5.3f\n"
          "| SLIDE AffiScore                     : %5.3f",
          /*features->hphob_score, */
          /*features->polar_score,*/
          /*features->unsat_polar_score, */
          features->contact_hphob_hphob,
          features->number_of_hbonds, 
          features->number_of_salt_bridges,
          features->number_of_metal_ligand_bonds,
          /*features->number_of_interfacial_unsatisfied_polar_atoms,*/
          /*features->number_of_interfacial_unsatisfied_charged_atoms,*/
          /*features->buried_carbons*100,*/
          features->orient_score,
          features->affiefficient_score,
          features->affi_score);
  if(docking)
    fprintf(fout, " %3d %5.3f\n", features->number_of_bumps,
            features->total_overlap);
  else fprintf(fout, "");
}

void init_features(dock_feats_pt features)
{
  features->ligand_name[0] ='\0';
  features->ligand_name_noconf[0] = '\0';
  features->ligand_conf[0] = '\0';
  features->binding_modes_counter = 0;
  memset(features->ligand_hbond_idz, 0,
         MAX_TOTAL_HBONDS*sizeof(*features->ligand_hbond_idz));
  memset(features->target_hbond_idz, 0,
         MAX_TOTAL_HBONDS*sizeof(*features->target_hbond_idz));
  memset(features->hbond_angles, 0,
         MAX_TOTAL_HBONDS*sizeof(*features->hbond_angles));
  memset(features->hbond_dists, 0,
         MAX_TOTAL_HBONDS*sizeof(*features->hbond_dists));
  memset(features->ligand_salt_bridge_idz, 0,
         MAX_TOTAL_SALT_BRIDGES*sizeof(*features->ligand_salt_bridge_idz));
  memset(features->target_salt_bridge_idz, 0,
         MAX_TOTAL_SALT_BRIDGES*sizeof(*features->target_salt_bridge_idz));
  memset(features->salt_bridge_dists, 0,
         MAX_TOTAL_SALT_BRIDGES*sizeof(*features->salt_bridge_dists));
  memset(features->score_component_terms, 0,
         NO_OF_SCORE_COMPONENT_TERMS*sizeof(*features->score_component_terms));
  memset(features->target_hphob_contacts, 0,
         MAX_PDB_ATOMS*sizeof(*features->target_hphob_contacts));
  features->flag = 0; /* indicates passed condition */
  features->orient_score = 1E20;
  features->affiefficient_score = 1E20;
  features->affi_score = 1E20;
  features->hphob_score = 0.0;
  features->polar_score = 0.0;
  features->unsat_polar_score = 0.0;
  features->contact_hphob_hphob = 0;
  features->number_of_hbonds = 0;
  features->number_of_salt_bridges = 0;
  features->number_of_metal_ligand_bonds = 0;
  features->number_of_interfacial_unsatisfied_polar_atoms = 0;
  features->number_of_interfacial_unsatisfied_charged_atoms = 0;
  features->buried_carbons = 0.0;
  features->number_of_bumps = 0;
  features->total_overlap = 0.0;
  features->num_anchor_corrections = 0;
  features->num_mean_field_optimizations = 0;
  features->num_side_chain_rotations = 0;
  features->num_ligand_side_chain_rotations = 0;
  features->num_target_side_chain_rotations = 0;
  memset(features->template_matches, 0, 3*sizeof(*features->template_matches));
  memset(features->ligand_matches, 0, 3*sizeof(*features->ligand_matches));
  memset(features->matched_type, 0,
         MAX_INTER_PAIRS*sizeof(*features->matched_type));
}

void write_features_header(dock_feats_pt features, FILE *fp, char *lig_fname,
                           interaction_pt template_interactions,
                           atom_pt ligand_atoms, atom_pt target_atoms)
{
  char time_stamp[80];
  int template_act;
  int template_ref;
  int lig_idx;
  int targ_idx;
  int i;

  compute_scores(features);

  /* May need to be updated to reflect which ligands are written out */
  fprintf(fp, "# Ligand [name]_[conformer]_[binding mode]: %s%s_%04d\n",
          features->ligand_name_noconf, features->ligand_conf, features->binding_modes_counter);
  fprintf(fp, "# HBIND Version: %s\n", VERSION);
  fprintf(fp, "# Time stamp: %s\n", hbind_get_local_time(time_stamp, 80));
  fprintf(fp, "# Affiscore (kcal/mol), for ranking ligand candidates: %5.3f\n",
          features->affi_score);
  fprintf(fp, "# Ligand Efficiency ((Affinity Score - Constant Term) / "
          "# non-H atoms): %5.3f\n", features->affiefficient_score);

  fprintf(fp, "# =============================================================="
          "=============\n");
  fprintf(fp, "# Weighted term values\n");
  fprintf(fp, "# hydrophobic complementarity term (kcal/mol)       : %5.3f \n",
          features->hphob_score);
  fprintf(fp, "# polar term (kcal/mol)                             : %5.3f \n",
          features->polar_score);
  fprintf(fp, "# unsatisfied polar term (kcal/mol)                 : %5.3f \n",
          features->unsat_polar_score);
  fprintf(fp, "# constant term (kcal/mol)                          : %5.3f \n",
          affinity_wts[AFFINITY_SCORE][TERM_0]);

  fprintf(fp, "# =============================================================="
          "=============\n");
  fprintf(fp, "# protein-ligand hydrophobic contacts         : %.0f\n",
          features->contact_hphob_hphob);
  fprintf(fp, "# protein-ligand H-bond count                 : %d\n",
          features->number_of_hbonds);
  fprintf(fp, "# protein-ligand salt-bridge count            : %d\n",
          features->number_of_salt_bridges);
  fprintf(fp, "# metal-ligand interactions count             : %d\n",
          features->number_of_metal_ligand_bonds );
  fprintf(fp, "# unsatisfied interfacial polar atom count    : %d\n",
          features->number_of_interfacial_unsatisfied_polar_atoms );
  fprintf(fp, "# unsatisfied interfacial charged atom count  : %d\n",
          features->number_of_interfacial_unsatisfied_charged_atoms );
  fprintf(fp, "# buried carbons (x 100%%)          : %5.3f\n",
          features->buried_carbons);
  fprintf(fp, "# remaining vdW collisions         : %d\n",
          features->number_of_bumps );
  fprintf(fp, "# total vdW overlap (A)            : %5.3f\n",
          features->total_overlap);
  fprintf(fp, "# anchor fragment translations     : %d\n",
          features->num_anchor_corrections );
  fprintf(fp, "# side-chain mean-field iterations : %d\n",
          features->num_mean_field_optimizations );
  fprintf(fp, "# ligand side-chain rotations      : %d\n",
          features->num_ligand_side_chain_rotations );
  fprintf(fp, "# protein side-chain rotations     : %d\n",
          features->num_target_side_chain_rotations );
  fprintf(fp, "# Orientscore (kcal/mol), used by HBIND to select binding "
          "mode: %5.3f\n", features->orient_score );

  fprintf(fp, "# =============================================================="
          "=============\n# Anchor fragment atoms matched to template:\n" );
  for(i = 0; i < 3; i++ ){
    lig_idx = features->ligand_matches[i];
    if(features->matched_type[i] != HYDROPHOB)
      fprintf(fp, "# ligand atom  %3d  %-3s  %-5s ",
              ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].name,
              ligand_atoms[lig_idx].type_str);
    else
      fprintf(fp, "# ligand hydrophobic center %2d ", lig_idx);
    /* Switch from 0 indexed to 1 indexed numbering for template points */
    template_ref = template_interactions[features->template_matches[i]].ref + 1;
    fprintf(fp, "to template point %2d ", template_ref);
    template_act = template_interactions[features->template_matches[i]].act;
    if(template_act == DONOR) fprintf(fp, "D\n");
    else if(template_act == ACCEPTOR) fprintf(fp, "A\n");
    else if(template_act == DONEPTOR) fprintf(fp, "AD\n");
    else if(template_act == HYDROPHOB) fprintf(fp, "H\n");
    else fprintf(fp, "ERROR\n");
  }

  fprintf(fp, "# =============================================================="
          "=============\n");
  fprintf(fp, "#          |Ligand Atom -- Protein Atom|\n");
  fprintf(fp, "#          |   #  type  -- RES   # type|\n");
  for(i = 0; i < features->number_of_hbonds; i++ ){
    lig_idx = features->ligand_hbond_idz[i];
    targ_idx = features->target_hbond_idz[i];
    if(targ_idx >= 0){
      if (features->hbond_angles[i] == 0.0){
	fprintf(fp, "# hbond %3d: %3d  %-5s -- %-3s %3s %-3s "
		"(dist: %5.3f, D-H-A angle: N/A)\n", i+1,
		ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
		target_atoms[targ_idx].residue,
		target_atoms[targ_idx].residue_num, target_atoms[targ_idx].name,
		features->hbond_dists[i]);
      }else{
	fprintf(fp, "# hbond %3d: %3d  %-5s -- %-3s %3s %-3s "
		"(dist: %5.3f, D-H-A angle: %3.1f)\n", i+1,
		ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
		target_atoms[targ_idx].residue,
		target_atoms[targ_idx].residue_num, target_atoms[targ_idx].name,
		features->hbond_dists[i], features->hbond_angles[i]);
      }
    }

#ifdef SCORING_WATERS

    else

    fprintf ( fp,
		"# hbond %3d: %3d  %-5s -- water %3d (dist: %5.3f, D-H-A angle: %3.1f)\n",
		      i+1,  
		ligand_atoms[score->hbonds[i]].atom_number,
		ligand_atoms[score->hbonds[i]].type_str,
		global->waters[abs(score->hbonds[i+MAX_TOTAL_HBONDS])-1].number,
		score->hbond_dists[i],
		score->hbond_angles[i] );
#endif
  }
  for(i = 0; i < features->number_of_salt_bridges; i++){
    lig_idx = features->ligand_salt_bridge_idz[i];
    targ_idx = features->target_salt_bridge_idz[i];
    fprintf(fp, "# salt bridge %3d : %3d  %-5s -- %-3s %3s %-3s (dist: %5.3f)"
            "\n", i+1,
            ligand_atoms[lig_idx].atom_number, ligand_atoms[lig_idx].type_str,
            target_atoms[targ_idx].residue, target_atoms[targ_idx].residue_num,
            target_atoms[targ_idx].name, features->salt_bridge_dists[i]);
  }
  fprintf(fp, "# =============================================================="
          "=============\n");
}

void compute_scores(dock_feats_pt features)
{
  features->affiefficient_score =
    (features->affi_score - affinity_wts[AFFINITY_SCORE][TERM_0]) /
     features->score_component_terms[NUMBER_OF_LIGAND_NONH_ATOMS];
  features->hphob_score = affinity_wts[AFFINITY_SCORE][TERM_1] *
    features->score_component_terms[TOTAL_HPHOB_HPHOB_COMP];
  features->polar_score = affinity_wts[AFFINITY_SCORE][TERM_2] *
    (features->score_component_terms[NUMBER_OF_HBONDS] +
     features->score_component_terms[NUMBER_OF_SALT_BRIDGES] +
     features->score_component_terms[NUMBER_OF_METAL_HBONDS]);
  features->unsat_polar_score = affinity_wts[AFFINITY_SCORE][TERM_3] *
    (features->score_component_terms[NUMBER_OF_UNSAT_POLAR] +
     features->score_component_terms[NUMBER_OF_UNSAT_CHARGE]);
  features->contact_hphob_hphob =
    features->score_component_terms[CONTACT_HPHOB_HPHOB];
}
