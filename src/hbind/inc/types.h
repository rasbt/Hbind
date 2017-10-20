#ifndef _TYPES_H
#define _TYPES_H

#include <stdio.h>
#include <distance_array.h>

typedef double              transform_matrix_t[3][4];
typedef transform_matrix_t  *transform_matrix_pt;

/* this is a data type which is handy when having some value associated
   with a list of indices, and one is interested in getting the order of
   indices based on the values */
typedef struct {
  int     index1;
  int     index2;
  double  value;
} sort_array_t;

typedef sort_array_t *sort_array_pt;

/* bond description */
typedef struct {
  int                number;        /* number of bond given in the mol2 file */
  int                atom1;              /* first atom (index in atom field) */
  int                atom2;                                   /* second atom */
  int                fragment1;                /* fragment of the first atom */
  int                fragment2;            /* the fragment of the other atom */
  int                type;                            /* SINGLE, DOUBLE, ... */
  char               type_str[5];         /* string describing the bond type */
} bond_t;

typedef bond_t  *bond_pt;

/* structure describing a complete molecule */
typedef struct {
  char               name[MAX_LEN_MOL2_COMPOUND_NAME];      /* compound name */
  char               name_noconf[MAX_LEN_MOL2_COMPOUND_NAME];/* compound name without conformer info */
  char               conf_number[MAX_LEN_MOL2_COMPOUND_NAME];/* compound name without conformer info */
  atom_pt            atoms;                           /* vector of all atoms */
  bond_pt            bonds;                           /* vector of all bonds */
  int                *atom_index;  /* the internal index for the mol2 number */
  int                *flexible_bonds;  /* list of indices of rotatable bonds */
  int                *fragment_locations;  /* marker for fragments in anchor */
  int                *anchor_dist; /* the number of flexible bonds to anchor */
  int                *way_to_anchor;/* index of the flex bond back to anchor */
  int                *bond_types;    /* marker for non-flexible anchor bonds */
  int                *bond_directions;       /* shows way to anchor fragment */
  int                **neighbors;      /* adjacency list for structure graph */
  int                **bond_to_neighbor;             /* list of bond indices */
  int                *number_of_neighbors;   /* number of neighbors of atoms */
  int                *bond_order;          /* total number of bonds of atoms */
  int                **fragment_neighbors;       /* adjacency list fragments */
  int                **bond_to_fragment_neighbor;    /* list of bond indices */
  int                *number_of_fragment_neighbors;   /* number of neighbors */
  int                *relations;   /* matrix that describes neighbored atoms */
  int                *orig_act;     /* the H-bonding activity before scoring */
  float              *orig_rad;     /* the original radius (covalent bonds!) */
  float              **carbon_ring_centers;     /* positions of ring centers */
  int                *carbon_ring_atom;      /* index of an atom on the ring */
  float              *atom_positions;          /* current positions of atoms */
  float              *orig_atom_positions;  /* the positions before rotating */
  int                number_of_atoms;     /* number of atoms in the compound */
  int                number_of_bonds;               /* total number of bonds */
  int                number_of_flexible_bonds;            /* rotatable bonds */
  int                number_of_carbon_rings;       /* number of carbon rings */
  int                number_of_added_hydrogens; /* number of artificial hydr */
  int                number_of_fragments;       /* number of rigid fragments */
  /** Added by RSK - Aug 13, 2002 **/
  int                number_of_substructures;     /* number of substructures */
  char               *substructure[MAX_SUBST_LEN];          /* subsctructure */

  /** Sameer 30/Nov/04 - for recording the positions during rotamer search */
  float              **pre_rot_atom_positions; /* positions before rotamer search */
  float              **max_rot_atom_positions; /* positions for max-scoring rotamer*/
  int                *pre_rot_act;  /* H-bonding activity before rotamer search*/
  int                *max_rot_act;  /* H-bonding activity for max-scoring rotamer*/

}*molecule_pt, molecule_t;

/* entry in a rule list in the DOCK rule syntax */
typedef struct rlist {
  int                number;                      /* minimal number of atoms */
  int                type;                                  /* type of atoms */
  int                orbit;                                /* orbit of atoms */
  struct rlist       *required;                /* list of required neighbors */
  struct rlist       *prohibited;            /* list of prohibited neighbors */
}*rule_list_pt, rule_list_t;

/* header of a rule list, descibes base atom, i.e. donor/acceptor or atom
   agjacent to flexible bond */
typedef struct {
  int                type;                                   /* type of atom */
  int                orbit;                                 /* orbit of atom */
  rule_list_pt       required[MAX_RULES_PER_BONDED_ATOM];  /* req. neighbors */
  rule_list_pt       prohibited[MAX_RULES_PER_BONDED_ATOM]; /* prohb. neigh. */
}*rule_pt, rule_t;

/* definition of a flexible bond, based on DOCK syntax */
typedef struct {
  char               name[32];           /* identifier for this type of bond */
  rule_t             atom[2];                /* rules for the adjacent atoms */
  char               flex_str[8];        /* identifies flexible bond or not */
  int                flex;               /* integer identifier of flex bond */
} flex_bond_defn_t;

typedef flex_bond_defn_t  *flex_bond_defn_pt;

/* definition of the rules for hydrogen-bond donors and acceptors */
typedef struct {
  rule_t             donor_rules[MAX_NUMBER_OF_HYD_RULES];    /* donor rules */
  rule_t             acceptor_rules[MAX_NUMBER_OF_HYD_RULES];  /* acc. rules */
  int                number_of_donor_rules;               /* number of rules */
  int                number_of_acceptor_rules;            /* number of rules */
} hyd_defn_t;

typedef hyd_defn_t  *hyd_defn_pt;

/* residue table */
typedef struct {
  int       start_atom;
  int       start_interaction;
  int       number_of_atoms;
  int       type;
  char      num[6];
  char      name[4];
} residue_t;

typedef residue_t  *residue_pt;


/* interaction points for H-bonds, hydrophobility */

typedef struct {
  int       type;           /* atom type */
  int       key_point;      /* is this template point a key point */
  int       key2_point;     /* is this template point a secondary key point */
  int       act;            /* H-bonding activity */
  int       ref;            /* toneroma 26SEP06*/
  int       index;          /* index of atom in vector of type atom_t */
  float     pos[3];         /* position of the interaction point */
} interaction_t;

typedef interaction_t  *interaction_pt;

typedef struct {
  int       index;
  float     pos[3];
  char      atom_name[5];        /* atom type */
  char      residue[4];     /* residue type */
  char      residue_num[5]; /* residue number */
} hbond_buddy_t;

typedef hbond_buddy_t  *hbond_buddy_pt;

/* triangle, described by three points and the set of distances */
typedef struct {
  float               dist[3];             /* lengths of the triangle sides */
  int                 act[3];      /* interaction types of the three points */
  int                 index[3];   /* the index of the end points for a side */
} triangle_t;

typedef triangle_t  *triangle_pt;

/* hash_table entry, contains array of triangles */
typedef struct {
  triangle_pt         *entries;                /* pointer to triangle array */
  int                 number_of_entries; /* counter for number of triangles */
} hash_table_t;

typedef hash_table_t  *hash_table_pt;

typedef struct {
  char               name[MAX_COMPOUNDNAME_LEN];
  char               name_noconf[MAX_COMPOUNDNAME_LEN];
  int                first_point;
  int                number_of_points;
} pts_compound_t;

typedef pts_compound_t  *pts_compound_pt;

typedef struct {
  float              *pos;
  int                *act;
  int                *atom_index;
  pts_compound_pt    compounds;
  int                number_of_compounds;
  int                number_of_interaction_points;
} pts_t;

typedef pts_t  *pts_pt;

typedef struct {
  double               probability;
  double               angle;
  double               mean_force;
  double               force;
  double               penalty;
  float                old_distance;
  transform_matrix_pt  matrix;
} unbump_data_t;

typedef unbump_data_t  *unbump_data_pt;

typedef struct {
  int            *hbonds;
  float          *hbond_angles;
  float          *hbond_dists;
  int            *salt_bridges;
  float          *salt_bridge_dists;
  float          *score_component_terms;
  int            number_of_hbonds;
  int            number_of_salt_bridges;
  int number_of_metal_ligand_bonds;  /*!< Metal:Lig bond count */
  int number_of_interfacial_unsatisfied_polar_atoms; /*! Count of unsatisfied
                                                polar atoms in the interface */
  int number_of_interfacial_unsatisfied_charged_atoms; /*! Count of unsatisfied 
                                              charged atoms in the interface */
  int            number_of_repulsive_charge_contacts;
  int            difference_intra_target_hbonds;
  int            number_of_water_mediated_hbonds;
  float          hphob;
  float          buried_carbons;
  float          orientation_score;
  float          affinity_score;
} score_data_t;
typedef score_data_t  *score_data_pt;

/*! Struct to hold features of a docking
 *
 * Added in an attempt to reduce code when messing with bump/nobump.
 *
 * NOTE: if the dock_feats_t structure is changed to contain pointers,
 * a copy routine will need to be implemented so that a deep copy is
 * performed (would be required in match_triangles.c)
 */
typedef struct{
  /* toneroma - 20090522 - Added molecule name and conformer entries
     to keep track of the best orientation
  */
  char ligand_name[FILENAME_MAX];
  char ligand_name_noconf[FILENAME_MAX];
  char ligand_conf[FILENAME_MAX];
  int  binding_modes_counter;
  /* Matched polar atoms & geometry -- values computed in score_complex */
  int ligand_hbond_idz[MAX_TOTAL_HBONDS];
  int target_hbond_idz[MAX_TOTAL_HBONDS];
  float hbond_angles[MAX_TOTAL_HBONDS];
  float hbond_dists[MAX_TOTAL_HBONDS];
  int ligand_salt_bridge_idz[MAX_TOTAL_SALT_BRIDGES];
  int target_salt_bridge_idz[MAX_TOTAL_SALT_BRIDGES];
  float salt_bridge_dists[MAX_TOTAL_SALT_BRIDGES];
  /* Hydrophobic contacts */
  int target_hphob_contacts[MAX_PDB_ATOMS];
  /* Score terms -- computed in score_complex */
  float score_component_terms[NO_OF_SCORE_COMPONENT_TERMS]; 
  float flag;                 /*!< Might not be used */
  float orient_score;          /*!< Orientation score */
  float affiefficient_score; /*!< Affinity score divided by N */
  float affi_score;          /*!< Affinity score */
  float hydro_score;         /*!< Hydrophobic component of score */
  float hphob_score;         /*!< Hydrophobic component of score */
  float polar_score;         /*!< Hydrophobic component of score */
  float unsat_polar_score;   /*!< Unsatisfied polar compoment of score */
  float contact_hphob_hphob; /*!< Count of hydrophobic contacts */
  int number_of_hbonds;              /*!< Prot:Lig hbond count */
  int number_of_salt_bridges;        /*!< Prot:Lig salt bridge count */
  int number_of_metal_ligand_bonds;  /*!< Metal:Lig bond count */
  int number_of_interfacial_unsatisfied_polar_atoms; /*! Count of unsatisfied
                                                polar atoms in the interface */
  int number_of_interfacial_unsatisfied_charged_atoms; /*! Count of unsatisfied 
                                              charged atoms in the interface */
  float buried_carbons;  /*!< Buried carbons score */
  /* Bump items -- determined by bump checks */
  /* NOTE except for total_overlap, the rest are still kept in the global
   * structure till the global structure is extricated from the bump check
   * routines */
  float total_overlap;   /*!< Total computed overlap */
  int number_of_bumps;   /*!< Number of bumps */
  int num_anchor_corrections;
  int num_mean_field_optimizations;
  int num_side_chain_rotations;
  int num_ligand_side_chain_rotations;
  int num_target_side_chain_rotations;
  /* Anchor to template matching info  -- set in match_triangles */
  int template_matches[3];
  int ligand_matches[3];
  int matched_type[MAX_INTER_PAIRS];
}*dock_feats_pt, dock_feats_t;

/* Used to save positions of atoms for best docked orientation */
typedef struct{
  int target_rotations[MAX_PDB_RESIDUES];
  float target_positions[3*MAX_PDB_ATOMS];
  /*int ligand_rotations[MAX_NUMBER_OF_MOL2_ATOMS]; */
  float ligand_positions[3*MAX_NUMBER_OF_MOL2_ATOMS];
  int water_states[MAX_BINDING_SITE_WATERS];
  float water_positions[3*MAX_BINDING_SITE_WATERS];
}*moved_positions_pt, moved_positions_t;

/* Sameer 08/Dec/04 - struct for retaining backbone
 * dependent rotamer specific information
 * during rotamer search.
 */
typedef struct {
  float 		score;
  int			phi;
  int			psi;
  float			**atom_pos;
} rotamer_t;

typedef rotamer_t	*rotamer_pt;

/* Sameer 26/Dec/04
 * Matrix element to be used in score-based MFM probability
 * calculations. This probabiliy helps in determining the
 * order of rotamer substitution.
 */
typedef struct unsat_rotamer{
  double 	score; /* stores score delta*/
	double  mean_score;
  /*	double 	score_gain;*/
  /*double 	score_loss;*/
	int 	no_of_atoms, phi, psi, rotamer_no;
  float	**pos; /*[MAX_RESIDUE_ATOMS][3];*/
	double	probability;
	double	bb_dep_probability;
	int		residue_index;

}unsat_rotamer_t;

typedef unsat_rotamer_t* unsat_rotamer_pt;

/* Holds the parameters used to determine valid triangles used to match ligand
 * interaction points with template interaction points.
 */
typedef struct{
  double min_short_side;
  double max_short_side;
  double min_long_side;   /*!< At present this changes depending on lig size*/
  double max_long_side;
  double min_perimeter;
  double max_perimeter;
}triangle_parameters_t;


/* global information */
typedef struct {
  hash_table_pt        ***hash_table;
  char                 *hbind_dir;
  char                 *data_root;
  char                 *compound_dir;
  char                 *compound_name;		/* singleton, name_of_mm2 etc */
  char                 *compound_name_noconf;		/* singleton, name_of_mm2 etc */
  char                 *database;
  char                 *protein;
  char                 *template;
  atom_pt              target_atoms;
  float*               target_atom_positions;
  dist_array_t         target_dists_array;
  atom_pt              target_water_atoms;   /* added by Litian for waters in protein file */
                                             /* these files are got from consolv */
  interaction_pt       template_interactions;
  pts_pt               compound_interactions;
  residue_pt           target_residues;
  triangle_pt          template_triangles;
  int                  number_of_template_triangles;
  int                  number_of_triangle_pointers;
  int                  max_number_of_triangles;
  int                  number_of_hash_buckets;
  int                  score_only;
  molecule_pt          ligand;
  atom_pt              waters;
  float                *orig_target_atom_positions;
  int                  *orig_target_atom_act;
  float                **orig_water_positions;
  char                 *ligand_file_name;	/* actual ligand filename */
  char                 *old_compound_name;
  char                 *old_ligand_name_noconf;
  char                 *best_affi_name;
  int                  binding_modes_counter;
  int                  anchor_fragments[3];
  int                  *target_bumps;
  int                  *ligand_bumps;
  transform_matrix_pt  unbump_matrices;
  unbump_data_pt       **unbump_data;
  unbump_data_pt       unbump_entries;
  int                  *unbump_dependencies;
  int                  **unbump_dependents;
  int                  *number_of_unbump_dependents;
  int                  **unbump_target_indices;
  int                  *unbump_indices;
  int                  *unbump_bonds;
  int                  *target_rotations;
  int                  *target_intra_overlap;
  int                  number_of_unbump_bonds;
  int                  number_of_bumps;
  int                  number_of_water_bumps;
  int                  number_of_unbump_entries;
  int                  template_size;
  int                  template_key_points;
  int                  template_key2_points;
  int                  database_num_conformers;
  int                  database_num_exist;
  int                  compound_index;
  int                  number_of_target_atoms;
  int                  number_of_target_interactions;
  int                  number_of_target_residues;
  int                  number_of_anchor_corrections;
  int                  number_of_mean_field_optimizations;
  int                  number_of_side_chain_rotations;
  int                  number_of_ligand_side_chain_rotations;
  int                  number_of_target_side_chain_rotations;
  int                  number_of_water_translations;
  int                  number_of_waters;
  int                  number_of_orig_intra_target_hbonds;
  int                  last;
  int                  max_template_triangles;
  int                  total_num_molecules_noconf;
  int                  total_num_molecules_conf;
  int                  total_num_output_molecules;
  float                dmetol;
  float                rmstol;
  float                anchor_translation;
  float                anchor_overlap;
  float                side_chain_overlap;
  float                intra_overlap;
  float                intermediate_overlap;
  float                finally_tolerated_overlap;
  float                finally_tolerated_max_bump;
  float                dist_matrix_err;
  float                rms_err;
  float                lowest_overlap;
  float                best_affiscore_so_far;
  float                affiscore_sums;
  float                affiscore_mean;
  float                affiscore_squared_sums;
  float                affiscore_svar;
  float                affiscore_stdd;
  float                affiscore_significance;
  dock_feats_t current_orientation;     /*!< Struct for current docking */
  dock_feats_t best_orientation;
  moved_positions_t best_orient_positions;
#if 0
  dock_feats_t best_orientation_bump;   /*!< Features of best docking */
  dock_feats_t best_orientation_nobump; /*!< Features of best bump free 
                                             docking */
#ifndef OUTPUT_ALL_MATCHES
  moved_positions_t positions_bump;
  moved_positions_t positions_nobump;
#endif
#endif 
  float                score_cutoff;
  flex_bond_defn_pt    flex_bond_rules;
  hyd_defn_t           hyd_atom_rules;
  int                  number_of_flex_bond_rules;
  int                  number_of_screened_compounds;
  int                  number_of_potential_ligands;
  int                  filter_counter[NUMBER_OF_FILTERS];
  float                *target_ligand_distances;
  float		       *target_ligand_sq_dists;
  float                *intra_target_distances;            /* added by Litian 07/23/2003 */
  float                *target_water_distances;            /* added by Litian 07/23/2003 */
  float                *ligand_water_distances;  /* Sameer - used in scoring 22/Nov/2003 */
  float                *intra_water_distances;   /* Sameer - used in scoring 22/Nov/2003 */
/***** Added by PCS -- 22-Mar-00 *****/
  int                  number_of_bumped_ligands;
/*************************************/

#ifdef CORRECT_SF_L_TERM  /* Sameer 16/May/05 - to be used for scoring func term 'L' */
  int 					*metal_atom_indices;    
  int 					number_of_metals; 
#endif

/********* Datastructures for rotamer sampling - Sameer 30/Nov/04 *************/
  int		           *unsat_residue_indices;
  int				   number_of_unsat_residues;      /* used to record unsat sidechains during scoring */
  /* 27/Dec/04 Sameer - Datastructures for storing high-scoring rotamers */
  unsat_rotamer_pt 		*unsat_residue_rotamers;
  /* # ifdef SCORE_MF_MAXIMIZATION*/
  int					no_of_unsat_SCMF_residues;
  /* # endif*/
  int					*no_of_top_rotamers;

  /* variables to store score at various points in time */
  float					pre_rot_srch_score;
  float 				pre_scmf_unbump_score;
  float					post_rot_srch_score;

  /* flag-arrays used the in the scoring funciton */
  short					*ligand_flag, 
					*target_flag;

/********* Datastructures for interaction opportunities- Sameer 28/Mar/05 *************/
#ifdef INTERACTION_OPPORTUNITIES


  int                   *target_interaction_opps;
  int                   *ligand_interaction_opps;
  int                   number_of_interaction_opps;
  int                   number_of_interaction_opps_bonds;
  int                   max_interaction_opportunity_distance;

#endif
/**************************************************************************/

  FILE                  *MM2_FILE;   /* file handle to current multimol2 */
  triangle_parameters_t triangle_parameters;
  int                   match_2_key_points;
  int                   group_conformers;
  int                   restart_molecule_check;
  char                  restart_molecule[MAX_COMPOUNDNAME_LEN];
  double                small_ligand_diameter;
}global_data_t;

typedef global_data_t  *global_data_pt;

typedef double   point_t[3];
typedef point_t  *point_pt;

typedef double   vector_t[3];
typedef vector_t *vector_pt;

typedef struct {
  double     protein_atom_pos[3];        /* atom coordinates */
  double     ligand_atom_pos[3];         /* atom coordinates */
  double     protein_atom_rad;           /* atom radius */
  double     ligand_atom_rad;            /* atom radius */
  double     overlap;
} bump_t;

typedef bump_t  *bump_pt;

#endif
