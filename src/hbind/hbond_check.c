#include <math.h>
#include <string.h>
#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "dist_fun.h"
#include "err_handle.h"
#include "trans_rotate.h"


#define STRAIGHT_IN_PLANE    0
#define TWO_POSSIBILITIES    1
#define POINTING_TO_ACCEPTOR 2
#define TETRAHEDRAL          3
#define TETRAHEDRAL_3_NEIGHBORS          4
#define TETRAHEDRAL_2_NEIGHBORS          5

float  compute_angle ( float  tail[],
		       float  joint[],
		       float  head[] )
{
  double vector1[3],
         vector2[3],
         angle;
  int    i;

  for ( i = 0; i < 3; i++ )
    {

      vector1[i] = (double) ( joint[i] - tail[i] );

      vector2[i] = (double) ( head[i] - joint[i] );
    }

  angle = acos ( ( vector1[0] * vector2[0] 
		   + vector1[1] * vector2[1]
		   + vector1[2] * vector2[2] )
		 / sqrt ( vector1[0] * vector1[0]
			  + vector1[1] * vector1[1]
			  + vector1[2] * vector1[2] )
		 / sqrt ( vector2[0] * vector2[0]
			  + vector2[1] * vector2[1]
			  + vector2[2] * vector2[2] ) );
  /* transform to deg */
  angle = 180.0 - ( angle * 180.0 / M_PI );
  return (float) angle;
}

float  compute_ligand_h_angle ( molecule_pt ligand,
				int         donor_index,
				float       *acceptor_position,
				int         hydrogen,
				global_data_pt global)
{
  transform_matrix_t transformation_matrix;
  atom_pt            atoms;
  float              xpoint[3],
                     ypoint[3],
                     hpoint[3],
                     hpoint2[3],
                     hpoint3[3],
                     pos[3];
  float              *hydrogen_position, *hydrogen_position2, *hydrogen_position3;
  double             help;
  float              bond_length,
                     angle;
  int                *neighbors;
  int                which_case,
                     number_of_hydrogens,
                     count;
  int                i, j;
  char               err_msg[FILENAME_MAX];

  bond_length = 1.0;

  atoms = ligand->atoms;
  hydrogen_position = atoms[hydrogen].pos;
  switch ( atoms[donor_index].type  )
    {
    case N:
      switch ( atoms[donor_index].orbit )
	{
	case AR:
	case AM:
	case PL3:
	case SP2:
	  angle = 120.0;
	  which_case = STRAIGHT_IN_PLANE;	  
	  break;
	case SP3:
	case SP4:
	  angle = 109.5;
	  which_case = TETRAHEDRAL;
	  break;
	default:
	  sprintf ( err_msg, 
		    "unknown orbit: %d (%s)\n", 
		    atoms[donor_index].orbit, ligand->name );
	  fprintf ( stderr, err_msg);
	  err_print ( err_msg);
	  err_warning2 ( "compute_ligand_h_angle",
			"unknown N donor orbit, results unreliable");
	  break;
	}
      bond_length = 1.01;
      break;
    case O:
      if ( atoms[donor_index].orbit != SP3 )
	{
	  sprintf ( err_msg, 
		    "atom %d wrong orbit: %d (%s)\n", 
		    donor_index + 1,
		    atoms[donor_index].orbit,
		    ligand->name );
	  fprintf ( stderr, err_msg);
	  err_print (err_msg);
	  err_warning2 ( "compute_ligand_h_angle",
			"no O donor orbit, results unreliable");
	}
      else
	{
	  angle = 109.5;
	  which_case = TETRAHEDRAL;
	  bond_length = 0.96;	  
	}
      break;
    default:
      sprintf ( err_msg, "WARNING: ligand %s\n", ligand->name);
      err_print (err_msg);
      err_warning2 ( "compute_ligand_h_angle",
		    "unknown hydrogen bond donor type, results unreliable");
      break;
    }
  neighbors = ligand->neighbors[donor_index];
  count = 0;
  number_of_hydrogens = 0;
  for(i = 0; i < 3; i++){ 
    ypoint[i] = 0.0;
    xpoint[i] = atoms[neighbors[0]].pos[i]; 
    hpoint[i] = 0.0;
    hpoint2[i] = 0.0;
    hpoint3[i] = 0.0;
  }
  for(i = 0; i < ligand->number_of_neighbors[donor_index]; i++ )
    if(atoms[neighbors[i]].type != H ){ 
      for(j = 0; j < 3; j++)ypoint[j] += atoms[neighbors[i]].pos[j];
      count++;
    }else{ 
      for(j = 0; j < 3; j++ ){ 
        if(number_of_hydrogens < 1) hpoint[j] = atoms[neighbors[i]].pos[j]; 
        else if (number_of_hydrogens < 2){ 
          hydrogen_position2 = atoms[neighbors[i]].pos; 
          hpoint2[j] = atoms[neighbors[i]].pos[j]; 
        }else if(number_of_hydrogens < 3){ 
          hydrogen_position3 = atoms[neighbors[i]].pos; 
          hpoint3[j] = atoms[neighbors[i]].pos[j]; 
        } 
      }
      number_of_hydrogens++;	
    }

  if(count > 2 && atoms[donor_index].orbit != SP4) /* toneroma 06MAR07 - a positively charged Nitrogen (N.4) may have more than 2 non-H neighbors*/
    {
      sprintf ( err_msg, "WARNING: ligand %s\n", ligand->name);
      err_print (err_msg);
      sprintf ( err_msg, "donor (#%d: %s) has more than two non-H neighbors",
		atoms[donor_index].atom_number, atoms[donor_index].name);
      err_warning2 ( "compute_ligand_h_angle", err_msg);
      which_case = TETRAHEDRAL_3_NEIGHBORS;
    }
  if(count > 1){ 
    for(i = 0; i < 3; i++){

      ypoint[i] -= xpoint[i];
      ypoint[i] /= (count - 1);
    }     
    which_case = TETRAHEDRAL_2_NEIGHBORS;


  }else{
    for ( i = 0; i < 3; i++ ) xpoint[i] = acceptor_position[i];
    which_case = POINTING_TO_ACCEPTOR;
  }

  compute_transformation_matrix ( xpoint,                  /* X-axis */
				  ypoint,                  /* Y-axis */
				  atoms[donor_index].pos,   /* origin */
				  transformation_matrix );
  switch ( which_case )
    {
    case POINTING_TO_ACCEPTOR:
      /* acceptor position is in the positive X-part of the XY-plane */
      if (number_of_hydrogens > 0)
	{
	  hydrogen_position[X] = 
	    bond_length * sin ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Y] = 
	    bond_length * cos ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Z] = 0.0;
	}
      if (number_of_hydrogens > 1)
	{
	  hydrogen_position2[X] = 
	    -( bond_length * sin ( ( 180.0 - angle ) / 180.0 * M_PI ) ) * cos ( ( 60.0 / 180.0 ) * M_PI );
	  hydrogen_position2[Y] = 
	    bond_length * cos ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position2[Z] =
	    ( bond_length * sin ( ( 180.0 - angle ) / 180.0 * M_PI ) ) * sin ( ( 60.0 / 180.0 ) * M_PI );
	}
      if (number_of_hydrogens > 2)
	{
	  hydrogen_position3[X] =
	    -( bond_length * sin ( ( 180.0 - angle ) / 180.0 * M_PI ) ) * cos ( ( 60.0 / 180.0 ) * M_PI );
	  hydrogen_position3[Y] =
	    bond_length * cos ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position3[Z] =
	    -( bond_length * sin ( ( 180.0 - angle ) / 180.0 * M_PI ) ) * sin ( ( 60.0 / 180.0 ) * M_PI );
	}
      break;
    case TETRAHEDRAL:
      for ( i = 0; i < 3; i++ ) pos[i] = acceptor_position[i];
      transform_point(pos, transformation_matrix );
      transform_point(xpoint,  transformation_matrix );
      
      hydrogen_position[X] = 0.0;
      hydrogen_position[Y] 
	= bond_length * sqrt ( xpoint[X] * xpoint[X] + xpoint[Y] * xpoint[Y] )
	* cos ( ( 180.0 - angle ) * M_PI / 180.0 ) / xpoint[Y];
      help = 1.0 
	- cos ( ( 180.0 - angle ) * M_PI / 180.0 )
	- ( xpoint[X] * xpoint[X]
	    * cos ( ( 180.0 - angle ) * M_PI / 180.0 )
	    / xpoint[Y] * xpoint[Y] );
      if ( help >= 0.0 )
	hydrogen_position[Z]  = bond_length * sqrt ( help );
      else
	{
	  sprintf ( err_msg, "WARNING: ligand %s\n", ligand->name);
	  err_print (err_msg);
	  err_warning2 ( "compute_ligand_h_angle",
			 "sqrt(0) -> computed hydrogen position incorrect");
	  hydrogen_position[Z]  = 0.0;
	}
      if ( pos[Z] < 0.0 )

	hydrogen_position[Z] *= (-1.0);
      break;
    case STRAIGHT_IN_PLANE:
      hydrogen_position[X] = 0.0;
      hydrogen_position[Y] = bond_length; 
      hydrogen_position[Z] = 0.0;
      break;
    case TETRAHEDRAL_2_NEIGHBORS:
      for ( i = 0; i < 3; i++ )
	pos[i] = acceptor_position[i];
      transform_point(pos, transformation_matrix);
      transform_point(xpoint, transformation_matrix);
      transform_point(hpoint, transformation_matrix);
      
      if (number_of_hydrogens > 0)
	{
	  hydrogen_position[X] = hpoint[X];
	  hydrogen_position[Y] = hpoint[Y];
	  hydrogen_position[Z] = hpoint[Z];
	  if ( pos[Z] < 0.0 )
	    hydrogen_position[Z] *= (-1.0);
	}
      if (number_of_hydrogens > 1)
	{
	  hydrogen_position2[X] = hpoint[X];
	  hydrogen_position2[Y] = hpoint[Y];
	  hydrogen_position2[Z] = 0.0 - hpoint[Z];
	  if ( pos[Z] < 0.0 )
	    hydrogen_position2[Z] *= (-1.0);
	}
      break;
    case TETRAHEDRAL_3_NEIGHBORS:
      hydrogen_position[X] = hpoint[X];
      hydrogen_position[Y] = hpoint[Y];
      hydrogen_position[Z] = hpoint[Z];
      break;
    default:
      break;
    }
  /* transform hydrogen position back to world coordinate system */
  if(number_of_hydrogens > 0)
    transform_point_back(hydrogen_position, transformation_matrix );      
  if(number_of_hydrogens > 1)
    transform_point_back(hydrogen_position2, transformation_matrix );      
  if(number_of_hydrogens > 2)
    transform_point_back(hydrogen_position3, transformation_matrix );      

  angle = compute_angle ( atoms[donor_index].pos,
			  hydrogen_position,
			  acceptor_position );
  return angle;
}

/* get positions of the atom defined by "*types" in the current residue */
int get_positions(atom_pt residue_atoms, int number_of_residue_atoms,
                  int *types, int number, float positions[3][3])
{
  int  found[3];
  int  i, index;
  int  matched = 0;
  
#ifdef TRACE
  printf ("get_positions():\n");
  printf("first_atom_in_residue=%d, number_of_residue_atoms=%d, types= %d,"
         " number=%d\n", first_atom_in_residue, number_of_residue_atoms, 
         *types, number); 
  printf("atom[%d].pos= ( %7.3f %7.3f %7.3f )\n", first_atom_in_residue,
         atoms[first_atom_in_residue].pos[X], 
         atoms[first_atom_in_residue].pos[Y], 
         atoms[first_atom_in_residue].pos[Z]);
#endif
  
  /* atoms might have multiple positions, we use the found flag to ensure we
     get only the first one */
  
  /* look for the number of atoms (number) with the given type (types) in the 
   * given residue (residue_atoms) */
  for(i = 0; i < 3; i++) found[i] = FALSE;
  for(index = 0; index <  number_of_residue_atoms; index++){ 
    for(i = 0; i < number; ++i){
      if(!found[i] && residue_atoms[index].type == types[i]){
        memcpy(positions[i], residue_atoms[index].pos, 3*sizeof(float));
        found[i] = TRUE;
        matched++;
      }
    }
  }

#ifdef TRACE
  printf("position = ( %7.3f %7.3f %7.3f )\n", positions[0][0], positions[0][1],
	 positions[0][2] );
#endif
  
  /* get error if didn't find  */
  if(matched == number){
    return SUCCESS;
  }else{

#ifdef TRACE  
    fprintf(stderr, "matched = %d\nlooking for types %d", matched, types[0]);
    if(number > 1) fprintf(stderr, ", %d", types[1]); 
    if(number > 2) fprintf(stderr, ", %d", types[2]);
    fprintf(stderr, "\n" ); 
    for(i = 0; i < number_of_residue_atoms; i++) 
      fprintf(stderr, "%d: %2d %3s %4s %4s %7.3f %7.3f %7.3f\n", i, 
              atoms[first_atom_in_residue+i].type, 
              atoms[first_atom_in_residue+i].residue,
              atoms[first_atom_in_residue+i].residue_num,
              atoms[first_atom_in_residue+i].name,
              atoms[first_atom_in_residue+i].pos[0],
              atoms[first_atom_in_residue+i].pos[1],
              atoms[first_atom_in_residue+i].pos[2]);
       err_panic("get_position", "target atom not found");
#endif

    return FAILURE;
  }
  return SUCCESS;
}

 /* get position of atoms that form the plane:
        for N-terminus, "CA", "N" and the acceptor
        for none n-terminus, "CA", N" and "C in the previous residue 
        for side-chain, ... */
float  compute_target_hydrogen_angle ( atom_pt    atoms,
				       residue_pt residues,
				       int        index,
				       float      *acceptor_position,
				       float      *hydrogen_position,
				       global_data_pt global)
{
  transform_matrix_t transformation_matrix;
  float              points[3][3];
  float              /*hydrogen_position[3],*/
                     pos[3];
  float              length,
                     angle;
  int                atom_types[3];
  int                which_case,
                     residue_index;
  int                i,
                     find_plane; 
  char               err_msg[FILENAME_MAX];


  if ( atoms[index].res == HETATM )
    /* this is a heteroatom, so we cannot compute the hydrogen position */
    return 175.0;
  residue_index = atoms[index].residue_index;
  find_plane = 1;

#ifdef TRACE
  printf("\ncompute_target_hydrogen_angle ():\n");
  printf("residue_index = %d, start_atom = %d, number_of_atoms = %d\n",
	 residue_index, residues[residue_index].start_atom,
	 residues[residue_index].number_of_atoms );
  printf("atom[%d].type=%d, atom_pos=[%7.3f, %7.3f, %7.3f]\n",
	 index, atoms[index].type, 
	 atoms[index].pos[X], atoms[index].pos[Y], atoms[index].pos[Z]);
#endif

  if ( atoms[index].type == N )
    /* main-chain nitrogen */
    {
      if ( residue_index == 0 )
	/* N-terminus, no previous residue, plane defined by CA, N, and
	   acceptor position */
	{
	  atom_types[0] = CA;
	  atom_types[1] = N;
	  /*printf("N-terminus\n");*/
	  /* get position of "CA" and "N" in the current residue, 
	     and save into positions[0][i], positions[1][i] */
	  get_positions(&atoms[residues[residue_index].start_atom], 
                        residues[residue_index].number_of_atoms, atom_types, 2,
                        points);

	  for ( i = 0; i < 3; i++ )
	    points[2][i] = acceptor_position[i];
	  angle = 114;
	  length = 1.02;
	  which_case = POINTING_TO_ACCEPTOR;
	}
      else
	/* no N-terminus, plane defined by N and CA in current residue
	   and C in previous residue */
	{
	  atom_types[0] = CA;
	  atom_types[1] = N;
	  find_plane = get_positions(&atoms[residues[residue_index].start_atom],
				       residues[residue_index].number_of_atoms,
				       atom_types,
				       2,
				       points );
	  if ( find_plane != FAILURE )
	      {
		  /* move CA to last position in vector points */
		  /*printf("no N-terminus\n");*/
		  for ( i = 0; i < 3; i++ )
		      points[2][i] = points[0][i];
		  /* search C in previous residue */
		  atom_types[0] = C;
		  find_plane = get_positions(&atoms[residues[residue_index-1].start_atom],
					       residues[residue_index-1].number_of_atoms,
					       atom_types,
					       1,
					       points );
	      }
	  else  
	      {
		  atom_types[0] = C;
		  atom_types[1] = O;
		  find_plane = get_positions(&atoms[residues[residue_index].start_atom],
					       residues[residue_index].number_of_atoms,
					       atom_types,
					       2,
					       points );
		  if ( find_plane != FAILURE )
		    {
		      /* move O to last position in vector points */
		      for ( i = 0; i < 3; i++ )
			points[2][i] = points[1][i]; 
		      /* move N to 2nd position in vector points */
		      for ( i = 0; i < 3; i++ )
			points[1][i] = atoms[index].pos[i]; 			   
		    }
	      }
	  length = 1.02;
	  which_case = STRAIGHT_IN_PLANE;
	}
    }
  else
    /* side-chain donor */
    switch ( atoms[index].res )
      {
      case ARG:
	switch ( atoms[index].type )
	  {
	  case NE:
	    /* axis is avg(CD,CZ),NE, point is CZ (just one of the
	       three, since everything is planar, hydrogen position
	       if independent of acceptor position) */
	    atom_types[0] = CD;
	    atom_types[1] = NE;
	    atom_types[2] = CZ;
	    find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			    residues[residue_index].number_of_atoms,
			    atom_types,
			    3,
			    points );
	    length = 1.02;
	    which_case = STRAIGHT_IN_PLANE;
	    break;
	  case NH1:
	    /* axis is CZ, NH1, point is NH2 */
	    atom_types[0] = CZ;
	    atom_types[1] = NH1;
	    atom_types[2] = NH2;
	    find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			    residues[residue_index].number_of_atoms,
			    atom_types,
			    3,
			    points );
	    length = 1.02;
	    angle = 117.0;
	    which_case = TWO_POSSIBILITIES;
	    break;
	  case NH2:
	    /* axis is CZ, NH2, point is NH1 */
	    atom_types[0] = CZ;
	    atom_types[1] = NH2;
	    atom_types[2] = NH1;
	    find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			    residues[residue_index].number_of_atoms,
			    atom_types,
			    3,
			    points );
	    angle = 117.0;
	    length = 1.02;
	    which_case = TWO_POSSIBILITIES;
	    break;
	  default:
	    break;
	  }
	break;
      case ASN:
	/* axis is CG, ND2, point is OD1 */
	atom_types[0] = CG;
	atom_types[1] = ND2;
	atom_types[2] = OD1;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			3,
			points );
	angle = 115.0;
	length = 1.02;
	which_case = TWO_POSSIBILITIES;
	break;
      case GLN:
	/* axis is CD, NE2, point is OE1 */
	atom_types[0] = CD;
	atom_types[1] = NE2;
	atom_types[2] = OE1;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			3,
			points );
	angle = 115.0;
	length = 1.02;
	which_case = TWO_POSSIBILITIES;
	break;
      case HIS:
	if ( atoms[index].type == ND1 )
	  /* axis is avg(CG,CE1), ND1, point is CE1 (could be any, since
	     ring is planar) */
	  {
	    atom_types[0] = CG;
	    atom_types[1] = ND1;
	    atom_types[2] = CE1;
	    find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			    residues[residue_index].number_of_atoms,
			    atom_types,
			    3,
			    points );
	  }
	else
	  /* axis is avg(CD1, CE1), NE2, point is CE1 (could be any, since
	     ring is planar) */
	  {
	    atom_types[0] = CE1;
	    atom_types[1] = NE2;
	    atom_types[2] = CD2;
	    find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			    residues[residue_index].number_of_atoms,
			    atom_types,
			    3,
			    points );
	  }
	length = 1.02;
	which_case = STRAIGHT_IN_PLANE;
	break;
      case LYS:
	/* axis is CE, NZ, point is acceptor position */
	atom_types[0] = CE;
	atom_types[1] = NZ;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			2,
			points );
	for ( i = 0; i < 3; i++ )
	  points[2][i] = acceptor_position[i];
	angle = 114.0;
	length = 1.02;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case SER:
	/* axis is CB, OG, point is acceptor position */
	atom_types[0] = CB;
	atom_types[1] = OG;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			2,
			points );
	for ( i = 0; i < 3; i++ )
	  points[2][i] = acceptor_position[i];
	angle = 107.0;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case THR:
	/* axis is CB, OG1, point is acceptor position */
	atom_types[0] = CB;
	atom_types[1] = OG1;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			2,
			points );
	for ( i = 0; i < 3; i++ )
	  points[2][i] = acceptor_position[i];
	angle = 106.0;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      case TRP:
	/* axis is avg(CD1, CE2), NE1, point is CE2 (could be any, since
	   ring is planar) */
	atom_types[0] = CD1;
	atom_types[1] = NE1;
	atom_types[2] = CE2;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			3,
			points );
	length = 1.02;
	which_case = STRAIGHT_IN_PLANE;
	break;
      case TYR:
	/* axis is CZ, OH, point is acceptor position */
	atom_types[0] = CZ;
	atom_types[1] = OH;
	find_plane = get_positions(&atoms[residues[residue_index].start_atom],
			residues[residue_index].number_of_atoms,
			atom_types,
			2,
			points );
	for ( i = 0; i < 3; i++ )
	  points[2][i] = acceptor_position[i];
	angle = 108.0;
	length = 0.96;
	which_case = POINTING_TO_ACCEPTOR;
	break;
      default:
	sprintf ( err_msg, 
		  "atom %d: residue=%s residue_num=%s atom_name=%s type=%d act=%d\n", 
		  index,
		  atoms[index].residue,
		  atoms[index].residue_num,
		  atoms[index].name,
		  atoms[index].type,
		  atoms[index].act );
	fprintf ( stderr, err_msg);
	err_print (err_msg);
	err_panic2 ( "compute_target_hydrogen_angle", "unknown residue");
	break;
      }

  if ( find_plane != FAILURE )
    {
#ifdef TRACE
      printf("angle = %f, length = %f, atom_types[0] = %d, atom_types[1] = %d\n",
	     angle, length, atom_types[0], atom_types[1] );
#endif

      /* get the transform matrix to tranform the world coordinates
     system to a new corrdinates using points[1] as the orgin */

      compute_transformation_matrix ( points[2],                /* X-axis */
				      points[0],                /* Y-axis */
				      points[1],                /* origin */
				      transformation_matrix );
      /*
      for(i = 0; i < 3; i++) printf("%f %f %f\n", points[i][0], points[i][1], points[i][2]);
      for(i = 0; i < 3; i++){
        memcpy(pos, points[i], 3*sizeof(points[0][0]));
        transform_point(pos, transformation_matrix ); 
        printf("%f %f %f\n", pos[0], pos[1], pos[2]);
      }
      */

      switch ( which_case )
	{
	case POINTING_TO_ACCEPTOR:
	  /* acceptor position is in the positive X-part of the XY-plane */
	  hydrogen_position[X] = length * sin ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Y] = length * cos ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Z] = 0.0;
	  break;
	case TWO_POSSIBILITIES:
	  /* two positions for the hydrogen are possible, both on the XY-plane
	     of the transformed system, one with positive X ordinate, the other
	     is it's mirror image with regard to the Y-axis */
	  /* better don't mess around with the original acceptor entry... */
	  for ( i = 0; i < 3; i++ )
	    pos[i] = acceptor_position[i];
	  transform_point ( pos,
			    transformation_matrix );
	  /* default is 'right' hydrogen position */
	  hydrogen_position[X] 
	    = length * sin ( ( 180.0 - angle ) * M_PI / 180.0 );
	  hydrogen_position[Y] 
	    = length * cos ( ( 180.0 - angle ) * M_PI / 180.0 );
	  if ( pos[X] < 0.0 )
	    /* projection of accceptor is in XY-plane left of Y-axis, so
	       choose 'left' hydrogen position */
	    hydrogen_position[X] *= (-1.0);
	  hydrogen_position[Z] = 0.0;
	  break;
	case STRAIGHT_IN_PLANE:
	  angle = compute_angle ( points[0],
				  points[1],
				  points[2] );
	  angle = ( 360.0 - angle ) / 2.0;
	  /* Hydrogen is located on the XY-plane just between points[0]
	     and points[2], points[1] is the origin of the transformed 
	     system */	 
	  hydrogen_position[X] = length * sin ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Y] = length * cos ( ( 180.0 - angle ) / 180.0 * M_PI );
	  hydrogen_position[Z] = 0.0;
	  for ( i = 0; i < 3; i++ )
	    pos[i] = points[2][i];
	  transform_point ( pos,
			    transformation_matrix );
	  if ( pos[X] > 0.0 )
	    /* projection of the third point is in XY-plane with a positive
	       X-ordinate, thus the hydrogen needs a position with negative
	       X-ordinate */
	    hydrogen_position[X] *= (-1.0);
	  break;
	default:
	  break;
	}
      /* transform hydrogen position back to world coordinate system */

      transform_point_back ( hydrogen_position, 
			     transformation_matrix ); 

      /*      printf("\ntoneroma2 Hx: %f Hy: %f Hz: %f\n", hydrogen_position[X], hydrogen_position[Y], hydrogen_position[Z]);  toneroma 29MAY06 - debugging*/
      /*      printf("\ntoneroma3 Dx: %f Dy: %f Dz: %f\n", atoms[index].pos[X], atoms[index].pos[Y], atoms[index].pos[Z]);  toneroma 29MAY06 - debugging*/
      /*      printf("\ntoneroma2 Ax: %f Ay: %f Az: %f\n", acceptor_position[X], acceptor_position[Y], acceptor_position[Z]);  toneroma 29MAY06 - debugging*/

      /* calculate the angle of N-H-Acceptor */
      angle = compute_angle(atoms[index].pos, hydrogen_position,
			    acceptor_position );
    }

  else angle = FAILURE_ANGLE;
  return angle;
}
				 

int  find_index ( atom_pt atoms,
		  int     first_atom_in_residue,
		  int     number_of_residue_atoms,
		  int     type,
		  global_data_pt global)
{
  int  i = 0;
  char err_msg[FILENAME_MAX];
  
  while ( i < number_of_residue_atoms 
	  && atoms[first_atom_in_residue+i].type != type )
    i++;
  if ( i == number_of_residue_atoms )
    {
      sprintf ( err_msg, "looking for type %d:\n", type );
      err_print (err_msg);
      fprintf ( stderr, err_msg);
      for ( i = 0; i < number_of_residue_atoms; i++ )
	{
	  sprintf ( err_msg,
		    "%d: %d %s %s %s\n", 
		    i, 
		    atoms[first_atom_in_residue+i].type,
		    atoms[first_atom_in_residue+i].residue,
		    atoms[first_atom_in_residue+i].residue_num,
		    atoms[first_atom_in_residue+i].name );
	  err_print(err_msg);
	  fprintf ( stderr, err_msg);
	}
      err_panic2 ( "find_index", "target atom not found");
    }
  return first_atom_in_residue + i;
}


int  check_hydrogen_angle ( global_data_pt global,
			    int            target_index,
			    int            ligand_index,
			    float          *angle )
{
  transform_matrix_t transformation_matrix;
  float              xpoint[3],
                     ypoint[3],
                     hydrogen_position[3];
  float              bond_length;
  int                hydrogen_indices[3];
  int                *neighbors;
  int                count,
                     number_of_hydrogens,
                     number_of_added_hydrogens;
  int                i, j;
  float              h_position[3], target_neighbor_pos[3];
  int                atom_types[3];
  float              points[3][3];

  int rv = FAILURE;
  int residue_index = global->target_atoms[target_index].residue_index;
  int preacceptor_pass = FALSE;
  float pre_angle = 0.0;
  int find_target_neighbor = 1;
  molecule_pt ligand = global->ligand;
  atom_pt target_atom = &global->target_atoms[target_index];

#ifdef TRACE
  printf("check_hydrogen_angle:\n");
  printf("ligand atom %d, target atom %d: ", ligand_index, target_index);
  printf("ligand->atoms[ligand_index].act = %d  ",
         ligand->atoms[ligand_index].act);
  printf("target_atom->act = %d\n", target_atom->act);
#endif
  
  for(i = 0; i < 3; i ++ ){
    hydrogen_indices[i] = -1;
    hydrogen_position[i] = 0.0;
    h_position[i] = 0.0;
    target_neighbor_pos[i] = 0.0;
  }
  
  if ( ligand->atoms[ligand_index].act == ACCEPTOR
       || ( ligand->atoms[ligand_index].act == DONEPTOR
	    && ( target_atom->act == DONOR 
		 || target_atom->act == METAL_1 
		 || target_atom->act == METAL_2 ) ) )
    /* the target is the donor, this is the easier case */
    {
      *angle = 
	compute_target_hydrogen_angle ( global->target_atoms,
					global->target_residues,
					target_index,
					ligand->atoms[ligand_index].pos,
					hydrogen_position,
					global);

      /* code to check pre-acceptor angle */
      neighbors = ligand->neighbors[ligand_index];
      for ( i = 0; i < ligand->number_of_neighbors[ligand_index]; i++ )
	{
	  if ( ligand->atoms[neighbors[i]].type != H )
	    {	
	      pre_angle = compute_angle ( ligand->atoms[neighbors[i]].pos,
					  ligand->atoms[ligand_index].pos,
					  /*global->target_atoms[target_index].pos );*/
					  hydrogen_position );
	      if (pre_angle >= MIN_PREACCEPTOR_ANGLE && pre_angle <= MAX_PREACCEPTOR_ANGLE)
		{
		  preacceptor_pass = TRUE;
		}
	      if ( preacceptor_pass == TRUE) break;
	    }     
	}

  }else{

    /* first get all hydrogens and check, if they already have a
       position, or if they were added */
    number_of_hydrogens = 0;
    number_of_added_hydrogens = 0;
    neighbors = ligand->neighbors[ligand_index];
#ifdef TRACE
    printf("----- ligand_index=%d (%7.3f, %7.3f, %7.3f) \n", ligand_index, 
           ligand->atoms[ligand_index].pos[X], 
           ligand->atoms[ligand_index].pos[Y], 
           ligand->atoms[ligand_index].pos[Z]);
#endif

    for(i = 0; i < ligand->number_of_neighbors[ligand_index]; i++){
#ifdef TRACE
      printf("       neighbors[%d]=%d type %s", i, neighbors[i], 
             ligand->atoms[neighbors[i]].type_str );
#endif
      if(ligand->atoms[neighbors[i]].type == H){		      
        hydrogen_indices[number_of_hydrogens++] = neighbors[i];
#ifdef TRACE
        printf("-- found neighbors[%d] = hydrogen, hydrogen_indices[%d] = %d\n",
               i, number_of_hydrogens-1, 
               hydrogen_indices[number_of_hydrogens-1] );
        printf("      number_of_hydrogens = %d\n", number_of_hydrogens );
#endif
        if(strcmp(ligand->atoms[neighbors[i]].type_str, "H_ADD" ) == 0) 
          number_of_added_hydrogens++; 
      }
#ifdef TRACE
      else printf( "-- found neighbors[%d] != hydrogen \n", i );
#endif
    }

    /* keep previously written hydrogen position for N.ar, am, pl3, and sp2 
     * orbitals */
    if(ligand->atoms[ligand_index].orbit == AR || 
       ligand->atoms[ligand_index].orbit == AM || 
       ligand->atoms[ligand_index].orbit == PL3 || 
       ligand->atoms[ligand_index].orbit == SP2){
      /* all hydrogens already have positions, so check, if there is
         one in a correct position to form a H-bond */
      /* check if the donor is connected via a rotatable bond */
      neighbors = ligand->neighbors[ligand_index];
      count = 0;
      for( i = 0; i < 3; i++ ){
        ypoint[i] = 0.0;
        xpoint[i] = ligand->atoms[neighbors[0]].pos[i];
      }
      for ( i = 0; i < ligand->number_of_neighbors[ligand_index]; i++ )
        if ( ligand->atoms[neighbors[i]].type != H ){
          /* Why are we adding? */
          for ( j = 0; j < 3; j++ ) 
            ypoint[j] += ligand->atoms[neighbors[i]].pos[j]; 
            count++; 
        }else{
        /* Presumedly if we have more than 1 H we don't care which one */
          for(j = 0; j < 3; j++)
            hydrogen_position[j] = ligand->atoms[neighbors[i]].pos[j];
        }

        /* only one neighbor, so we just assume that the hydrogen 
           is pointing towards the acceptor */
        if ( count == 1 ){
          for ( i = 0; i < 3; i++ ) xpoint[i] = target_atom->pos[i];
          bond_length = dist_fun(hydrogen_position, 
                                 ligand->atoms[ligand_index].pos);
          *angle = compute_angle(ypoint, ligand->atoms[ligand_index].pos, 
                                 hydrogen_position);
          compute_transformation_matrix(xpoint,              /* X-axis */
			  ypoint,              /* Y-axis */
			  ligand->atoms[ligand_index].pos, /* origin */
			  transformation_matrix);

          /* acceptor position is in the positive X-part of the XY-plane */
          hydrogen_position[X] = 
            bond_length * sin ( ( 180.0 - *angle ) / 180.0 * M_PI );
          hydrogen_position[Y] = 
            bond_length * cos ( ( 180.0 - *angle ) / 180.0 * M_PI );
          hydrogen_position[Z] = 0.0;
          transform_point_back(hydrogen_position, transformation_matrix);      
          *angle = compute_angle(ligand->atoms[ligand_index].pos, 
                                 hydrogen_position, target_atom->pos);
        }else{ 
          i = 0; 
          if(number_of_hydrogens > 0) 
            do{ 
              *angle = compute_angle (ligand->atoms[ligand_index].pos,
		        ligand->atoms[hydrogen_indices[i]].pos,
		        target_atom->pos );		  
              i++; 
            }while(i < number_of_hydrogens && (*angle < MIN_HYDROGEN_ANGLE || 
                                               *angle > MAX_HYDROGEN_ANGLE));
#ifdef TRACE 
        else			  
          printf("! warning: In %s, no H is connected to  ligand atom %d.\n", 
                 global->ligand->name, ligand_index );
#endif	 
      } 
      /* toneroma 06JUN06 - end of code for N orbitals*/

      /* there are added hydrogens, the position of these has to be computed */
      /* toneroma 06JUN06 - changed this code from number_of_added_hydrogens to
       * number_of_hydrogens to handle ALL hydrogens of a ligand for position 
       * computation (except for when N and orbital of ar, am, pl3, or sp2) */
      }else if( number_of_hydrogens > 0){
        *angle = compute_ligand_h_angle(ligand, ligand_index, target_atom->pos,
				       hydrogen_indices[0], global); 
        hydrogen_position[X] = ligand->atoms[hydrogen_indices[0]].pos[X];
        hydrogen_position[Y] = ligand->atoms[hydrogen_indices[0]].pos[Y];
        hydrogen_position[Z] = ligand->atoms[hydrogen_indices[0]].pos[Z];
      }	

      /* HETATM is acting as donor, just allowing it to pass */
      if( target_atom->res == HETATM ){
#ifdef TRACE
        printf("HETATM is donor\n");
#endif
        pre_angle = 160.0;
        preacceptor_pass = TRUE; 
      }else{
        rv = target_preacc_angle(target_atom, hydrogen_position, 
                                 &(global->target_atoms[global->target_residues[residue_index].start_atom]),
                                 global->target_residues[residue_index].number_of_atoms,
                                 &pre_angle);
        if(rv == SUCCESS){
          if(pre_angle >= MIN_PREACCEPTOR_ANGLE && 
             pre_angle <= MAX_PREACCEPTOR_ANGLE) preacceptor_pass = TRUE;
        }
      }
    }

#if 0
  /* This is kept here for printing if one desires to see the computed values
   * for hbond distance and angles */
  printf("H: %f %f %f\n", hydrogen_position[0], hydrogen_position[1],
         hydrogen_position[2]);
  printf("%s %s %s | %s %s %s | %f %f %f\n", 
         target_atom->residue, target_atom->residue_num, target_atom->name,
         ligand->atoms[ligand_index].residue, 
         ligand->atoms[ligand_index].residue_num, 
         ligand->atoms[ligand_index].name, 
         dist_fun(target_atom->pos, ligand->atoms[ligand_index].pos), *angle,
         pre_angle);
#endif

  if(( *angle > MIN_HYDROGEN_ANGLE && *angle < MAX_HYDROGEN_ANGLE ) && 
     preacceptor_pass == TRUE){
    return TRUE;
  }else{
    return FALSE;
  }
}

/* Compute the preacceptor angle for an acceptor in a target residue
 *
 * acceptor:  Pointer to the target acceptor atom
 * hydrogen_position:  Pointer to float array of length 3 holding the hydrogen 
 *                     position
 * residue_atoms:  Pointer to the first atom in the residue that contains
 *                 the acceptor atom
 * number_of_residue_atoms:  Number of atoms in the residue
 * angle:  Pointer to a float that will return the value of the preacceptor 
 *         angle 
 * Returns FAILURE if atom pointed to by the acceptor variable is not a
 *         known acceptor atom -- SUCCESS if the preacceptor angle was computed
 *         or the preacceptor atom could not be found
 */
int
target_preacc_angle(atom_pt acceptor, float *hydrogen_position,
                    atom_pt residue_atoms, int number_of_residue_atoms,
                    float *angle)
{
  int atom_types[3];
  float points[3][3];
  int found_nbr = 0;
  atom_types[0] = 0;
 
  if(acceptor->type == O) atom_types[0] = C; 
  else if(acceptor->type == OXT) atom_types[0] = C;
  else if(acceptor->res == ASN){
    if(acceptor->type == OD1) atom_types[0] = CG;
  }else if(acceptor->res == ASP){
    if(acceptor->type == OD1 || acceptor->type == OD2) atom_types[0] = CG;
  }else if(acceptor->res == GLN){
    if(acceptor->type == OE1) atom_types[0] = CD;
  }else if(acceptor->res == GLU){
    if(acceptor->type == OE1 || acceptor->type == OE2) atom_types[0] = CD;
  /* His not an acceptor in scoring at present */
  }else if(acceptor->res == HIS){
    if(acceptor->type == ND1){
      atom_types[0] = CG;
      atom_types[1] = CE1;  /* Not used presently */
    }else if(acceptor->type == NE2){
      atom_types[0] = CD2;
      atom_types[1] = CE1;  /* Not used presently */
    }
  }else if(acceptor->res == SER){
    if(acceptor->type == OG) atom_types[0] = CB;
  }else if(acceptor->res == THR){
    if(acceptor->type == OG1) atom_types[0] = CB;
  }else if(acceptor->res == TYR){
    if(acceptor->type == OH) atom_types[0] = CZ;
  /* CYS && MET were moved to the end of the if statements since they are
   * unlikely to be used anytime soon  -- no Sulfur hbonds in HBIND at present
   */
  }else if(acceptor->res == CYS){
    if(acceptor->type == SG) atom_types[0] = CB;
  }else if(acceptor->res == MET){
    if(acceptor->type == SD){
      atom_types[0] = CG;
      atom_types[1] = CE;
    }
  }
  /* Not a known protein acceptor atom */
  if(atom_types[0] == 0){
    fprintf(stderr, "(%c)%s %s %s is not a known acceptor atom (hbond_check.c)",
            acceptor->chain_id, acceptor->residue, acceptor->residue_num, 
            acceptor->name);
    return FAILURE;
  }

  found_nbr = get_positions(residue_atoms, number_of_residue_atoms, atom_types,
                            1, points);
  if(!found_nbr){
    fprintf(stderr, "Unable to get the preacceptor atom for "
            "(%c)%s %s %s", acceptor->chain_id, acceptor->residue,
            acceptor->residue_num, acceptor->name);
    fprintf(stderr, "\tIgnoring the preacceptor angle for this hbond\n");
    *angle = 160.0;
  }else
    *angle = compute_angle(points[0], acceptor->pos, hydrogen_position);

  return SUCCESS;
}


int  compute_ligand_hydrogen_angle ( global_data_pt global,
				     int            ligand_donor_index,
				     float          *acceptor_position,
				     float          *angle )
{
  transform_matrix_t transformation_matrix;
  molecule_pt        ligand;
  float              xpoint[3],
    ypoint[3],
    hydrogen_position[3];
  float              bond_length;
  int                hydrogen_indices[3];
  int                *neighbors;
  int                count,
    number_of_hydrogens,
    number_of_added_hydrogens;
  int                i, j;


  ligand = global->ligand;
  for ( i = 0; i < 3; i ++ )
    {
      hydrogen_indices[i] = -1;
      hydrogen_position[i] = 0.0;
    }

  /* first get all hydrogens and check, if they already have a
	 position, or if they were added */
  number_of_hydrogens = 0;
  number_of_added_hydrogens = 0;
  neighbors = ligand->neighbors[ligand_donor_index];
#ifdef TRACE
  printf("----- ligand_donor_index=%d (%7.3f, %7.3f, %7.3f) \n", 
	 ligand_donor_index,
	 ligand->atoms[ligand_donor_index].pos[X],
	 ligand->atoms[ligand_donor_index].pos[Y],
	 ligand->atoms[ligand_donor_index].pos[Z] );
#endif
  for ( i = 0; i < ligand->number_of_neighbors[ligand_donor_index]; i++ )
    {
#ifdef TRACE
      printf("       neighbors[%d]=%d type %s", 
	     i, neighbors[i], ligand->atoms[neighbors[i]].type_str );
#endif
      if ( ligand->atoms[neighbors[i]].type == H )
	{		      
	  hydrogen_indices[number_of_hydrogens++] = neighbors[i];
#ifdef TRACE
	  printf( "-- found neighbors[%d] = hydrogen, hydrogen_indices[%d] = %d\n", 
		  i, number_of_hydrogens-1, hydrogen_indices[number_of_hydrogens-1] );
	  printf( "      number_of_hydrogens = %d\n", number_of_hydrogens );
#endif
	  if ( strcmp ( ligand->atoms[neighbors[i]].type_str,
			"H_ADD" ) == 0 )
	    number_of_added_hydrogens++;
	}
#ifdef TRACE
      else
	printf( "-- found neighbors[%d] != hydrogen \n", i );
#endif
    }
  if ( number_of_added_hydrogens > 0 )
    /* there are added hydrogens, the position of these has to be
	   computed */
    {
      *angle = 
	compute_ligand_h_angle ( ligand,
				 ligand_donor_index,
				 acceptor_position,
				 hydrogen_indices[0],
				 global);
      sprintf ( ligand->atoms[hydrogen_indices[0]].name,
		"H%d",
		ligand_donor_index + 1 );
      sprintf ( ligand->atoms[hydrogen_indices[0]].type_str,
		"H" );
    }	
  else
    /* all hydrogens already have positions, so check, if there is
	   one in a correct position to form a H-bond */
    {
      /* check if the donor is connected via a rotatable bond */
      neighbors = ligand->neighbors[ligand_donor_index];
      count = 0;
      for ( i = 0; i < 3; i++ )
	{
	  ypoint[i] = 0.0;
	  xpoint[i] = ligand->atoms[neighbors[0]].pos[i];
	}
      for ( i = 0; i < ligand->number_of_neighbors[ligand_donor_index]; i++ )
	if ( ligand->atoms[neighbors[i]].type != H )
	  {
	    for ( j = 0; j < 3; j++ )
	      ypoint[j] += ligand->atoms[neighbors[i]].pos[j];
	    count++;
	  }
	else
	  for ( j = 0; j < 3; j++ )
	    hydrogen_position[j] = ligand->atoms[neighbors[i]].pos[j];
      if ( count == 1 )
	{
	  /* only one neighbor, so we just assume that the hydrogen 
	     in pointing towards the acceptor */
	  for ( i = 0; i < 3; i++ )
	    xpoint[i] = acceptor_position[i];
	  bond_length = dist_fun ( ypoint, 
				   ligand->atoms[ligand_donor_index].pos );
	  *angle = compute_angle ( ypoint,
				   ligand->atoms[ligand_donor_index].pos,
				   hydrogen_position );
	  compute_transformation_matrix ( xpoint,              /* X-axis */
					  ypoint,              /* Y-axis */
					  ligand->atoms[ligand_donor_index].pos,
					  /* origin */
					  transformation_matrix );
	  /* acceptor position is in the positive X-part of the XY-plane */
	  hydrogen_position[X] = 
	    bond_length * sin ( ( 180.0 - *angle ) / 180.0 * M_PI );
	  hydrogen_position[Y] = 
	    bond_length * cos ( ( 180.0 - *angle ) / 180.0 * M_PI );
	  hydrogen_position[Z] = 0.0;
	  transform_point_back ( hydrogen_position, 
				 transformation_matrix );      
	  *angle = 
	    compute_angle ( ligand->atoms[ligand_donor_index].pos,
			    hydrogen_position,
			    acceptor_position );
	}
      else
	{
	  i = 0;
	  if ( number_of_hydrogens > 0 )
	    do 
	      {
		*angle = 
		  compute_angle ( ligand->atoms[ligand_donor_index].pos,
				  ligand->atoms[hydrogen_indices[i]].pos,
				  acceptor_position );		  
		i++;
	      }
	    while ( i < number_of_hydrogens
		    && ( *angle < MIN_HYDROGEN_ANGLE 
			 || *angle > MAX_HYDROGEN_ANGLE ) );
#ifdef TRACE
	  else			  
	    printf("! warning: In %s, no H is connected to  ligand atom %d.\n", 
		   global->ligand->name, ligand_donor_index );
#endif	
	}
    }
  if ( *angle > MIN_HYDROGEN_ANGLE && *angle < MAX_HYDROGEN_ANGLE )
    return TRUE;
  else
    return FALSE;
}
