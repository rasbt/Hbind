#include <math.h>
#include "defs.h"
#include "types.h"

/* This routine computes the matrix for the transformation of points
   after rotations around the X-, Z-, and Y-axis in *exact* this order.
   This matrix can be used to transform coordinates based on:
     new(x,y,z) = Rot ( Y ( Rot ( Z ( Rot X ( Trans (x,y,z) ) ) ) ) )
   with all rotation angles have to consider the previous rotations.
*/
void  compute_yzx_matrix ( double             angle[3],
			   double             trans[3],
			   transform_matrix_t matrix )
{
  double  c[3];
  double  s[3];
  int     i;

  for ( i = 0; i < 3; i++ )
    {
      s[i] = sin ( angle[i] );
      c[i] = cos ( angle[i] );
      matrix[i][3] = trans[i];
    }
  matrix[0][0] = c[Y] * c[Z];
  matrix[0][1] = s[X] * s[Y] - c[Y] * s[Z] * c[X];
  matrix[0][2] = s[X] * s[Z] * c[Y] + c[X] * s[Y];

  matrix[1][0] = s[Z];
  matrix[1][1] = c[X] * c[Z];
  matrix[1][2] = (-1.0) * s[X] * c[Z];
  
  matrix[2][0] = (-1.0) * s[Y] * c[Z];
  matrix[2][1] = c[X] * s[Y] * s[Z] + s[X] * c[Y];
  matrix[2][2] = c[X] * c[Y] - s[X] * s[Y] * s[Z];
}

void  compute_transformation_matrix ( float              n[3],
				      float              ca[3],
				      float              cb[3],
				      transform_matrix_t matrix )
{
  double axis[3],
         angle[3],
         trans[3],
         help[3];
  int    i;
  
  for ( i = 0; i < 3; i++ ) axis[i] = (double) cb[i] - (double)  ca[i];

  /* take care that all coordinates are relative of the position of
     CB, thus this is the origin of the new system */
  for ( i = 0; i < 3; i++ ) trans[i] = (double) cb[i];


  if ( axis[Y] == 0.0 ) angle[X] = M_PI/2;  
  /* angle should not be set to 0.0, as was previous.  
  Did not want to pass axix[Y]==0.0 into the next function in order to avoid a divide by 0.  
  atan of (a number / 0) is PI/2.*/
  else{
    angle[X] = (-1) * atan ( axis[Z] / axis[Y] );
    if ( axis[Y] < 0.0 ) angle[X] -= M_PI;
  }

  /* then rotate CA-CB-axis to the Y-axis, here it is important to
     use the 'new' (after X-rotation) length of axis[Y] */
  if ( axis[Y] * cos ( angle[X] ) - axis[Z] * sin ( angle[X] ) == 0.0 )
    angle[Z] = M_PI/2;
  else {
    angle[Z] = atan ( axis[X] / ( axis[Y] * cos ( angle[X] ) 
                      - axis[Z] * sin ( angle[X] ) ) );
    if ( axis[Y] * cos ( angle[X] ) - axis[Z] * sin ( angle[X] ) < 0.0 ) 
      angle[Z] += M_PI;
  }
  

  /* Last a rotation around the Y-axis to align the N to the X-axis,
     help defines the original position of N, for the computation
     of the angle its position after the two previou rotations
     has to be computed. Note that here the transformation has to
     be considered when computing the new position of the atom after
     the previous two rotations. This has not been necessary for
     the axis when computing angle[Z], since the head of axis is the
     position of CB, which is the origin of the new coordinate system. */
  for ( i = 0; i < 3; i++ ) help[i] = n[i] - trans[i];
  if ( help[X] * cos ( angle[Z] )
       - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
       + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) == 0 )
    angle[Y] = M_PI/2;  
  else
    {
      angle[Y] = atan ( ( help[Y] * sin ( angle[X] )
			  + help[Z] * cos ( angle[X] ) )
			/ ( help[X] * cos ( angle[Z] )
			  - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
			  + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) ) );
      if ( help[X] * cos ( angle[Z] )
	   - help[Y] * cos ( angle[X] ) * sin ( angle[Z] )
	   + help[Z] * sin ( angle[X] ) * sin ( angle[Z] ) < 0.0 )
	angle[Y] += M_PI;
    }

  /* now compute the matrix which transforms the current residue in a way
     that CA and CB are aligned to the Y-axis, N is aligned to the X-axis,
     and CB is in the origin of the coordinate-system */
  compute_yzx_matrix ( angle, trans, matrix );
}



float cut (float f)
{
    int temp;
    temp = (int) ( (f + 0.00005) * 10000 );
    return ( (float) temp / 10000 );
}

void  transform_point ( float              pos[3],
			transform_matrix_t matrix )
{
  int   i, j;
  double x[3];


  for ( i = 0; i < 3; i++ )
    x[i] = (double) pos[i] - matrix[i][3];
  for ( i = 0; i < 3; i++ )
    {
      pos[i] = 0.0;
      for ( j = 0; j < 3; j++ )
	pos[i] += (float) ( matrix[i][j] * x[j] );
      cut( pos[i] );
    }
 
  /* PDB-coordinates are only given with three digits accuracy */
  for ( i = 0; i < 3; i++ )
    if ( pos[i] > -0.0005 && pos[i] < 0.0005 )
      pos[i] = 0.0;
}

void  transform_point_back ( float              pos[3],
			     transform_matrix_t matrix )
{
  double x[3];
  int   i, j;


  for ( i = 0; i < 3; i++ )
    x[i] = (double) pos[i];
  for ( i = 0; i < 3; i++ )
    {
      pos[i] = 0.0;
      for ( j = 0; j < 3; j++ )
	pos[i] += (float) ( matrix[j][i] * x[j] );
      cut( pos[i] );
    }
  /* PDB-coordinates are only given with three digits accuracy */
  for ( i = 0; i < 3; i++ )
    {
      pos[i] += (float) matrix[i][3];
      if ( pos[i] > -0.0005 && pos[i] < 0.0005 )
	pos[i] = 0.0;
    }

}

/* This routine takes all atoms with level higher or equal than 'level' and
   rotates them around the Y-axis.
   */
void  rotate_around_y_axis ( atom_pt atoms,
			     int     number,
			     int     level,
			     double  angle )
{
  double  help[3];
  int     i, j;

  for ( i = 0; i < number; i++ )
    if ( atoms[i].level >= level )
      {
	for ( j = 0; j < 3; j++ )
	  help[j] = (double) atoms[i].pos[j];
	atoms[i].pos[X] = (float) ( help[X] * cos ( angle ) 
				    + help[Z] * sin ( angle ) );
	if ( atoms[i].pos[X] > -0.0001 && atoms[i].pos[X] < 0.0001 )
	  atoms[i].pos[X] = 0.0;
	atoms[i].pos[Z] = (float) ( help[Z] * cos ( angle ) 
				    - help[X] * sin ( angle ) );
	if ( atoms[i].pos[Z] > -0.0001 && atoms[i].pos[Z] < 0.0001 )
	  atoms[i].pos[Z] = 0.0;
      }
}

/* this routine rotates a single atom around the Y-axis */
void  rotate_single_atom_around_y_axis ( float   *position,
					 double  angle )
{
  double  help[3];
  int     i;

  for ( i = 0; i < 3; i++ )
    help[i] = (double) position[i];
  position[X] = (float) ( help[X] * cos ( angle ) 
			  + help[Z] * sin ( angle ) );
  if ( position[X] > -0.0001 && position[X] < 0.0001 )
    position[X] = 0.0;
  position[Z] = (float) ( help[Z] * cos ( angle ) 
			  - help[X] * sin ( angle ) );
  if ( position[Z] > -0.0001 && position[Z] < 0.0001 )
    position[Z] = 0.0;
}

