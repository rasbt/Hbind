#include <stdio.h>
#include "read_template.h"
#include "err_handle.h"

/*
 *  Function 'read_template()' reads the coordinates and the types
 *  of the interaction-points from the template-file 'filename'.
 *  It returns the number of points read. If the file contains more
 *  than 'MAX_TEMPLATE_POINTS' then the function is exited via 'err_panic()'
 */
int read_template(char           *filename, global_data_pt global )
{
  FILE           *in;
  char           act[3];
  char           err_msg[FILENAME_MAX];
  char           line[MAX_TEMPLATE_LINELENGTH];
  interaction_pt inter, old_inter;
  int            number_of_points, /* this is the internal number of points */
                 number_of_read_points, /* count doneptors as single points */
                 number_of_key_points,
                 number_of_key2_points,
                 i;
  
  in = fopen ( filename, "r" );
  if ( in == NULL ) {
    perror ( "ERROR  " );
    err_panic2 ( "read_template", "unable to open file");
  }
  
  number_of_points = 0;
  number_of_key_points = 0;
  number_of_key2_points = 0;
  number_of_read_points = 0;
  inter = global->template_interactions;

  /* read interactions data for all interaction-points in the template */
  printf ("\nTemplate Parameters:\n");
  while ( fgets ( line, MAX_TEMPLATE_LINELENGTH, in ) != NULL ){
    if ( line[0] == '#' ){
      printf ("%s", line);
    }
    /* comment line or empty line, skip */
    if ( line[0] == '#' || line[0] == '\0' || line[0] == '\n' ) continue;
      
    /* read activity and position of this template point */
    sscanf(line, "%s %f %f %f", act, &inter->pos[X], &inter->pos[Y],
           &inter->pos[Z] );
    if ( act[1] == '*' ) {
      inter->key_point = TRUE;
      if ( act[2] == '*' ) inter->key_point = COVALENT;
      else number_of_key_points++;
      if (act[2] == '^' )
	{
	  inter->key2_point = TRUE;
	  number_of_key2_points++;
	}
    }
    else if ( act[1] == '^')
      {
	inter->key2_point = TRUE;
	number_of_key2_points++;
	if (act[2] == '*' )
	  {
	    inter->key_point = TRUE;
	    number_of_key_points++;
	  }
	
      }
    else 
      {
	inter->key_point = FALSE;
	inter->key2_point = FALSE;
      }

    /* determine activity of this point */
    inter->act = UNKNOWN;
    if ( act[0] == 'A' ) {
      inter->act = ACCEPTOR;
      inter->ref = number_of_read_points;  
    } if ( act[0] == 'D' ) {
      inter->act = DONOR;
      inter->ref = number_of_read_points; 
    } if ( act[0] == 'H' ) {
      inter->act = HYDROPHOB;
      inter->ref = number_of_read_points; 
    } if ( act[0] == 'N' ){
    /* this is a doneptor point, i.e., a template point where either
       a donor or an acceptor, or a doneptor of a ligand can be placed,
       this point is internally represented as two separate points,
       a donor and an acceptor point, which share the same position */
      inter->act = ACCEPTOR;
      inter->ref = number_of_read_points; 
      inter->type = UNKNOWN;
      inter->index = UNKNOWN;
      number_of_points++;
      old_inter = inter;
      inter++;
      inter->act = DONOR;
      inter->ref = number_of_read_points; 
      inter->pos[X] = old_inter->pos[X];
      inter->pos[Y] = old_inter->pos[Y];
      inter->pos[Z] = old_inter->pos[Z];
      inter->key_point = old_inter->key_point;
      inter->key2_point = old_inter->key2_point;
    }
    if ( inter->act == UNKNOWN ) {
      sprintf ( err_msg, "Wrong activity definition: %s", act );
      err_panic2 ( "read_template", err_msg);
    }
    inter->type = UNKNOWN;
    inter->index = UNKNOWN;

    inter++;
    number_of_points++;
    number_of_read_points++;
    if ( number_of_points >= MAX_TEMPLATE_POINTS ) {
      /* keep in mind that the internal number of template points is larger
         than the number of points in the template file, since we represent
         each doneptor point by two points */
      sprintf(err_msg, "Number of (A + D + 2*N + H) template points more than "
              "MAX_TEMPLATE_POINTS (%d)", MAX_TEMPLATE_POINTS);
      err_panic2 ( "read_template", err_msg);
    }
  }

  fclose ( in );
  global->template_size = number_of_points;
  if ((number_of_key_points == 0) && (number_of_key2_points == 0))
    err_panic2 ( "read_template", "no template key points");
  if ((global->match_2_key_points) && (number_of_key2_points == 0))
    err_panic2 ( "read_template", "no secondary key points");
  global->template_key_points = number_of_key_points;
  global->template_key2_points = number_of_key2_points;

  printf ("\nCHECK: template file read\n");
  for (i = 0; i < 3; i++)
    printf ("   %d: %s %8.3f %8.3f %8.3f\n",
	    i, 
	    global->template_interactions[i].act == DONOR       ? "DONOR    " 
	    : global->template_interactions[i].act == ACCEPTOR  ? "ACCEPTOR "  
	    : global->template_interactions[i].act == HYDROPHOB ? "HYDROPHOB" 
	    : "NOTHING  ",
	    global->template_interactions[i].pos[X],
	    global->template_interactions[i].pos[Y],
	    global->template_interactions[i].pos[Z] );
  
  return number_of_read_points;
}
       
