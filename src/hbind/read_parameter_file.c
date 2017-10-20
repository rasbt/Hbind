/*
 * $Source: /psa/share/repository/hbind/src/hbind/read_parameter_file.c,v $
 * $Revision: 1.13 $
 * $Author: toneroma $
 * $Date: 2009/06/29 18:57:24 $
 *
 * $Log: read_parameter_file.c,v $
 * Revision 1.13  2009/06/29 18:57:24  toneroma
 * added the ability to track molecule & conformer correctly so that the group_conformers will work properly (and so naming is consistent no matter the options set)
 *
 * Revision 1.12  2009/02/26 21:04:15  vanvoor4
 * no need to send global to err functions
 *
 * Revision 1.11  2008/09/08 14:21:02  vanvoor4
 * Removed string issue when reading "None" for the restart molecule
 * option
 *
 * Revision 1.10  2008/09/02 15:15:06  toneroma
 * changes to allow restarting of runs
 *
 * Revision 1.9  2007/10/09 21:33:31  toneroma
 * fixed problem with match_2_key_points and group_conformers where hbind.parameters in the /in/ directory wouldn't be used
 *
 * Revision 1.8  2007/09/28 18:33:48  toneroma
 * *** empty log message ***
 *
 * Revision 1.7  2007/02/08 19:51:33  toneroma
 * added error messages for using the incorrect hbind.parameters file
 *
 * Revision 1.6  2007/01/29 17:41:51  toneroma
 * Changed formatting of parameters being printed to output file
 *
 * Revision 1.5  2007/01/17 20:42:33  vanvoor4
 * Moved the more delicate parameters to inc/defs.H from hbind.parameters.
 *
 * Revision 1.4  2006/09/20 14:27:41  vanvoor4
 * Need to be really careful when using the *scanf family of functions or
 * any other method of scanning that depends on whitespace for identifying
 * fields.  For example in the hbind.parameters file the unit field may not
 * have any embedded spaces--that is (A) or (_A_) or _MY_UNITS_ are acceptable
 * where as ( A) or (A ) or ( A ) are not.
 *
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"
#include "basics.h"

#define MAX_PARAMETER_VARIABLE_LEN 32
#define NUM_PARAMETER_TYPES 16 

int  scan_parameter_file_line ( char *line,
				char *value,
				global_data_pt global)
{
  char variable[MAX_PARAMETER_LINELENGTH];
  char old_value[MAX_PARAMETER_LINELENGTH];
  char unit[MAX_PARAMETER_LINELENGTH];
  char junk[MAX_PARAMETER_LINELENGTH]; /* to catch bad lines */
  char *help;
  int len, i, rv;
  char err_msg[FILENAME_MAX];


  static char param_labels[NUM_PARAMETER_TYPES][MAX_PARAMETER_VARIABLE_LEN + 1];
  static int initialized = 0;
  if(!initialized){
    strcpy(param_labels[0], "FLEX_BOND_DEFN_FILE");
    strcpy(param_labels[1], "DATA_ROOT");
    strcpy(param_labels[2], "DME_THRESHOLD");
    strcpy(param_labels[3], "RMS_THRESHOLD");
    strcpy(param_labels[4], "ANCHOR_TRANSLATION");
    strcpy(param_labels[5], "ANCHOR_OVERLAP");
    strcpy(param_labels[6], "SIDE_CHAIN_OVERLAP");
    strcpy(param_labels[7], "INTRA_OVERLAP");
    strcpy(param_labels[8], "INTERMEDIATELY_TOLERATED_OVERLAP");
    strcpy(param_labels[9], "FINALLY_TOLERATED_OVERLAP");
    strcpy(param_labels[10], "FINALLY_TOLERATED_MAX_BUMP");
    strcpy(param_labels[11], "SCORE_CUTOFF");
    strcpy(param_labels[12], "MAX_TEMPLATE_TRIANGLES");
    strcpy(param_labels[13], "MATCH_2_KEY_POINTS");
    strcpy(param_labels[14], "GROUP_CONFORMERS");
    strcpy(param_labels[15], "RESTART_MOLECULE");
    initialized = 1;
  }

  if ( line[0] == '#' || line[0] == ' ' || line[0] == '\0' || line[0] == '\n' )
    return IGNORE;
  rv = sscanf(line, "%s %s %s %s %*s", variable, unit, value, junk);
  if(rv == EOF)
    {
      sprintf(err_msg, "Error when scanning the line:\n\t%sin the "
	      "hbind.parameters file: %s"
	      "Make sure the parameters file matches the format of the default parameters file in:\n"
	      "%s/params/hbind.parameters\n\n", line, strerror(errno), getenv("HBIND_DIR"));
      err_print(err_msg);
      fprintf(stderr, err_msg);
      fprintf(stdout, err_msg);
      exit(1);
    }
  else if(rv > 3)
    {
      sprintf(err_msg, "Error: The line\n\t%sin the hbinds.parameter file has "
	      "too many fields.  Each line must have exactly\n 3 fields\n"
	      "Make sure the parameters file matches the format of the default parameters file in:\n"
	      "%s/params/hbind.parameters\n\n", line, getenv("HBIND_DIR"));
      err_print(err_msg);
      fprintf(stderr, err_msg);
      printf(err_msg);
      exit(1);
    }
  else if(rv < 3)
    {
      sprintf(err_msg, "Error: The line\n\t%sin the hbinds.parameter file has "
	      "too few fields.  Each line must have exactly\n 3 fields\n"
	      "Make sure the parameters file matches the format of the default parameters file in:\n"
	      "%s/params/hbind.parameters\n\n", line, getenv("HBIND_DIR"));
      err_print(err_msg);
      fprintf(stderr, err_msg);
      printf(err_msg);
      exit(1);
    }
  else if(rv == 3){
    /* this token starts with an environment variable */
    if ( *value == '$' ) {
      help = value;
      while ( *help != '/' && *help != '\0' ) help++;
      *help = '\0';
      help++;
      len = strlen ( line );
      strncpy ( old_value, help, len - ( help - line ) );
      sprintf ( value, "%s/%s", getenv ( value + 1 ), old_value );
    }
    
    for(i = 0; i < NUM_PARAMETER_TYPES; i++)
      if(strcmp(param_labels[i], variable) == 0) return i + 1;
  }

  return INVALID_ENTRY_OR_ERROR;
}

void  check_environment_parameter_variables ( global_data_pt  global )
{
  char   *value;

  value = getenv ( "DME_THRESHOLD" );
  if ( value != NULL )
    {
      printf ( "environment variable: DME_THRESHOLD = %s\n", value );
      global->dmetol = atof ( value );
    }
  value = getenv ( "RMS_THRESHOLD" );
  if ( value != NULL ) 
    {
      printf ( "environment variable: RMS_THRESHOLD = %s\n", value );
      global->rmstol = atof ( value );
    }
  value = getenv ( "ANCHOR_TRANSLATION" );
  if ( value != NULL )
    {
      printf ( "environment variable: ANCHOR_TRANSLATION = %s\n", value );
      global->anchor_translation = atof ( value );
    }
  value = getenv ( "ANCHOR_OVERLAP" );
  if ( value != NULL )
    {
      printf ( "environment variable: ANCHOR_OVERLAP = %s\n", value );
      global->anchor_overlap = atof ( value );
    }
  value = getenv ( "SIDE_CHAIN_OVERLAP" );
  if ( value != NULL )
    {
      printf ( "environment variable: SIDE_CHAIN_OVERLAP = %s\n", value );
      global->side_chain_overlap = atof ( value );
    }
  value = getenv ( "INTRA_OVERLAP" );
  if ( value != NULL )
    {
      printf ( "environment variable: INTRA_OVERLAP = %s\n", value );
      global->intra_overlap = atof ( value );
    }
  value = getenv ( "INTERMEDIATELY_TOLERATED_OVERLAP" );
  if ( value != NULL )
    {
      printf ( "environment variable: INTERMEDIATELY_TOLERATED_OVERLAP = %s\n",
	       value );
      global->intermediate_overlap = atof ( value );
    }
  value = getenv ( "FINALLY_TOLERATED_OVERLAP" );
  if ( value != NULL )
    {
      printf("environment variable: FINALLY_TOLERATED_OVERLAP = %s\n", value);
      global->finally_tolerated_overlap = atof ( value );
    }
  /***** Added by PCS -- 14-May-01 *****/
  value = getenv ( "FINALLY_TOLERATED_MAX_BUMP" );
  if ( value != NULL )
    {
      printf("environment variable: FINALLY_TOLERATED_MAX_BUMP = %s\n", value );
      global->finally_tolerated_max_bump = atof ( value );
    }
  /*************************************/
  value = getenv ( "SCORE_CUTOFF" );
  if ( value != NULL )
    {
      printf ( "environment variable: SCORE_CUTOFF = %s\n", value );
      global->score_cutoff = atof ( value );
    }
  value = getenv ( "MAX_TEMPLATE_TRIANGLES" );
  if ( value != NULL )
    {
      printf ( "environment variable: MAX_TEMPLATE_TRIANGLES = %s\n", value );
      global->max_template_triangles = atoi ( value );
    }
  value = getenv ( "MATCH_2_KEY_POINTS" ); /* NOTE: Environment variables aren't used for this function, and so this code is "broken, and I will fix when I can - toneroma 06MAR08 */
  if ( value != NULL )
    {
      printf ( "environment variable: MATCH_2_KEY_POINTS = %s\n", value );
      global->max_template_triangles = atoi ( value );
    }
  value = getenv ( "GROUP_CONFORMERS" ); /* NOTE: Environment variables aren't used for this function, and so this code is "broken, and I will fix when I can - toneroma 06MAR08 */
  if ( value != NULL )
    {
      printf ( "environment variable: GROUP_CONFORMERS = %s\n", value );
      global->max_template_triangles = atoi ( value );
    }
  value = getenv ( "RESTART_MOLECULE" ); /* NOTE: Environment variables aren't used for this function, and so this code is "broken, and I will fix when I can - toneroma 06MAR08 */
  if ( value != NULL )
    {
      printf ( "environment variable: RESTART_MOLECULE = %s\n", value );
      global->max_template_triangles = atoi ( value );
    }
}

int read_parameter_file(global_data_pt global, char *file )
{
  FILE  *fp;
  char line[MAX_PARAMETER_LINELENGTH];
  char variable[MAX_PARAMETER_LINELENGTH];
  char unit[MAX_PARAMETER_LINELENGTH];
  char value[MAX_PARAMETER_LINELENGTH];    
  int i;
  char err_msg[FILENAME_MAX];

  fp = fopen ( file, "r" );
  if ( fp == NULL ) return FAILURE;

  printf ( "\nreading parameter-file %s\n", file );  
  while ( fgets ( line, MAX_PARAMETER_LINELENGTH, fp ) )
    switch ( scan_parameter_file_line ( line, value, global ) ){
    case DME_THRESHOLD:
      global->dmetol = atof ( value );
      printf ("> DME_THRESHOLD (A):              %4.2f\n", global->dmetol);
      break;
    case RMS_THRESHOLD:
      global->rmstol = atof ( value );
      printf ("> RMS_THRESHOLD (A):              %4.2f\n", global->rmstol);
      break;
    case ANCHOR_TRANSLATION:
      global->anchor_translation = atof ( value );
      printf ("> ANCHOR_TRANSLATION (A):         %4.2f\n",
	      global->anchor_translation);
      break;
    case ANCHOR_OVERLAP:
      global->anchor_overlap = atof ( value );
      printf ("> ANCHOR_OVERLAP (A):             %4.2f\n",
	      global->anchor_overlap);
      break;
    case SIDE_CHAIN_OVERLAP:
      global->side_chain_overlap = atof ( value );
      printf ("> SIDE_CHAIN_OVERLAP (A):         %4.2f\n",
	      global->side_chain_overlap);
      break;
    case INTRA_OVERLAP:
      global->intra_overlap = atof ( value );
      printf ("> INTRA_OVERLAP (A):              %4.2f\n",
	      global->intra_overlap);
      break;
    case INTERMEDIATELY_TOLERATED_OVERLAP:
      global->intermediate_overlap = atof ( value );
      printf ("> INTERMEDIATE_OVERLAP (A):       %4.2f\n",
	      global->intermediate_overlap);
      break;
    case FINALLY_TOLERATED_OVERLAP:
      global->finally_tolerated_overlap = atof ( value );
      printf ("> FINALLY_TOLERATED_OVERLAP (A):  %4.2f\n",
	      global->finally_tolerated_overlap);
      break;
      /***** Added by PCS -- 14-May-01 *****/
    case FINALLY_TOLERATED_MAX_BUMP:
      global->finally_tolerated_max_bump = atof ( value );
      printf ("> FINALLY_TOLERATED_MAX_BUMP (A): %4.2f\n",
	      global->finally_tolerated_max_bump);
      break;
      /*************************************/
    case SCORE_CUTOFF:
      global->score_cutoff = atof ( value );
      printf ("> SCORE_CUTOFF (kcal/mol):        %4.2f\n",
	      global->score_cutoff);
      break;
    case MAX_TEMPLATE_TRIANGLES:
      global-> max_template_triangles = atoi ( value );
      printf ("> MAX_TEMPLATE_TRIANGLES (#):     %d\n",
	      global->max_template_triangles);
      break;
    case MATCH_2_KEY_POINTS:
      for(i = 0; i < strlen(value); i++) value[i] = toupper(value[i]);
      if(!strcmp(value, "TRUE") || !strcmp(value, "YES")){
        global->match_2_key_points = TRUE;
        printf("NOTE: Each match requires 2 key template points.\n");
      }
      else global->match_2_key_points = FALSE;
      break;
#ifndef OUTPUT_ALL_MATCHES
    case GROUP_CONFORMERS:
      for(i = 0; i < strlen(value); i++) value[i] = toupper(value[i]);
      if(!strcmp(value, "TRUE") || !strcmp(value, "YES")){
        global->group_conformers = TRUE;
        printf("NOTE: Only the best conformer of a ligand will be output.\n");
      }
      else {
	global->group_conformers = FALSE;
        printf("NOTE: All conformers of a ligand will be output.\n");
      }
     break;
#endif
    case RESTART_MOLECULE:
      /*      for(i = 0; i < strlen(value); i++) value[i] = toupper(value[i]);*/
      if((!strcmp(value, "NONE")) || (!strcmp(value, "None")) || (!strcmp(value, "none"))) global->restart_molecule_check = FALSE;
      else{
        strcpy(global->restart_molecule, value);
        printf("NOTE: HBIND run restarting at molecule: %s.\n", value);
        global->restart_molecule_check = TRUE;
      }
      break;
    case INVALID_ENTRY_OR_ERROR:
      sscanf(line, "%s %s %s %*s", variable, unit, value);
      sprintf(err_msg, "\n * INVALID HBIND.PARAMETERS LINE * \nparameter: %s\n"
              "unit: %s\nvalue: %s\n"
	    "Make sure the parameters file matches the format of the default parameters file in:\n"
	    "%s/params/hbind.parameters\n\n", variable, unit, value, getenv("HBIND_DIR") );
      err_print(err_msg);
      fprintf(stderr, err_msg);
      printf(err_msg);
      exit(1);
      return FAILURE;
      break;
    case IGNORE:
    default:
      break;
    }
  
  if ( scan_parameter_file_line ( line, value, global) != NUM_PARAMETER_TYPES )
    {
      sprintf(err_msg, "\n * INCORRECT # OF LINES IN hbind.paramters * \n"
	      "Make sure the parameters file matches the format of the default parameters file in:\n"
	      "%s/params/hbind.parameters\n\n", getenv("HBIND_DIR") );
      err_print(err_msg);
      fprintf(stderr, err_msg);
      printf(err_msg);
      exit(1);
    }
  fclose ( fp );
  return SUCCESS;
}

