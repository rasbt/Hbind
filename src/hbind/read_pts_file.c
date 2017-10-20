#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"

int  read_pts_file ( global_data_pt  global,
		      char            *filename )
{
  FILE            *in;
  pts_compound_pt compounds;
  char            linebuffer[128],
                  name[MAX_COMPOUNDNAME_LEN],
                  oldname[MAX_COMPOUNDNAME_LEN];
  float           *pos;
  int             *act,
                  *atom_index;
  char            type_char;
  int             number_of_compounds,
                  points,
                  count,
                  ind;
  char            err_msg[FILENAME_MAX];

  struct stat file_stat_buf;
  if(stat(filename, &file_stat_buf)){
    sprintf(err_msg, "FATAL ERROR:  Unable to open the file %s\n%s\n\n", 
            filename, strerror(errno));
    err_print(err_msg);
    fprintf(stderr, err_msg);
    printf("FATAL ERROR: Unable to open the interaction points file %s -- Please see err file\n", filename);
    return FAILURE;
  }else if(file_stat_buf.st_size == 0){
    sprintf(err_msg, "\nFATAL ERROR: Database file %s is empty\n\n", filename);
    err_print(err_msg);
    fprintf(stderr, err_msg);
    fprintf(stdout, err_msg);
    return FAILURE;
  }
  in = fopen ( filename, "r" );
  if( in == NULL){
    perror ( "ERROR  " );
    sprintf (err_msg, "unable to open file %s", filename );
    err_warning2 ("read_pts_file", err_msg);
    return FAILURE;
  }

  compounds = global->compound_interactions->compounds;
  act = global->compound_interactions->act;
  pos = global->compound_interactions->pos;
  atom_index = global->compound_interactions->atom_index;
  count = 0; 
  points = 0;
  number_of_compounds = 0;

  while ( fgets ( linebuffer, sizeof (linebuffer), in ) != NULL ){
    if(linebuffer[0] == '\0' || linebuffer[0] == '\n' || linebuffer[0] == ' ' )
      continue;
    ind = 3 * count;
    sscanf(linebuffer, "%s %c %f %f %f %d", name, &type_char, &pos[ind+X],
           &pos[ind+Y], &pos[ind+Z], &atom_index[count] );
    if ( type_char == 'A' ) act[count] = ACCEPTOR;
    else if ( type_char == 'D' ) act[count] = DONOR;
    else if ( type_char == 'H' ) act[count] = HYDROPHOB;
    else if ( type_char == 'N' ) act[count] = DONEPTOR;
    else {
      printf ( "%s", linebuffer );
      err_panic2 ( "read_pts_file", "illegal interaction type");
    }
    if ( count == 0 ){
      strcpy ( oldname, name ); 	
      strcpy ( compounds[0].name, name ); 
      compounds[0].first_point = 0;
    }
    if ( points > 0 && strcmp ( name, oldname ) != 0 ) {
      compounds[number_of_compounds].number_of_points = points;
      number_of_compounds++;
      strcpy ( compounds[number_of_compounds].name, name ); 
      strcpy ( oldname, name ); 
      compounds[number_of_compounds].first_point = count;
      points = 0;
    }
    count++;
    points++;
         
    if ( count > MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS )
      err_panic2("read_pts_file", 
                 "more than MAX_NUMBER_OF_TOTAL_COMPOUND_INTERACTION_POINTS "
                 " interaction_points");
    if ( number_of_compounds+1 > MAX_NUMBER_OF_PTS_COMPOUNDS )
      err_panic2 ( "read_pts_file",
		  "more than MAX_NUMBER_OF_PTS_COMPOUNDS compounds");
  }
  compounds[number_of_compounds].number_of_points = points;
  global->compound_interactions->number_of_compounds = number_of_compounds + 1;
  global->compound_interactions->number_of_interaction_points = count;
  fclose ( in );
  if ( count == 0 ) return FAILURE;
  return SUCCESS;
}
  
