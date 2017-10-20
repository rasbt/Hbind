#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <basics.h>
#include <err_handle.h>

static char err_fname[FILENAME_MAX + 1];
static size_t err_fname_len = 0;


void set_err_filename(char *fname)
{
  if(strlen(fname) > FILENAME_MAX)
    err_panic("set_err_filename", "Path to error file is too long"); 
  strcpy(err_fname, fname);
  err_fname_len = strlen(fname);
}

void err_print(char *message)
{ 
  FILE *log_fp = stderr;
  if(err_fname_len) log_fp = hbind_fopen(err_fname,"a+");
  if(log_fp){
    fprintf ( log_fp, message );
    if(err_fname_len) fclose(log_fp);
  }
}

void err_panic(char  *function, char  *message)
{
  fprintf(stdout, "FATAL ERROR: %s (%s)\n", message, function );
  fprintf(stderr, "FATAL ERROR: %s (%s)\n", message, function );
  exit(EXIT_FAILURE);
}

void err_panic2(char  *function, char  *message)
{
  FILE *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf(log_fp, "FATAL ERROR: %s (%s)\n", message, function );
      fclose(log_fp);
    }
  }
  fprintf(stderr, "FATAL ERROR: %s (%s)\n", message, function );
  printf ( "FATAL ERROR: %s (%s) -- Please see err file\n", message, function );
  exit ( EXIT_FAILURE);
}

void err_error( char  *function, char  *message )
{
  fprintf ( stderr, "ERROR: %s (%s)\n", message, function );
  printf ( "ERROR: %s (%s) -- Please see err file\n", message, function );
}

void err_error2(char  *function, char  *message)
{
  FILE *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf ( log_fp, "ERROR: %s (%s)\n", message, function );
      fclose(log_fp);
    }
  }
  fprintf ( stderr, "ERROR: %s (%s)\n", message, function );
  printf ( "ERROR: %s (%s) -- Please see err file\n", message, function );
}

void  err_warning ( char  *function, char  *message )
{
  fprintf ( stderr, "WARNING: %s (%s)\n", message, function );
  printf ( "WARNING: %s (%s) -- Please see err file\n", message, function );
}

void  err_warning2 ( char  *function, char  *message)
{
  FILE           *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf ( log_fp, "WARNING: %s (%s)\n", message, function );
      fclose(log_fp);
    }
  }
  fprintf ( stderr, "WARNING: %s (%s)\n", message, function );
  printf ( "WARNING: %s (%s) -- Please see err file\n", message, function );
}

void err_usage(char  *message)
{
  fprintf ( stderr, "USAGE: %s\n", message );
  exit ( EXIT_FAILURE );
}

void  err_unknown_atom ( char  *name, char  *residue, char  *residue_num )
{
  fprintf(stdout, "ERROR: unknown atom type: %s %s %s\n", name, residue,
	  residue_num );
  fprintf(stderr, "ERROR: unknown atom type: %s %s %s\n", name, residue,
	  residue_num );
  /* whenever we see an unknown atom type, residue, or level in the PDB file,
     HBIND exits, since this is something the user should have fixed prior
     to screening */
  exit ( EXIT_FAILURE );
}

void err_unknown_atom2(char  *name, char  *residue, char  *residue_num)
{
  FILE           *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf(log_fp, "ERROR: unknown atom type: %s %s %s\n", name, residue,
              residue_num );
      fclose(log_fp);
    }
  }
  fprintf ( stdout, "ERROR: unknown atom type: %s %s %s\n", name, residue,
	    residue_num );
  fprintf ( stderr, "ERROR: unknown atom type: %s %s %s\n", name, residue,
	    residue_num );

  exit ( EXIT_FAILURE );
}
	    
void  err_unknown_residue ( char  *residue, char  *residue_num )
{
  fprintf(stdout, "ERROR: unknown residue type: %s %s\n", residue, residue_num);
  fprintf(stderr, "ERROR: unknown residue type: %s %s\n", residue, residue_num);
  exit ( EXIT_FAILURE );
}

void  err_unknown_residue2 ( char  *residue, char  *residue_num)
{
  FILE           *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf(log_fp, "ERROR: unknown residue type: %s %s\n", residue, 
              residue_num);
      fclose(log_fp);
    }
  }
  fprintf(stdout, "ERROR: unknown residue type: %s %s\n", residue, residue_num);
  fprintf(stderr, "ERROR: unknown residue type: %s %s\n", residue, residue_num);
  exit ( EXIT_FAILURE );
}

void  err_unknown_level ( char  *name, char  *residue, char  *residue_num )
{
  fprintf ( stdout, "ERROR: unknown atom level: %s %s %s\n", name, residue,
	    residue_num );
  fprintf ( stderr, "ERROR: unknown atom level: %s %s %s\n", name, residue,
	    residue_num );
  exit ( EXIT_FAILURE );
}  

void  err_unknown_level2 ( char  *name, char  *residue, char  *residue_num)
{
  FILE           *log_fp;
  if(err_fname_len){
    log_fp = hbind_fopen(err_fname,"a+");
    if(log_fp){
      fprintf(log_fp, "ERROR: unknown atom level: %s %s %s\n", name, residue,
	      residue_num );
      fclose(log_fp);
    }
  }
  fprintf ( stdout, "ERROR: unknown atom level: %s %s %s\n", name, residue,
	    residue_num );
  fprintf ( stderr, "ERROR: unknown atom level: %s %s %s\n", name, residue,
	    residue_num );
  exit ( EXIT_FAILURE );
}  

