#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <err_handle.h>

int hbind_strtod(const char *nptr, double *rv)
{
  char *end;
  *rv = strtod(nptr, &end);
  if(end == nptr){
    char msg[200];
    snprintf(msg, 200, "Unable to perform the conversion of \"%s\" to double."
             "\n", nptr);
    err_error("hbind_strtod", msg);
    return 0;
  } 
  return 1;
}

FILE *hbind_fopen(const char *path, const char *mode)
{
  FILE *rv;
  int errsv;
  rv = fopen(path, mode);
  if(rv == 0){
    errsv = errno;
    fprintf(stderr, "Unable to open the file: %s\nReason: %s\n",
            path, strerror(errsv));
  }
  return rv;
}
