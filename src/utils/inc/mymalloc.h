#ifndef _MYMALLOC_H
#define _MYMALLOC_H

#include <stdlib.h>
#include <err_handle.h>
#include <mymalloc.h>

static inline void*
mymalloc ( size_t  size )
{
  void  *adr;
  adr = malloc ( size );
  if ( !adr ) err_panic( "mymalloc", "malloc failed" );
  return ( adr );
}

static inline void*
myrealloc ( void    *old, size_t  size )
{
  void  *adr;

  adr = realloc ( old, size );
  if ( !adr ) err_panic( "myrealloc", "realloc failed" );
  return ( adr );
}

static inline void
my_free(void *ptr)
{
  if(ptr) free(ptr);
  ptr = NULL;
}
#endif
