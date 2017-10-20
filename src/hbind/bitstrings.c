#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "types.h"
#include "mymalloc.h"
 
int  min_int ( int  a, 
	       int  b )
{
  if ( a < b )
    return a;
  else 
    return b;
}

int  bitstring_get_bit ( bitstring_pt  vec, 
			 unsigned int  pos )
{
  if ( pos > vec->len )
    return -1;
  else
    return ( ( vec->bits[pos/8] >> ( 7 - pos % 8 ) ) & 1 );
}

int  bitstring_all_clear ( bitstring_pt  vec )
{
  unsigned int  c, len;
  
  len = vec->len;
  for ( c = 0; c < len / 8; c++ ) 
    if ( vec->bits[c] != 0 )
      return FALSE;
  for ( c = len - len % 8; c < len; c++ ) 
    if ( bitstring_get_bit ( vec, c ) != 0 )
      return FALSE;
  return TRUE;
}

int  bitstring_all_set ( bitstring_pt  vec )
{
  unsigned int  c, len;
  
  len = vec->len;
  for ( c = 0; c < len / 8; c++ ) 
    if ( vec->bits[c] != 255 )
      return FALSE;
  for ( c = len - len % 8; c < len; c++ )
    if ( bitstring_get_bit ( vec, c ) != 1 )
      return FALSE;
  return TRUE;
}

void  bitstring_and_all ( bitstring_pt  vec, 
			  bitstring_pt  vec2 )
{
  unsigned int  c, len_c;
  
  len_c = min_int ( vec->len / 8 + 1, vec2->len / 8 + 1 );
  for ( c = 0; c < len_c; c++ )
    vec->bits[c] &= vec2->bits[c];
}

void  bitstring_clear_all ( bitstring_pt  vec )
{
  unsigned int  c;
  
  for ( c = 0; c < vec->len / 8 + 1; c++ )
    vec->bits[c] = 0;
}

int  bitstring_clear_bit ( bitstring_pt  vec, 
			   unsigned int  pos )
{
  if ( bitstring_get_bit ( vec, pos ) == 1 )
    {
      vec->bits[pos/8] -= 1 << ( 7 - pos % 8 );
      return 0;
    }
  else
    return -1;
}

int  bitstring_cmp_all ( bitstring_pt  vec, 
			 bitstring_pt  vec2 )
{
  unsigned int  c, len;
  
  len = min_int ( vec->len, vec2->len );
  for ( c = 0; c < len / 8; c++ )
    if ( vec->bits[c] != vec2->bits[c] )
      return 1;
  for ( c = len - len % 8; c < len; c++ )
    if ( bitstring_get_bit ( vec, c ) != bitstring_get_bit ( vec2, c ) )
      return 1;
  return 0;
}

int  bitstring_cmp_and_0 ( bitstring_pt  vec1, 
			   bitstring_pt  vec2 )
{
  unsigned int  c, len;
  
  len = min_int ( vec1->len, vec2->len );
  for ( c = 0; c < len / 8; c++ )
    if ( vec1->bits[c] & vec2->bits[c] )
      return FALSE;
  for ( c = len - len % 8; c < len; c++ )
    if ( bitstring_get_bit ( vec1, c ) && bitstring_get_bit ( vec2, c ) )
      return FALSE;
  return TRUE;
}

int  bitstring_cmp_or_1 ( bitstring_pt  vec1, 
			  bitstring_pt  vec2 )
{
  unsigned int  c, len;
  
  len = min_int ( vec1->len, vec2->len );
  for ( c = 0; c < len / 8; c++ )
    if ( ( vec1->bits[c] | vec2->bits[c] ) != 255 )
      return FALSE;
  for ( c = len - len % 8; c < len; c++ )
    if ( ! bitstring_get_bit ( vec1, c ) && ! bitstring_get_bit ( vec2, c ) )
      return FALSE;
  return TRUE;
}

int bitstring_cmp_xor_1 ( bitstring_pt  vec1, 
			  bitstring_pt  vec2 )
{
  unsigned int  c, len;
  
  len = min_int ( vec1->len, vec2->len );
  for ( c = 0; c < len / 8; c++ )
    if ( ( vec1->bits[c] ^ vec2->bits[c] ) != 255 )
      return FALSE;
  for ( c = len - len % 8; c < len; c++ )
    if ( ( bitstring_get_bit ( vec1, c ) 
	   && bitstring_get_bit ( vec2, c ) )
	 || ( ! bitstring_get_bit ( vec1, c ) 
	      && ! bitstring_get_bit ( vec2, c ) ) )
      return FALSE;
  return TRUE;
}

unsigned int  bitstring_count_bits ( bitstring_pt  vec )
{
  static unsigned char	no_of_bits[] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2,
					 2, 3, 2, 3, 3, 4 };
  int		        i, 
                        count = 0, 
                        nr;

  for ( i = 0; i < vec->len/8; i++ )
    {
      nr = vec->bits[i];
      count += (no_of_bits[(nr & 0xf0) >> 4] + no_of_bits[(nr & 0x0f)]);
    }
  nr = vec->bits[vec->len / 8];
  if ( vec->len % 8 > 0 )
    {
      nr = vec->bits[vec->len / 8];
      nr &= ~( (128 >> (vec->len % 8 - 1)) - 1 );
      count += (no_of_bits[(nr & 0xf0) >> 4] + no_of_bits[(nr & 0x0f)]);
    }
  
  return count;
}

void  bitstring_cpy_all ( bitstring_pt  vec, 
			  bitstring_pt  vec2 )
{
  unsigned int  c, len_c;
  
  len_c = min_int ( vec->len / 8 + 1, vec2->len / 8 + 1 );
  for ( c = 0; c < len_c; c++ )
    vec->bits[c] = vec2->bits[c];
}

bitstring_pt  bitstring_create ( unsigned int  len )
{
  bitstring_pt  vec;

  vec = (bitstring_pt) mymalloc ( sizeof (bitstring_t) );
  vec->bits = (unsigned char *) mymalloc ( len / 8 + 1);
  vec->len  = len;
  return vec;
}

void  bitstring_free ( bitstring_pt  vec )
{
  free ( vec->bits );
  free ( vec );
}

void  bitstring_or_all ( bitstring_pt  vec, 
			 bitstring_pt  vec2 )
{
  unsigned int  c, len_c;
  
  len_c = min_int ( vec->len / 8 + 1, vec2->len / 8 + 1 );
  for ( c = 0; c < len_c; c++ )
    vec->bits[c] |= vec2->bits[c];
}

void  bitstring_print ( bitstring_pt  vec )
{
  int  i;
  
  for ( i = 0; i < vec->len; i++ )
    if ( bitstring_get_bit ( vec, i ) == 1 )
      printf ( "1" );
    else
      printf ( "0" );
  printf ( "\n" );
}

void  bitstring_print_part ( bitstring_pt  vec,
			     int           length )
{
  int  i;
  
  for ( i = 0; i < length; i++ )
    if ( bitstring_get_bit ( vec, i ) == 1 )
      printf ( "1" );
    else
      printf ( "0" );
  printf ( "\n" );
}

void  bitstring_set_all ( bitstring_pt  vec )
{
  unsigned int  c;
  
  for ( c = 0; c < vec->len / 8 + 1; c++ )
    vec->bits[c] = 255;
}

int  bitstring_set_bit ( bitstring_pt  vec, 
			 unsigned int  pos )
{
  if ( bitstring_get_bit ( vec, pos ) == 0 )
    {
      vec->bits[pos/8] += 1 << ( 7 - pos % 8 );
      return 1;
    }
  else
    return -1;
}

void  bitstring_toggle_all ( bitstring_pt  vec )
{
  unsigned int  c;

  for ( c = 0; c < vec->len / 8 + 1; c++ )
    vec->bits[c] ^= 255;
}

int  bitstring_toggle_bit( bitstring_pt vec, 
			   unsigned int pos )
{
  int bit;

  bit = bitstring_get_bit ( vec, pos );
  if ( bit == 0 )
    vec->bits[pos/8] += 1 << ( 7 - pos % 8 );
  else 
    if ( bit == 1 )
      vec->bits[pos/8] -= 1 << ( 7 - pos % 8 );
  return bit;
}

void  bitstring_xor_all ( bitstring_pt  vec, 
			  bitstring_pt  vec2 )
{
  unsigned int  c, len_c;

  len_c = min_int ( vec->len / 8 + 1, vec2->len / 8 + 1 );
  for ( c = 0; c < len_c; c++ )
    vec->bits[c] ^= vec2->bits[c];
}

