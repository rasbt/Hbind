#ifndef  _BITSTRING_H
#define  _BITSTRING_H
 
extern int  min_int ( int  a, 
		      int  b );

extern int  bitstring_get_bit ( bitstring_pt  vec, 
				unsigned int  pos );
 
extern int  bitstring_all_clear ( bitstring_pt  vec );

extern int  bitstring_all_set ( bitstring_pt vec );

extern void  bitstring_and_all ( bitstring_pt  vec, 
				 bitstring_pt  vec2 );

extern void  bitstring_clear_all ( bitstring_pt  vec );

extern int  bitstring_clear_bit ( bitstring_pt  vec, 
				  unsigned int  pos );

extern int  bitstring_cmp_all ( bitstring_pt  vec, 
				bitstring_pt  vec2 );

extern int  bitstring_cmp_and_0 ( bitstring_pt  vec1, 
				  bitstring_pt  vec2 );

extern int  bitstring_cmp_or_1 ( bitstring_pt  vec1, 
				 bitstring_pt  vec2 );

extern int bitstring_cmp_xor_1 ( bitstring_pt  vec1, 
				 bitstring_pt  vec2 );

extern unsigned int  bitstring_count_bits ( bitstring_pt  vec );

extern void  bitstring_cpy_all ( bitstring_pt  vec, 
				 bitstring_pt  vec2 );

extern bitstring_pt  bitstring_create ( unsigned int  len );

extern void  bitstring_free ( bitstring_pt  vec );

extern void  bitstring_or_all ( bitstring_pt  vec, 
				bitstring_pt  vec2 );

extern void  bitstring_print ( bitstring_pt  vec );

extern void  bitstring_print_part ( bitstring_pt  vec,
				    int           number );

extern void  bitstring_set_all ( bitstring_pt  vec );

extern int  bitstring_set_bit ( bitstring_pt  vec, 
				unsigned int  pos );

extern void  bitstring_ptoggle_all ( bitstring_pt  vec );

extern int  bitstring_ptoggle_bit( bitstring_pt  vec, 
				   unsigned int  pos );

extern void  bitstring_xor_all ( bitstring_pt  vec, 
				 bitstring_pt  vec2 );

#endif
