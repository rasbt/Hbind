#include <stdlib.h>
#include <stdio.h>
#include "defs.h"
#include "types.h"
#include "basics.h"

void change_sort_array ( sort_array_pt item1, 
			 sort_array_pt item2)
{
    item1->index1 = item2->index1;
    item1->index2 = item2->index2;
    item1->value = item2->value;
}
 
void  insertion_sort ( sort_array_pt array,
		       int           length )
{
  sort_array_t  help_item;
  int           i = 0, j = 0;

  if ( length <= 1)
    return;

#ifdef TRACE
  printf("   insertion_sort: before sorting \n");
  for ( i = 0; i < length; i++ )
	printf("insertion_sort: i=%d, value=%.16f, index1=%d, index2=%d\n", i, array[i].value, array[i].index1, array[i].index2);
  printf("\n");
#endif

    for (i = 1; i < length; i++)
	{
	    change_sort_array( &help_item, &array[i] );
	    j = i;
	    while ( j > 0 && compare_double( array[j-1].value, help_item.value ) == -1 )
		{
		    change_sort_array( &array[j], &array[j-1] );
		    j--;
		}
	    change_sort_array( &array[j], &help_item );
	}

#ifdef TRACE
  for ( i = 0; i < length; i++ )
	printf("insertion_sort: i=%d, value=%.16f, index1=%d, index2=%d\n", i, array[i].value, array[i].index1, array[i].index2);
#endif

}

