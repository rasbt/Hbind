#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>

#include "defs.h"
#include "types.h"
#include "basics.h"
#include "quicksort.h"

void swap( int in_array[], int i, int j);

void quicksort_int( int in_array[], int left, int right )
{
	int current=0, last=0;

	if (left >= right)
        	return;
	swap(in_array, left, (left+right)/2);
	last = left;
	for (current = left + 1; current <= right; current ++)
	    if ( in_array[current] < in_array[left] )
		swap(in_array, ++last, current);
	swap(in_array, left, last);
	quicksort_int(in_array, left, last-1);
	quicksort_int(in_array, last+1, right);
}

/** The function swap interchanges the values in two positions of an array. **/ 

void swap( int in_array[], int i, int j )
{
	int temp;

	temp = in_array[i];
	in_array[i] = in_array[j];
	in_array[j] = temp;
}


