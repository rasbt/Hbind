#include <string.h>
#include <math.h>
#include <types.h>
#include <mymalloc.h>
#include <err_handle.h>
#include <dist_fun.h>


void
distance_array(dist_array_pt darray, const atom_pt atoms, int num_atoms, 
               const float box_width, const float max_dist)
{
  int i, j, k;
  float min_pt[3], max_pt[3];
  atom_pt a;
  const atom_pt atoms_end = atoms + num_atoms;
  dist_bin_pt bin_p;

  darray->box_width = box_width;

  /* Spin through the positions to find the min and max values */
  memcpy(min_pt, atoms->pos, 3*sizeof(min_pt[0]));
  memcpy(max_pt, atoms->pos, 3*sizeof(min_pt[0]));
  for(a = atoms; a < atoms_end; ++a){
    for(i = 0; i < 3; ++i){
      if(max_pt[i] < a->pos[i]) max_pt[i] = a->pos[i];
      if(a->pos[i] < min_pt[i]) min_pt[i] = a->pos[i];
    }
  }

  /* Make a box that has a buffer as well as side lengths divisible by 
   * box_width */
  for(i = 0; i < 3; ++i){
    max_pt[i] = box_width * ceil(max_pt[i] / box_width) + 2*box_width;
    min_pt[i] = box_width * floor(min_pt[i] / box_width) - 2*box_width;
    darray->num_bins[i] = (int) ((max_pt[i] - min_pt[i]) / box_width);
  }
  memcpy(darray->min_corner, min_pt, 3 * sizeof(min_pt[0]));
  memcpy(darray->max_corner, max_pt, 3 * sizeof(max_pt[0]));
 
  /* Populate the arrays */ 
  darray->bins = 
    (dist_bin_pt) mymalloc(darray->num_bins[0] * darray->num_bins[1] * 
                           darray->num_bins[2] * sizeof(dist_bin_t));
  bin_p = darray->bins;
  min_pt[0] = darray->min_corner[0];
  max_pt[0] = darray->min_corner[0] + box_width;
  for(i = 0; i < darray->num_bins[0]; ++i){
    min_pt[1] = darray->min_corner[1];
    max_pt[1] = darray->min_corner[1] + box_width;
    for(j = 0; j < darray->num_bins[1]; ++j){
      min_pt[2] = darray->min_corner[2];
      max_pt[2] = darray->min_corner[2] + box_width;
      for(k = 0; k < darray->num_bins[2]; ++k, ++bin_p){
        /* Intialize the bin variables */
        bin_p->atoms = 0;
        bin_p->num_atoms = 0;
        bin_p->size = 0;
        /* Populate the bin */
        for(a = atoms; a < atoms_end; ++a){
          if(min_pt[0] - max_dist <= a->pos[0] &&
             a->pos[0] <= max_pt[0] + max_dist &&
             min_pt[1] - max_dist <= a->pos[1] &&
             a->pos[1] <= max_pt[1] + max_dist && 
             min_pt[2] - max_dist <= a->pos[2] &&
             a->pos[2] <= max_pt[2] + max_dist){
            if(bin_p->size <= bin_p->num_atoms){
              bin_p->size = (bin_p->size ? 2*bin_p->size : 10);
              bin_p->atoms = (atom_pt *) myrealloc(bin_p->atoms, bin_p->size * 
                                                   sizeof(atom_pt));
            }
            bin_p->atoms[bin_p->num_atoms] = a;
            ++(bin_p->num_atoms);
          }
        }
        min_pt[2] += box_width;
        max_pt[2] += box_width;
      }
      min_pt[1] += box_width; 
      max_pt[1] += box_width; 
    }
    min_pt[0] += box_width;
    max_pt[0] += box_width;
  }
}

/*
 Need to fix this as it causes a segfault -- a bit confusing since my_free
 doesn't attempt to free null pointers.
*/
void
free_distance_array(dist_array_pt d)
{
  int i;
  int num_bins = d->num_bins[0] * d->num_bins[1] * d->num_bins[2];
  for(i = 0; i < num_bins; ++i)
    my_free(d->bins[i].atoms);
  my_free(d->bins);
}

dist_bin_pt
get_bin(dist_array_pt darray, float *pt)
{
  int i;
  int idx[3];

  /* Check if point falls outside the range of the bins */
  for(i = 0; i < 3; ++i)
    if(pt[i] < darray->min_corner[i] || darray->max_corner[i] < pt[i]) return 0;

  /* Compute the index */
  for(i = 0; i < 3; ++i)
    idx[i] = (int) floor((pt[i] - darray->min_corner[i]) / darray->box_width);
 
  return darray->bins + (idx[0] * darray->num_bins[1]*darray->num_bins[2] +
                         idx[1] * darray->num_bins[2] + idx[2]); 
}
