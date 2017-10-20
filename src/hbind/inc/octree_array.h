#ifndef OCTREE_ARRAY_HEADER_FILE_INCLUDED
#define OCTREE_ARRAY_HEADER_FILE_INCLUDED

#include <string.h>
#include <atom.h>


typedef enum{
  NO_INTERSECTION = 0,       /*!< Cube and sphere do not intersect nor is the
			          cube fully contained in the sphere */
  CENTER_IN_CUBE,            /*!< Sphere center is in the cube */
  CENTROID_IN_SPHERE,        /*!< Cube centroid is in the sphere */
  CUBE_AND_SPHERE_INTERSECT,  /*!< Cube and sphere intersect */
  POSSIBLE_INTERSECTION
}cube_sphere_intersect_t;

struct octree_node_struct; /*!< Forward declaration to allow pointer to 
			        octree_node_struct in octree_data_t */

typedef struct{
  float *position;
  atom_pt atom;
  struct octree_node_struct* orig_node; /*!< Node based on original position */
  struct octree_node_struct* cnode;     /*!< Node based on current position */
}*octree_data_p, octree_data_t;

typedef struct octree_node_struct{
  struct octree_node_struct* parent;
  struct octree_node_struct* children[8];
  float centroid[3];
  float half_width;  
  octree_data_p *data;
  size_t num_el;     /*!< Number of pointers stored in the data pointer array */
  octree_data_p *added_data;
  size_t num_added_el; /*!< Number of pointers stored in added_data */
  size_t added_size;   /*!< Max # of pters that can be stored in added_data */
  int leaf_node;       /*!< True imples leaf node, false is not a leaf node */
}*octree_node_p, octree_node_t;

/*! Simple Octree for any datatype with an associated 3D position
 *  It is reasonably efficient when compared with brute force, but a few more
 * checks and the number of returned points, etc would be reduced further.
 */
typedef struct{
  octree_node_t root_node;  /*!< Root node of the tree */
  float min_half_width;
  size_t max_numel;
  octree_data_p data;
  octree_data_p data_end;
}*octree_p, octree_t;

void init_octree_node(octree_node_p node);

void build_octree(octree_p tree, atom_pt atoms, float *positions, 
                  size_t num_atoms, float min_half_width, size_t max_el_in_bin);

void subdivide(octree_p tree, octree_node_p cnode, int level);

void free_octree(octree_p tree);

void free_octree_node(octree_node_p cnode);
/*
void get_close_atoms(const octree_p octree, const float *pt, 
                     const float pt_tol, atom_pt **close_atoms, 
                     size_t *num_atoms);

void get_close_points(const octree_p octree, const float *pt, 
		      const float pt_tol, float **close_pts_p, size_t *npts);
		      */

void find_close_atoms(const octree_node_p cnode, const float *pt, 
                      const float pt_tol, atom_pt **close_atoms, 
                       size_t *num_atoms, int level);

void find_close_points(const octree_node_p cnode, const float *pt, 
                       const float pt_tol, float ***close_points_p, 
		       size_t *npts, int level);

int point_in_cube(const float *centroid, const float half_width, 
                  const float *pt);

int cube_and_sphere_intersect(const float *centroid, const float half_width,
                              const float *center, const float radius);

void print_bin_trace(const octree_p otree, const octree_node_p cnode, 
		     int level);

void check_position(octree_p tree, const atom_pt atom);

/*! This method is ok for side chain unbumping since we do not have much
 * movement back and forth.  If this method is to be used for other 
 * purposes that have multiple steps with sidechains moving in and out of
 * boxes it would bog things down since it does not remove the point from
 * boxes it was in previously (this is to keep things easier at this time
 * to not have to support removing points from their original bin
 */
void _update_bin(octree_node_p cnode, octree_data_p atom_data);

void reset_positions(octree_p tree);

void zero_out_moved_positions(octree_node_p cnode);

/* Assumes *close_atoms_array_ptr does not point to previously allocated 
 * memory 
 *
 * Returns a dynamically allocated array of pointers to nearby atoms
 */ 
static inline void 
get_close_atoms(const octree_p octree, const float *pt, const float pt_tol, 
                atom_pt **close_atoms_array_ptr, size_t *num_atoms)
{
  *close_atoms_array_ptr = 0;
  *num_atoms = 0;
  find_close_atoms(&octree->root_node, pt, pt_tol, close_atoms_array_ptr, 
                   num_atoms, 0);
}

/*
 * Returns a dynamically allocated array of pointers to nearby points 
 */
static inline void 
get_close_points(const octree_p octree, const float *pt, const float pt_tol,
                 float ***close_pts_p, size_t *npts)
{
  *close_pts_p = 0;
  *npts = 0;
  find_close_points(&octree->root_node, pt, pt_tol, close_pts_p, npts, 0);
}

static inline int 
pt_in_cube(const float *centroid, const float half_width, const float *center, 
           const float radius)
{
  float dist[3];
  size_t i;

  for(i = 0; i < 3; ++i){
    dist[i] = centroid[i] - center[i];
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }

  /* If center in box, return true */
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width)
    return CENTER_IN_CUBE;

  if(dist[0] <= half_width + radius && dist[1] <= half_width + radius && 
     dist[2] <= half_width + radius) return POSSIBLE_INTERSECTION;

  /* If center not in box + tol, return false */
  return NO_INTERSECTION;
}

static inline void
update_bin(octree_p tree, const atom_pt atom)
{
#if 0
  fprintf(stderr, "\n\ntime to update the bin\n");
  fprintf(stderr, "current atom pos: %f %f %f\n", atom->pos[0], 
          atom->pos[1], atom->pos[2]);
#endif
  octree_data_p data_ptr = tree->data + (atom - tree->data->atom);
  //fprintf(stderr, "passed the pointer voodoo line\n");
  _update_bin(data_ptr->cnode, data_ptr);
}

#if 0
void
init_octree_node(octree_node_p node)
{
  node->parent = 0;
  memset(node->children, 0, 8 * sizeof(*node->children));
  memset(node->centroid, 0, 3 * sizeof(*node->centroid));
  node->half_width = 0;
  node->data = 0;
  node->num_el = 0;
  node->added_data = 0;
  node->num_added_el = 0;
  node->added_size = 0;
}

void 
build_octree(octree_p tree, atom_pt atoms, float *positions, size_t num_atoms, 
             float min_half_width, size_t max_el_in_bin)
{
  float max_pt[3], min_pt[3];
  const float *p;
  size_t i;
  float longest_side;
  float tmp;
  float *pos = 0;
  octree_data_p data_ptr = 0;
  octree_data_p *root_data_ptr = 0;
  octree_node_p root_node = &tree->root_node;
  float *positions_end = &positions[3*num_atoms];
  atom_pt atom = 0;
  /*float cube_width;*/

  tree->min_half_width = min_half_width;
  tree->max_numel = max_el_in_bin;
  init_octree_node(root_node);
  

  /* Spin through the positions to find the min and max values */
  memcpy(max_pt, positions, 3*sizeof(*positions));
  memcpy(min_pt, positions, 3*sizeof(*positions));
  for(p = positions; p < positions_end; p += 3)
    for(i = 0; i < 3; ++i){
      if(max_pt[i] < p[i]) max_pt[i] = p[i];
      if(min_pt[i] > p[i]) min_pt[i] = p[i];
    }

  /*
   * Use a square at the cost of more storage (larger tree) and possibly
   * a longer search time.  We do want a reasonable tolerance on the added
   * to each side so that points do not easily move outside
   */

  /* Compute the centroid and halfwidths */
  longest_side = 0;
  for(i = 0; i < 3; ++i){
    root_node->centroid[i] = 0.5 * (max_pt[i] + min_pt[i]);
    tmp = max_pt[i] - min_pt[i];
    longest_side = (longest_side < tmp ? tmp : longest_side);
  }

  /*
   * Should be an input variable at some point in time -- a Lys or Arg
   * sidechain can easily move more than 10.0 (A)
   */

  longest_side += 40.0;
  root_node->half_width = 0.5 * longest_side;
#if 0
  /* Tune size to be divisible by min half width otherwise we get 
     a less tight set of boxes -- need to test on a number of different
     point sets in order to get an idea of its impact*/
  cube_width = min_half_width;
  for( ; longest_side > min_half_width; cube_width *= 2.0, longest_side /= 2.0);
  root_node->half_width = 0.5 * cube_width;
#endif

  /* Fill the data pointer array */
  tree->data = (octree_data_p) mymalloc(num_atoms * sizeof(octree_data_t));
  tree->data_end = &tree->data[num_atoms];
  atom = atoms;
  pos = positions;
  for(data_ptr = tree->data; data_ptr < tree->data_end; ++data_ptr, ++atom){
    data_ptr->position = pos;
    data_ptr->atom = atom;
    data_ptr->orig_node = 0;
    data_ptr->cnode = 0;
    pos += 3;
  }

  /* Fill the array for the head node */
  root_node->data = (octree_data_p*)mymalloc(num_atoms * sizeof(octree_data_p));
  root_node->num_el = num_atoms; 
  root_data_ptr = root_node->data;
  for(data_ptr = tree->data; data_ptr < tree->data_end; ++data_ptr){
    *root_data_ptr = data_ptr;
    ++root_data_ptr;
  }

  subdivide(tree, root_node, 0);
}

void
subdivide(octree_p tree, octree_node_p cnode, int level)
{
  size_t idx;
  size_t j;
  size_t dim_mask;
  int i;
  const float *centroid = cnode->centroid;
  octree_data_p *data;
  octree_data_p *data_end = &(cnode->data[cnode->num_el]);
  const float child_h_width = 0.5 * cnode->half_width;
  octree_node_p child = 0;
  int octant_count[8];
  int *partition;
  int *octant;

/*
  printf("subdivision level %d\n", level);
  printf("half width %f\n", cnode->half_width);
  printf("centroid: %f %f %f\n", cnode->centroid[0], cnode->centroid[1], 
         cnode->centroid[2]);
  printf("number of positions: %d\n", cnode->num_el);
*/

  if(cnode->num_el <= tree->max_numel || child_h_width < tree->min_half_width){
    for(data = cnode->data; data < data_end; ++data) (*data)->orig_node = cnode;
    return;
  }

  memset(octant_count, 0, 8 * sizeof(*octant_count));
  partition = (int*) mymalloc(cnode->num_el * sizeof(int));
  octant = partition;

  /* partition based on the position of each point */
  for(data = cnode->data; data < data_end; ++data, ++octant){
    idx = 0;
    for(i = 2; i > -1; --i){
      idx <<= 1;
      if((*data)->position[i] >= centroid[i]) ++idx;
    }
    *octant = idx;
    ++octant_count[idx];
  }

  /* Allocate the children for the octants with data and the childrens' 
   * data arrays */
  for(idx = 0; idx < 8; ++idx){
    if(octant_count[idx] == 0) continue;

    child = (octree_node_p) mymalloc(sizeof(octree_node_t));
    init_octree_node(child);
    child->parent = cnode;
    child->half_width = child_h_width;
    dim_mask = 1; 
    for(j = 0; j < 3; ++j, dim_mask <<= 1){
      if(idx & dim_mask) child->centroid[j] = centroid[j] + child_h_width;
      else child->centroid[j] = centroid[j] - child_h_width;
    }
    child->data = 
      (octree_data_p*) mymalloc(octant_count[idx] * sizeof(octree_data_p));
    cnode->children[idx] = child;
  }

  /* Partition the data based on the partion array */
  octant = partition;
  for(data = cnode->data; data < data_end; ++data, ++octant){
    child = *(cnode->children + *octant);
    child->data[child->num_el] = *data;
    ++child->num_el;
  }
  free(partition);
  partition = 0;
  octant = 0;
  free(cnode->data);
  cnode->data = 0;
  cnode->num_el = 0;

  /* Recursively subdivide children */
  for(j = 0; j < 8; ++j)
    if(octant_count[j]) subdivide(tree, cnode->children[j], level + 1);
}

void free_octree(octree_p tree)
{
  free_octree_node(&tree->root_node);
}

void free_octree_node(octree_node_p cnode)
{
  int leaf_node = 0;
  int i;

  if(cnode->num_el){
    if(cnode->data) free(cnode->data);
    cnode->data = 0;
    cnode->num_el = 0;
    leaf_node = 1;
  }
  if(cnode->num_added_el){
    if(cnode->added_data) free(cnode->added_data);
    cnode->added_data = 0;
    cnode->num_added_el = 0;
    cnode->added_size = 0;
    leaf_node = 1;
  }
  if(leaf_node) return;

  for(i = 0; i < 8; ++i)
    if(cnode->children[i]){
      free_octree_node(cnode->children[i]);
      free(cnode->children[i]);
      cnode->children[i] = 0;
    }
}

/* Assumes *close_atoms_array_ptr is NULL and *num_atoms is 0 */
void 
get_close_points(const octree_p octree, const float *pt, const float pt_tol, 
                 atom_pt **close_atoms_array_ptr, size_t *num_atoms)
{
  find_close_points(&octree->root_node, pt, pt_tol, close_atoms_array_ptr, 
                    num_atoms);
}

void 
find_close_points(const octree_node_p cnode, const float *pt, 
                  const float pt_tol, atom_pt **close_atoms_array_ptr, 
		  size_t *num_atoms)
{
  int i;
  size_t new_num_atoms;
  octree_data_p *data_p;
  atom_pt *atom_p;
/*
    std::cout << "\nposition to find: " << pt[0] << " " << pt[1] << " " << pt[2] 
              << std::endl;
    std::cout << "level: " << level + 1 << "\n"
              << "size: " << cnode.positions.size() << "\n"
              << "h_width: " << cnode.half_width << " \n"
              << "cendroid: " << cnode.centroid[0] << " "
              << cnode.centroid[1] << " " << cnode.centroid[2] << "\n";
*/
  /* 
   * We could check distance and only push back those meeting the 
   * distance criterion here but that might cause the distance 
   * calculations for the points close to pt to be computed at least twice
   */
  if(cnode->num_el + cnode->num_added_el){
    new_num_atoms = *num_atoms + cnode->num_el + cnode->num_added_el;
    *close_atoms_array_ptr = 
      (atom_pt*) myrealloc(*close_atoms_array_ptr, 
                           new_num_atoms * sizeof(atom_pt));

    if(cnode->num_el){
      atom_p = &((*close_atoms_array_ptr)[*num_atoms]);
      data_p = cnode->data;
      for(; data_p < &(cnode->data[cnode->num_el]); ++data_p, ++atom_p)
        *atom_p = (*data_p)->atom; 
      *num_atoms += cnode->num_el;
    }
    if(cnode->num_added_el){
      atom_p = &((*close_atoms_array_ptr)[*num_atoms]);
      data_p = cnode->added_data;
      for(;data_p < &(cnode->added_data[cnode->num_added_el]); ++data_p){
        *atom_p = (*data_p)->atom; 
	++atom_p;
      }
      *num_atoms += cnode->num_added_el;
    }
    return;
  }

  for(i = 0; i < 8; ++i)
    if(cnode->children[i] && 
       cube_and_sphere_intersect(cnode->children[i]->centroid, 
                                 cnode->children[i]->half_width, pt, pt_tol))
      find_close_points(cnode->children[i], pt, pt_tol, close_atoms_array_ptr, 
                        num_atoms);
}

int
point_in_cube(const float *centroid, const float half_width, const float *pt)
{
  float dist[3];
  size_t i;

  for(i = 0; i < 3; ++i){
    dist[i] = centroid[i] - pt[i];     
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }

  /* If pt in box, return true */
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width) 
    return 1;

  return 0;
}

int 
cube_and_sphere_intersect(const float *centroid, const float half_width,
                          const float *center, const float radius)
{
  float dist[3];
  float pt_on_sphere[3];
  float pt2centroid[3];
  float mag_pt2cent;
  size_t i;

  for(i = 0; i < 3; ++i){
    dist[i] = centroid[i] - center[i];     
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }

  /* If center in box, return true */
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width) 
    return CENTER_IN_CUBE;

  /* If center not in box + tol, return false */
  if(dist[0] > half_width + radius && dist[1] > half_width + radius &&
     dist[2] > half_width + radius) 
    return NO_INTERSECTION;

  /* Project cube centroid to the sphere */
  mag_pt2cent = 0;
  for(i = 0; i < 3; ++i){
    pt2centroid[i] = centroid[i] - center[i];
    mag_pt2cent += pt2centroid[i]*pt2centroid[i];
  }
  mag_pt2cent = sqrt(mag_pt2cent);
  for(i = 0; i < 3; ++i){
    pt2centroid[i] /= mag_pt2cent;
    pt_on_sphere[i] = center[i] + radius * pt2centroid[i];
  }

  /* If projected centroid is in the box, return true, else return false  */
  for(i = 0; i < 3; ++i){  
    dist[i] = centroid[i] - pt_on_sphere[i];
    dist[i] = (dist[i] < 0 ? -1.0 * dist[i] : dist[i]);
  }
  if(dist[0] <= half_width && dist[1] <= half_width && dist[2] <= half_width) 
    return CUBE_AND_SPHERE_INTERSECT;
  return NO_INTERSECTION;
}

void
print_bin_trace(const octree_node_p cnode, int level)
{
  size_t i;
  if(cnode->num_el){
    printf("\nBin at level %d", level);
    printf("\nCentroid: %f %f %f", cnode->centroid[0], cnode->centroid[1], 
           cnode->centroid[2]);
    printf("\nHalf Width: %f", cnode->half_width);
    printf("\nPositions:\n");
    
    for(i = 0; i < cnode->num_el; ++i){
      printf("\t%f %f %f\n", cnode->data[i]->position[0],
             cnode->data[i]->position[1], cnode->data[i]->position[2]);
    }
    return;
  }

  for(i = 0; i < 8; ++i)
    if(cnode->children[i]) print_bin_trace(cnode->children[i], level + 1);
}
#endif

#endif
