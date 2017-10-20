/* NOTE: This is the file that contains user adjustable parameters and definitions.  Please refer to the User Guide for recommended settings */


#undef  OUTPUT_ALL_MATCHES /* WARNING: turning on this option will result in MANY extra dockings, and may be unnecessary unless attempting to rescore dockings with an alternate scoring method.  When #undef it will output only the best scoring orientation for overlap free dockings and for dockings with tolerated overlap.  When #define it will output all dockings for overlap free dockings and for dockings with tolerated overlap */

#undef VERBOSE_OUTPUT /* When OUTPUT_ALL_MATCHES is undefined it is possible to still see the score output of the potential dockings in the stdout (or redirected file), if necessary. To turn on this option, simply set to #define.  This may, however, create a very large output file */

#define  FILTER_BURIED_CARBONS /* When #define it will output and score only dockings that have at least 50% of the ligand's carbons buried in the protein-ligand interface */

/* Define the ranges of edge lengths for the triangles (units are Angstroms) */
#define TRIANGLE_MIN_PERIMETER             7.0
#define TRIANGLE_MAX_PERIMETER            30.0
#define TRIANGLE_MIN_SHORTEST_SIDE         2.0
#define TRIANGLE_MAX_SHORTEST_SIDE        10.0
#define TRIANGLE_MAX_LONGEST_SIDE         12.0
/* Support smaller triangle sides for smaller ligands */
#define SMALL_TRIANGLE_MIN_LONGEST_SIDE    3.0
#define LARGE_TRIANGLE_MIN_LONGEST_SIDE    5.0

/* Define in Angstroms the diameter of the sphere in which a small ligand 
 * must fit
 */
#define SMALL_LIGAND_DIAMETER  10.0
