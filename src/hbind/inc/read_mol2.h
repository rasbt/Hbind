#ifndef _READ_MOL2_H
#define _READ_MOL2_H

/*! /brief Assigns the type, orbit and vdw radius for given atom string
 * 
 * Assign the atom type, orbit (hybridization) and vdw radius for the given
 * SYBIL atom string.  
 *
 * @param str Input: string holding the atom name (PDB or mol2)
 * @param type Output: pointer to type
 * @param orbit Output: pointer to orbit
 * @param radius Output: pointer to radius
 * @param global_data_t Input: pointer to global stew
 * @return SUCCESS or FAILURE
 */
int assign_type_orbit_radius(char *str, int *type, int *orbit, float *radius);

FILE* open_mol2(char *filename);

/*! /brief Given a file pointer, if needle, get the molecule that matches
 *         needle; otherwise get the first molecule.
 *
 *  Start at the current position of the mol2 file in.  If needle is not null,
 *  do a linear scan (from current position) until we find a molecule whose
 *  name matches needle.  Otherwise, take the first molecule (from current
 *  position).
 *
 * @param in Input: pointer to file to read
 * @param filename Input: string holding the name of the open file (for error
 *        messages
 * @param global Input/Output: pointer to the global stew
 * @param needle Input: Null or or name of molecule to find
 * @return SUCCESS, FAILURE or FATAL_FAILURE
 */
int read_mol2(FILE* in, char *filename, global_data_pt global, char *needle);

#endif
