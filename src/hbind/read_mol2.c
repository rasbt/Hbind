/*
 * $Source: /psa/share/repository/hbind/src/hbind/read_mol2.c,v $
 * $Revision: 1.11 $
 * $Author: toneroma $
 * $Date: 2009/06/29 18:57:23 $
 *
 * $Log: read_mol2.c,v $
 * Revision 1.11  2009/06/29 18:57:23  toneroma
 * added the ability to track molecule & conformer correctly so that the group_conformers will work properly (and so naming is consistent no matter the options set)
 *
 * Revision 1.10  2009/05/19 17:23:43  vanvoor4
 * Logic error -- forgot to set is_conf to FALSE when !conf -- Doh!
 *
 * Revision 1.9  2009/05/13 14:04:40  vanvoor4
 * Repaired bugs dealing with reduction in mirror variables in the global
 * structure.
 *
 * Revision 1.8  2009/03/19 14:52:55  vanvoor4
 * Added initial values to variables that were used later but not
 * always initialized
 *
 * Revision 1.7  2009/02/26 20:55:29  vanvoor4
 * Added dummy values to atom.residue, atom.residue_num, etc.
 * NO need to pass global to err_msg, etc
 *
 * Revision 1.6  2008/09/02 15:14:28  toneroma
 * changes to allow more atom types that start with the same letter (ie B* can be B, Be, Br, etc)
 * changes to remove unused variables
 * changes to initialize variables
 * group conformers changes
 * changes for restarting runs
 *
 * Revision 1.5  2007/10/12 16:49:42  toneroma
 * Added support for S.o type atoms (previously only accepted it if it was S.O)
 *
 * Revision 1.4  2007/10/09 21:32:20  toneroma
 * Now saving conf number so when group_conformers not selected each one will be output separately
 *
 * Revision 1.3  2007/09/28 18:33:48  toneroma
 * *** empty log message ***
 *
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"
#include "vdwrad.h"
#include <errno.h>

/* #define TRACE  */

void error_orbit(char *code)
{
  char err_msg[256];
  sprintf(err_msg, "WARNING: unknown orbit type; atom definition: %s\n", code);
  fprintf(stderr, "%s", err_msg);
  printf("%s", err_msg);
  err_print(err_msg);
  err_warning2 ( "assign_type_orbit_radius", "unknown atom.orbit type");
}

void error_mol2_section(char *expected, char *found, char *filename)
{
  char *help;
  char msg[256];
  help = strrchr(found, '\n');
  if(*help == '\n') *help = '\0';
  sprintf(msg, "Expected a(an) %s section in (%s)\nFound \"%s\" instead\n",
          expected, filename, found);
  err_warning2("read_mol2", msg);
}

void error_mol2_line(char *expected, char *found, char *filename)
{
  char *help;
  char msg[256];
  help = strrchr(found, '\n');
  if(*help == '\n') *help = '\0';
  sprintf(msg, "Expected a(an) %s record in (%s)\nFound \"%s\" instead\n",
          expected, filename, found);
  err_warning2("read_mol2", msg);
}

char *mol2_getline(char *s, int size, FILE *stream)
{
  char *rv;
  rv = fgets(s, size, stream);
  if(rv) return rv;

  if(feof(stream))
    err_warning2("mol2_getline", "Unexpected end of file");
  else err_warning2("mol2_getline", "Error reading a line from file");
  return NULL;
}


int assign_type_orbit_radius(char *str, int *type, int *orbit, float *radius)
{
  char  *orbit_str;
  char  err_msg[FILENAME_MAX];
  char tmp = 0;
  int i;

  *orbit = ANY;
  *radius = -1.0;

  /* check if there is a '.' that separates both parts */
  orbit_str = strchr ( str, '.' );
  if ( orbit_str == str || *str == '*' )
    /* no atom type specified */
    *type = ANY;
  else{
    if(orbit_str)
      for(i = orbit_str - str; i < strlen(str); ++i)
        orbit_str[i] = tolower(orbit_str[i]);

    /* determine atom type */
    /* ADD MOL2 ATOM: Add the parsing code for the new atom based on
       the atom name. Also, add the vdW radius to inc/vdwrad.h (at
       'ADD MOL2 ATOM'), the radius parsing code in read_mol2.c
       (at 'ADD MOL2 ATOM'), and a type definition in inc/defs.h
       (at 'ADD MOL2 ATOM') */
      if ( strncmp ( str, "Ag", 2 ) == 0 ){
        *type = AG;
        *radius = MRAD_AG;
      }else if ( strncmp ( str, "As", 2 ) == 0 ){
        *type = AS;
        *radius = MRAD_AS;
      }else if ( strncmp ( str, "Al", 2 ) == 0 ){
        *type = AL;
        *radius = MRAD_AL;
      }else if ( strncmp ( str, "Au", 2 ) == 0 ){
        *type = AU;
        *radius = MRAD_AU;
      }else if ( strncmp ( str, "Ba", 2 ) == 0 ){
        *type = BA;
        *radius = MRAD_BA;
      }else if ( strncmp ( str, "Be", 2 ) == 0 ){
        *type = BE;
        *radius = MRAD_BE;
      }else if ( strncmp ( str, "Bi", 2 ) == 0 ){
        *type = BI;
        *radius = MRAD_BI;
      }else if ( strncmp ( str, "Br", 2 ) == 0 ){
        *type = BR;
        *radius = MRAD_BR;
      }else if ( *str == 'B' ){
        *type = B;
        *radius = MRAD_B;
      }else if ( strncmp ( str, "Ca", 2 ) == 0 ){
        *type = CA;
        *radius = MRAD_CA;
      }else if ( strncmp ( str, "Cd", 2 ) == 0 ){
        *type = CD;
        *radius = MRAD_CD;
      }else if ( strncmp ( str, "Ce", 2 ) == 0 ){
        *type = CE;
        *radius = MRAD_CE;
      }else if ( strncmp ( str, "Co", 2 ) == 0 ){
        *type = CO;
        *radius = MRAD_CO;
      }else if ( strncmp ( str, "Cl", 2 ) == 0 ){
        *type = CL;
        *radius = MRAD_CL;
      }else if ( strncmp ( str, "Cs", 2 ) == 0 ){
        *type = CS;
        *radius = MRAD_CS;
      }else if ( strncmp ( str, "Cu", 2 ) == 0 ){
        *type = CU;
        *radius = MRAD_CU;
      }else if ( strncmp ( str, "Cr", 2 ) == 0 ){
        *type = CR;
        *radius = MRAD_CR;

        /* Assumes we are reading from a mol2 file */
        if(orbit_str){
          if(strncmp(orbit_str + 1, "th", 2) == 0) *orbit = TH;
          else if(strncmp(orbit_str + 1, "oh", 2) == 0) *orbit = OH;
          else{
            error_orbit(str);
            return FAILURE;
          }
        }
      }else if ( *str == 'C' ){
        *type = C;
        *radius = MRAD_CDEF;

        /* Assumes we are reading from a mol2 file */
        if(orbit_str){
          tmp = orbit_str[1];
          if(tmp == '1'){
            *orbit = SP1;
            *radius = MRAD_CSP1;
          }else if(tmp == '2'){
            *orbit = SP2;
            *radius = MRAD_CSP2;
          }else if(tmp == '3'){
            *orbit = SP3;
            *radius = MRAD_CSP3;
          }else if(tmp == '4'){
            *orbit = SP4;
          }else if(strncmp(orbit_str + 1, "ar", 2) == 0){
            *orbit = AR;
            *radius = MRAD_CAR;
          }else if(strncmp(orbit_str + 1, "cat", 3) == 0){
            *orbit = CAT;
            *radius = MRAD_CCAT;
          }else{
            error_orbit(str);
            return FAILURE;
          }
        }
      }else if ( strncmp ( str, "Er", 2 ) == 0 ){
        *type = ER;
        *radius = MRAD_ER;
      }else if ( strncmp ( str, "Fe", 2 ) == 0 ){
        *type = FE;
        *radius = MRAD_FE;
      }else if ( *str == 'F' ){
        *type = F;
        *radius = MRAD_F;
      }else if ( strncmp ( str, "Ga", 2 ) == 0 ){
        *type = GA;
        *radius = MRAD_GA;
      }else if ( strncmp ( str, "Ge", 2 ) == 0 ){
        *type = GE;
        *radius = MRAD_GE;
      }else if ( strncmp ( str, "Hg", 2 ) == 0 ){
        *type = HG;
        *radius = MRAD_HG;
      }else if ( *str == 'H' ){
        *type = H;
        *radius = MRAD_H;
      }else if ( strncmp ( str, "In", 2 ) == 0 ){
        *type = IN;
        *radius = MRAD_IN;
      }else if ( *str == 'I' ){
        *type = I;
        *radius = MRAD_I;
      }else if ( *str == 'K' ){
        *type = K;
        *radius = MRAD_K;
      }else if ( strncmp ( str, "La", 2 ) == 0 ){
        *type = LA;
        *radius = MRAD_LA;
      }else if ( strncmp ( str, "Li", 2 ) == 0 ){
        *type = LI;
        *radius = MRAD_LI;
      }else if ( strncmp ( str, "Mg", 2 ) == 0 ){
        *type = MG;
        *radius = MRAD_MG;
      }else if ( strncmp ( str, "Mn", 2 ) == 0 ){
        *type = MN;
        *radius = MRAD_MN;
      }else if ( strncmp ( str, "Mo", 2 ) == 0 ){
        *type = MO;
        *radius = MRAD_MO;
      }else if ( strncmp ( str, "Na", 2 ) == 0 ){
        *type = NA;
        *radius = MRAD_NA;
      }else if ( strncmp ( str, "Ni", 2 ) == 0 ){
        *type = NI;
        *radius = MRAD_NI;
      }else if ( *str == 'N' ){
        *type = N;
        *radius = MRAD_NDEF;

        /* Assumes we are reading from a mol2 file */
        if(orbit_str){
          tmp = orbit_str[1];
          if(tmp == '3'){
            *orbit = SP3;
            *radius = MRAD_NSP3;
          }else if(tmp == '4'){
            *orbit = SP4;
            *radius = MRAD_NSP4;
          }else if(tmp == '1' ) *orbit = SP1;
          else if(tmp == '2' ) *orbit = SP2;
          else if(!strncmp(orbit_str + 1, "ar", 2)) *orbit = AR;
          else if(!strncmp(orbit_str + 1, "am", 2)) *orbit = AM;
          else if(!strncmp(orbit_str + 1, "pl3", 3)) *orbit = PL3;
          else{
            error_orbit(str);
            return FAILURE;
          }
        }
      }else if ( strncmp ( str, "Os", 2 ) == 0 ){
        *type = OS;
        *radius = MRAD_OS;
      }else if ( *str == 'O' ){
        *type = O;
        *radius = MRAD_ODEF;

        /* Assumes we are reading from a mol2 file */
        if(orbit_str){
          tmp = orbit_str[1];
          if(tmp == '3'){
            *orbit = SP3;
            *radius = MRAD_OSP3;
          }else if(tmp == '2') *orbit = SP2;
          else if(!strncmp(orbit_str + 1, "co2", 2)) *orbit = CO2;
          else{
            error_orbit(str);
            return FAILURE;
          }
        }
      }else if ( strncmp ( str, "Pd", 2 ) == 0 ){
        *type = PD;
        *radius = MRAD_PD;
      /* Assume the the orbit for mol2 Phosphorus is always sp3 ("P.3") */
      }else if ( *str == 'P' ){
        *type = P;
        *radius = MRAD_P;
        *orbit = SP3;
      }else if ( strncmp ( str, "Rb", 2 ) == 0 ){
        *type = RB;
        *radius = MRAD_RB;
      }else if ( strncmp ( str, "Re", 2 ) == 0 ){
        *type = RE;
        *radius = MRAD_RE;
      }else if ( strncmp ( str, "Ru", 2 ) == 0 ){
        *type = RU;
        *radius = MRAD_RU;
      }else if ( strncmp ( str, "Rh", 2 ) == 0 ){
        *type = RH;
        *radius = MRAD_RH;
      }else if ( strncmp ( str, "Sb", 2 ) == 0 ){
        *type = SB;
        *radius = MRAD_SB;
      }else if ( strncmp ( str, "Se", 2 ) == 0 ){
        *type = SE;
        *radius = MRAD_SE;
      }else if ( strncmp ( str, "Si", 2 ) == 0 ){
        *type = SI;
        *radius = MRAD_SI;
      }else if ( strncmp ( str, "Sn", 2 ) == 0 ){
        *type = SN;
        *radius = MRAD_SN;
      }else if ( strncmp ( str, "Sr", 2 ) == 0 ){
        *type = SR;
        *radius = MRAD_SR;
      }else if ( *str == 'S' ){
        *type = S;
        *radius = MRAD_S;

        /* Assumes we are reading from a mol2 file */
        if(orbit_str){
          tmp = orbit_str[1];
          if(tmp == '3') *orbit = SP3;
          else if(tmp == '2') *orbit = SP2;
          else if(!strncmp(orbit_str + 1, "o2", 2)) *orbit = O2;
          else if((tmp == 'O')||(tmp == 'o')) *orbit = O;
          else{
            error_orbit(str);
            return FAILURE;
          }
        }
      }else if ( strncmp ( str, "Te", 2 ) == 0 ){
        *type = TE;
        *radius = MRAD_TE;
      }else if ( *str == 'T' ){
        *type = TL;
        *radius = MRAD_TL;
      }else if ( *str == 'U' ){
        *type = U;
        *radius = MRAD_U;
      }else if ( *str == 'V' ){
        *type = V;
        *radius = MRAD_V;
      }else if ( *str == 'Y' ){
        *type = Y_;
        *radius = MRAD_Y_;
      }else if ( strncmp ( str, "Zn", 2 ) == 0 ){
        *type = ZN;
        *radius = MRAD_ZN;
      }else if ( strncmp ( str, "Zr", 2 ) == 0 ){
        *type = ZR;
        *radius = MRAD_ZR;
      /* don't panic if it is just a dummy */
      }else if(strncmp(str, "Du", 2) == 0){
        *type = DUMMY;
        if(orbit_str && orbit_str[1] == 'C') *orbit = C;
      }else{
        fprintf(stderr, "WARNING: unknown atom type; atom definition: %s\n",
                str);
        printf("WARNING: unknown atom type; atom definition: %s\n",
               str);
        sprintf(err_msg, "WARNING: unknown atom type; atom definition: %s\n",
                str);
        err_print(err_msg);
        err_warning2("assign_type_orbit_radius", "unknown atom.orbit type");
        return FAILURE;
      }
    }

  return SUCCESS;
}

FILE* open_mol2(char *filename)
{
  FILE  *in;
  char  linebuffer[MAX_MOL2_LINELENGTH];

  if((in = fopen(filename, "r")) == NULL ){
    sprintf(linebuffer, "unable to open file %s; %s", filename,
            strerror(errno));
    err_warning2("open_mol2", linebuffer);
    return NULL;
  }
  return in;
}

int read_mol2(FILE *in, char *filename, global_data_pt global, char *needle)
{
  char  linebuffer[MAX_MOL2_LINELENGTH];
  char  mol_name[MAX_LEN_MOL2_COMPOUND_NAME];
  char  mol_name_noconf[MAX_LEN_MOL2_COMPOUND_NAME];
  char  *conf;
  char  msg[80];
  int   i, j;
  char *help = NULL;
  int counter = 0;
  int num_read = 0;
  int num_atoms = 0;
  int num_bonds = 0;
  int num_subst = 0;
  int length_molname = 0;
  int length_conf = 0;
  int is_conf;

  molecule_pt molecule = NULL;;
  linebuffer[0] = '\0';
  msg[0] = '\0';
  mol_name[0] = '\0';
  mol_name_noconf[0] = '\0';
  is_conf = TRUE;


  if(!in){
    sprintf(linebuffer, "The input file was not opened");
    err_warning2("read_mol2", linebuffer);
    return FATAL_FAILURE;
  }
  molecule = global->ligand;
  molecule->name[0] = '\0';
  molecule->number_of_atoms = 0;
  molecule->number_of_bonds = 0;
  molecule->conf_number[0] = '\0';
  molecule->number_of_substructures = 0;

  /* Scan from the current position in the file till we get a molecule section.
     If needle is null, take the first such section as the one we want.
     Otherwise do linear scan till molecule name matches needle */
  while(mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in)){
    if(linebuffer[0] == '@'
       && strncmp(linebuffer, "@<TRIPOS>MOLECULE", 17 ) == 0){
      if(!mol2_getline(mol_name, MAX_LEN_MOL2_COMPOUND_NAME, in))
        return FATAL_FAILURE;

      /* Change the newline char to NULL or if it doesn't exist, bail out */
      if((help = strrchr(mol_name, '\n'))) *help = '\0';
      else{
        sprintf(msg, "Molecule name is longer than MAX_MOL2_LINELENGTH (%s)",
                filename );
        err_warning2("read_mol2", msg);
        return FAILURE;
      }

      if(!needle || !strcmp(needle, mol_name)) break;
      mol_name[0] = '\0';
    }
  }

  if(mol_name == '\0') return FATAL_FAILURE;

  /* In the future this block of code should be moved out since it is not
   * directly associated with reading of mol2 files
   */

  /* begin conf */
#ifndef OUTPUT_ALL_MATCHES
  /* We will not be grouping conformers */
  if(!global->group_conformers){
    strcpy (mol_name_noconf, mol_name);
    strcpy (molecule->name_noconf, mol_name);
    global->total_num_molecules_noconf++;
    /*    printf ("NG S :: %s : %s : %s : %s\n", mol_name, mol_name_noconf, molecule->name_noconf, global->ligand_file_name);*/
    /* We have a molecule name starting with "singleton" AND we want to group
     * conformers
     */
    /* !!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!
     * Molecules that contain multiple conformers MUST be part of the same
     * Multi-mol2 file, otherwise they will be considered "singleton"
     * molecules, and will not be grouped together if group_conformers is
     * set to 'True'.  (not sure about this.  Checking.  toneroma 20090601)
     */
    /* toneroma 20090601 - temporarly handling singleton and multimol2 the same
      }else if(!strncmp(global->compound_name, "singleton", 9)){
    sprintf(mol_name, "%s", global->ligand_file_name );
    sprintf(molecule->name_noconf, "%s", global->ligand_file_name );
    global->total_num_molecules_noconf++;
    printf ("singNG S  :: %s : %s : %s : %s ::\n", mol_name, mol_name_noconf, molecule->name_noconf, global->ligand_file_name);

  /* We have a molecule name that does not start with "singleton" AND
   * is therefore a Multimol2 file AND we want to group conformers
   */
  }else{
    /* !!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!
     * Molecules that contain multiple conformers MUST be part of the same
     * Multi-mol2 file, otherwise they will be considered "singleton"
     * molecules, and will not be grouped together if group_conformers is
     * set to 'True'.
    */

    /* We check for a '_' in the molecule name -- if there is one it is a
     * conformer, if there is not a '_' in the molecule name it is not a
     * conformer.  The characters after the '_' should be either strictly
     *  numeric up to the ".mol2" extension.  If any alphabetic characters
     * occur between the last '_' and the .mol2 extension, this will be
     * considered the 'xtal' conformer.  All molecules with the same
     * characters prior to the last '_' will be considered conformers of
     * the same molecule, so the user must be careful to follow the
     * conformer formatting of:
     * [Molecule name]_[conformer name/number]
     *
     *  !!!!!!!!!!!CAUTION!!!!!!!!!!!
     * If for some reason group_conformers is set to true, and there are
     * multiple ligands in a row in the same multi-mol2 file that have the
     * same characters before the last '_', but they are not part of the
     * same molecule (eg. ligand_1abc, ligand_1a3b, ligand_1vr1), these
     * will erroneously be considered to be the same molecule, named
     * "ligand", and they will overwrite each other so only one molecule
     * is written out.  Therefore, it is important that if group_conformers
     * is set to "True" that every molecule in the database has a conformer
     * designation, even if it is only '0' or 'xtal' or 'molecule' or
     * anything else, as long as there is SOMETHING after the final '_'
     * that is NOT part of the molecule name/designation
     */
    conf = strrchr(mol_name, '_');
    if(!conf){
      /* Multimol2 file, group conformers on, no conformer '_' detected
       */
      is_conf = FALSE;
      strcpy (mol_name_noconf, mol_name);
      strcpy (molecule->name_noconf, mol_name);
      global->total_num_molecules_noconf++;
      /*      printf ("GY S nonconf  :: %s : %s : %s : %s ::\n", mol_name, mol_name_noconf, molecule->name_noconf, global->ligand_file_name);*/
    /* We have a conformer -- (found a '_' in the file name at least) */
    }else{
      length_molname = strlen(mol_name);
      length_conf = strlen(conf);
      for(j = 0; j < length_conf; j++){
        if(isalpha(conf[j])){
          is_conf = FALSE;
          strcpy(mol_name_noconf, mol_name);
          strcpy(molecule->name_noconf, mol_name);
          strcpy(molecule->conf_number, "_xtal");
          global->total_num_molecules_noconf++;
	  /*	  printf ("GY S yesnconf  :: %s : %s : %s : %s ::\n", mol_name, mol_name_noconf, molecule->name_noconf, global->ligand_file_name);*/
          break;
        }
      }
    }

    if(is_conf == TRUE){
      /* Multimol2 file, group conformers on, conformer '_' IS detected
       */
      strncpy(molecule->name_noconf, mol_name, length_molname - length_conf);
      molecule->name_noconf[length_molname - length_conf] = '\0';
      strncpy(mol_name_noconf, mol_name, length_molname - length_conf);
      mol_name_noconf[length_molname - length_conf] = '\0';
      strncpy(molecule->conf_number, conf, length_conf);
      molecule->conf_number[length_conf] = '\0';

      /* The current molecule is NOT a conformer of the previous molecule */
      /* This probably should be checked in match_triangles, rather than here */
      if(strcmp(molecule->name_noconf, global->old_ligand_name_noconf)){
        global->total_num_molecules_noconf++;
	printf ("New Molecule\n");
      }
      /*      printf ("MM2 GY S yesnconf  :: %s : %s : %s : %s ::\n", mol_name, mol_name_noconf, molecule->name_noconf, global->ligand_file_name);*/
    }
  }
#endif /* end ifndef OUTPUT_ALL_MATCHES */

#ifdef OUTPUT_ALL_MATCHES
  strcpy (mol_name_noconf, mol_name);
  strcpy (molecule->name_noconf, mol_name);
  global->total_num_molecules_noconf++;
#endif
  /* end conf */

  /* Check if we are restarting a run */
  if(global->restart_molecule_check == TRUE){
    /*    printf ("restart? %s %s\n", molecule->name_noconf, global->restart_molecule);*/
    if(!strcmp(molecule->name_noconf, global->restart_molecule))
      global->restart_molecule_check = FALSE;
    else return RESTART_SKIP_MOL;
  }

  strcpy(molecule->name, mol_name);
  num_read = fscanf(in, "%d %d %d %*d", &num_atoms, &num_bonds, &num_subst);
#ifdef TRACE
  fprintf(stdout, "Number of items read from molecule line is %d\n", num_read);
#endif

  if(num_atoms > MAX_NUMBER_OF_MOL2_ATOMS ){
    sprintf(msg, "number of atoms is greater than "
            "MAX_NUMBER_OF_MOL2_ATOMS (%s)", filename );
    err_warning2("read_mol2", msg);
    return FAILURE;
  }else if(num_bonds > MAX_NUMBER_OF_MOL2_BONDS ){
    sprintf(msg, "number of bonds is greater than "
            "MAX_NUMBER_OF_MOL2_BONDS (%s)", filename );
    err_warning2( "read_mol2", msg);
    return FAILURE;
  }else if(num_subst > MAX_SUBSTS){
    sprintf(msg, "number of substructures is greater than "
            "MAX_SUBSTS (%s)", filename );
    err_warning2("read_mol2", msg);
    return FAILURE;
  }

  /* Get the atom section; if another section is found, backup and return */
  while(mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in)
        && linebuffer[0] != '@');
  if(!strncmp(linebuffer, "@<TRIPOS>ATOM", 13 ) == 0){
    fseek(in, strlen(linebuffer), SEEK_CUR);
    error_mol2_section("ATOM", linebuffer, filename);
    return FAILURE;
  }
  /* read num_atoms lines of atoms defs */
  for(i = 0, counter = 0; i < num_atoms; i++){
    /* we cannot use fscanf() here, since we only want to grep the first six
       items out of a line that can have six or more entries */
    if(!mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in))
      return FATAL_FAILURE;
    else if(linebuffer[0] == '@'){
      fseek(in, strlen(linebuffer), SEEK_CUR);
      error_mol2_line("ATOM", linebuffer, filename);
      return FAILURE;
    }
    num_read = sscanf(linebuffer, "%d %s %f %f %f %s %d %s %f",
                      &molecule->atoms[counter].atom_number,
                      molecule->atoms[counter].name,
                      &molecule->atoms[counter].pos[X],
                      &molecule->atoms[counter].pos[Y],
                      &molecule->atoms[counter].pos[Z],
                      molecule->atoms[counter].type_str,
                      &molecule->atoms[counter].subst_id,
                      molecule->atoms[counter].subst_name,
                      &molecule->atoms[counter].charge);
    molecule->atoms[counter].chain_id = ' ';
    molecule->atoms[counter].residue_num[0] = '\0';
    molecule->atoms[counter].residue[0] = '\0';

    /*    printf("%g %g %g\n",
	   molecule->atoms[counter].pos[X],
	   molecule->atoms[counter].pos[Y],
	   molecule->atoms[counter].pos[Z]); */

#ifdef TRACE
    fprintf(stdout, "Number of items read from atom line %d is %d\n",
            i, num_read);
#endif

    /* this is sort of confusing: there are three different numbers
       assigned to each atom in the ligand:
       - 'atom_number' is the number of the atom in the original mol2-file,
         usually the atoms are numbered consecutively, starting with 1, and
         in the bond entries the connected atoms are identified via this
         variable
       - 'number' is here initialized with the atom_number, but it will later
         be overwritten when we renumber the atoms during the fragmentation
         of the ligand, there all bonds have a direction, which will be used
         to decide on which side of a bond a fragment is located, and this
         direction is based on the 'number' of each atom
       - 'atom_index' is the reference from 'atom_number' to the actual index
         in the 'atoms[]' array, we need this to figure out what atoms are
         connected by a bond
       generally, we should be fine with only 'number', but it is likely that
       we cannot assume that atoms are numbered consecutively, and with these
       three variables we can handle mol2-files with random numbering.
       Furthermore, many of the tools used in the lab puke on
       non-consecutively numbered atoms or bonds or ones that do not start
       counting from 1.  Such as InsightII, pymol, etc */
    molecule->atoms[counter].number = molecule->atoms[counter].atom_number;
    molecule->atom_index[molecule->atoms[counter].number] = counter;
    molecule->atoms[counter].act = NOTHING;
    if(assign_type_orbit_radius(molecule->atoms[counter].type_str,
                                &molecule->atoms[counter].type,
                                &molecule->atoms[counter].orbit,
                                &molecule->atoms[counter].rad)
      == FAILURE){
      /* there is either an unknown atom type in this molecule, or it
         contains a dummy atom, in either case we just skip this compound */
         return FAILURE;
    }
    switch ( molecule->atoms[counter].type ){
    case O:
      molecule->atoms[counter].hydro = 530;
    break;
    case N:
      molecule->atoms[counter].hydro = 350;
      break;
    case C:
    case S:
      molecule->atoms[counter].hydro = 80;
      break;
    default:
      molecule->atoms[counter].hydro = 317;
      break;
    }
    counter++;
  }
  /* reassign the number of atoms, since there might have been some dummy
     atoms */
  molecule->number_of_atoms = counter;

  /* Get the bond section; if another section is found, backup and return */
  while(mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in)
        && linebuffer[0] != '@');
  if(!strncmp(linebuffer, "@<TRIPOS>BOND", 13 ) == 0){
    fseek(in, strlen(linebuffer), SEEK_CUR);
    error_mol2_section("BOND", linebuffer, filename);
    return FAILURE;
  }
  /* read num_bonds lines of bonds defs */
  for(i = 0, counter = 0; i < num_bonds; i++){
    if(!mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in))
      return FATAL_FAILURE;
    else if(linebuffer[0] == '@'){
      fseek(in, strlen(linebuffer), SEEK_CUR);
      error_mol2_line("BOND", linebuffer, filename);
      return FAILURE;
    }
    num_read = sscanf(linebuffer, "%d %d %d %s",
                      &molecule->bonds[counter].number,
                      &molecule->bonds[counter].atom1,
                      &molecule->bonds[counter].atom2,
                      molecule->bonds[counter].type_str);
#ifdef TRACE
    fprintf(stdout, "Number of items read from bond line %d is %d\n",
            i, num_read);
#endif
    /* convert file numbering into internal numbering */
    molecule->bonds[counter].atom1 =
      molecule->atom_index[molecule->bonds[counter].atom1];
    molecule->bonds[counter].atom2 =
      molecule->atom_index[molecule->bonds[counter].atom2];
    /* make shure that this is not a bond to a dummy atom */
    if(molecule->bonds[counter].atom1 != UNKNOWN
       && molecule->bonds[counter].atom2 != UNKNOWN){
      if ( molecule->bonds[counter].type_str[0] == '1' )
        molecule->bonds[counter].type = SINGLE;
      else if ( molecule->bonds[counter].type_str[0] == '2' )
        molecule->bonds[counter].type = DOUBLE;
      else if ( molecule->bonds[counter].type_str[0] == '3' )
        molecule->bonds[counter].type = TRIPLE;
      else if ( molecule->bonds[counter].type_str[0] == '4' )
        molecule->bonds[counter].type = QUADRUPLE;
      else if ( molecule->bonds[counter].type_str[0] == '7' )
        molecule->bonds[counter].type = DELOCALIZED;
      else if ( molecule->bonds[counter].type_str[0] == 'a' ||
         molecule->bonds[counter].type_str[0] == 'A' ) {
        if ( molecule->bonds[counter].type_str[1] == 'm' ||
             molecule->bonds[counter].type_str[1] == 'M')
          molecule->bonds[counter].type = AMIDE;
        else
          molecule->bonds[counter].type = AROMATIC;
      }
    /* it is  an unknown bond type, set the bond type
       to DOUBLE, so that this bond will not be rotated */
    else if(molecule->bonds[counter].type_str[0] == 'u' ||
            molecule->bonds[counter].type_str[0] == 'U' )
        molecule->bonds[counter].type = DOUBLE;
      else{
        sprintf(msg, "unknown bond type: %s (%s)",
                molecule->bonds[counter].type_str, filename );
        err_print(linebuffer);
        err_warning2("read_mol2", msg);
        return FAILURE;
      }
      counter++;
    }
  }
  /* reassign the number of bonds, since there might have been some dummy
      atoms */
  molecule->number_of_bonds = counter;

  if(num_subst){
    while(mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in)
          && linebuffer[0] != '@');
    if(!strncmp(linebuffer, "@<TRIPOS>SUBSTRUCTURE", 21 ) == 0){
      fseek(in, strlen(linebuffer), SEEK_CUR);
      error_mol2_section("SUBSTRUCTURE", linebuffer, filename);
      return FAILURE;
    }

    for(i=0; i < num_subst; i++)
      {
	if(!mol2_getline(linebuffer, MAX_MOL2_LINELENGTH, in))
          return FATAL_FAILURE;
	else if(linebuffer[0] == '@'){
	  fseek(in, strlen(linebuffer), SEEK_CUR);
	  error_mol2_line("SUBSTRUCTURE", linebuffer, filename);
	  return FAILURE;
	}
	strcpy (molecule->substructure[i], linebuffer);
	molecule->number_of_substructures = i+1;
      }
  }
  return SUCCESS;
}
