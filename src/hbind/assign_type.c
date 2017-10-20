#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "err_handle.h"

int  assign_type_and_orbit ( char  *str, int   *type, int   *orbit)
{
  char  *orbit_str;
  char  err_msg[FILENAME_MAX];

  /* check if there is a '.' that separates both parts */
  orbit_str = strchr ( str, '.' );
  if ( orbit_str == str || *str == '*' )
    /* no atom type specified */
    *type = ANY;
  else
    /* determine atom type */
    /* ADD MOL2 ATOM: Add the parsing code for the new atom based on
       the atom name. Also, add the vdW radius to inc/vdwrad.h (at 
       'ADD MOL2 ATOM'), the radius parsing code in read_mol2.c 
       (at 'ADD MOL2 ATOM'), and a type definition in inc/defs.h 
       (at 'ADD MOL2 ATOM') */
    {
      if ( strncmp ( str, "Ag", 2 ) == 0 ) *type = AG;
      else if ( strncmp ( str, "As", 2 ) == 0 ) *type = AS;
      else if ( strncmp ( str, "Al", 2 ) == 0 ) *type = AL;
      else if ( strncmp ( str, "Au", 2 ) == 0 ) *type = AU;
      else if ( strncmp ( str, "Ba", 2 ) == 0 ) *type = BA;
      else if ( strncmp ( str, "Be", 2 ) == 0 ) *type = BE;
      else if ( strncmp ( str, "Bi", 2 ) == 0 ) *type = BI;
      else if ( strncmp ( str, "Br", 2 ) == 0 ) *type = BR;
      else if ( *str == 'B' ) *type = B;
      else if ( strncmp ( str, "Ca", 2 ) == 0 ) *type = CA;
      else if ( strncmp ( str, "Cd", 2 ) == 0 ) *type = CD;
      else if ( strncmp ( str, "Ce", 2 ) == 0 ) *type = CE;
      else if ( strncmp ( str, "Co", 2 ) == 0 ) *type = CO;
      else if ( strncmp ( str, "Cl", 2 ) == 0 ) *type = CL;
      else if ( strncmp ( str, "Cs", 2 ) == 0 ) *type = CS;
      else if ( strncmp ( str, "Cu", 2 ) == 0 ) *type = CU;
      else if ( strncmp ( str, "Cr", 2 ) == 0 ) *type = CR;
      else if ( *str == 'C' ) *type = C;
      else if ( strncmp ( str, "Er", 2 ) == 0 ) *type = ER;
      else if ( strncmp ( str, "Fe", 2 ) == 0 ) *type = FE;
      else if ( *str == 'F' ) *type = F;
      else if ( strncmp ( str, "Ga", 2 ) == 0 ) *type = GA;
      else if ( strncmp ( str, "Ge", 2 ) == 0 ) *type = GE;
      else if ( strncmp ( str, "Hg", 2 ) == 0 ) *type = HG;
      else if ( *str == 'H' ) *type = H;
      else if ( strncmp ( str, "In", 2 ) == 0 ) *type = IN;
      else if ( *str == 'I' ) *type = I;
      else if ( *str == 'K' ) *type = K;
      else if ( strncmp ( str, "La", 2 ) == 0 ) *type = LA;
      else if ( strncmp ( str, "Li", 2 ) == 0 ) *type = LI;
      else if ( strncmp ( str, "Mg", 2 ) == 0 ) *type = MG;
      else if ( strncmp ( str, "Mn", 2 ) == 0 ) *type = MN;
      else if ( strncmp ( str, "Mo", 2 ) == 0 ) *type = MO;
      else if ( strncmp ( str, "Na", 2 ) == 0 ) *type = NA;
      else if ( strncmp ( str, "Ni", 2 ) == 0 ) *type = NI;
      else if ( *str == 'N' ) *type = N;
      else if ( strncmp ( str, "Os", 2 ) == 0 ) *type = OS;
      else if ( *str == 'O' ) *type = O;
      else if ( strncmp ( str, "Pd", 2 ) == 0 ) *type = PD;
      else if ( *str == 'P' ) *type = P;
      else if ( strncmp ( str, "Rb", 2 ) == 0 ) *type = RB;
      else if ( strncmp ( str, "Re", 2 ) == 0 ) *type = RE;
      else if ( strncmp ( str, "Ru", 2 ) == 0 ) *type = RU;
      else if ( strncmp ( str, "Rh", 2 ) == 0 ) *type = RH;
      else if ( strncmp ( str, "Sb", 2 ) == 0 ) *type = SB;
      else if ( strncmp ( str, "Se", 2 ) == 0 ) *type = SE;
      else if ( strncmp ( str, "Si", 2 ) == 0 ) *type = SI;
      else if ( strncmp ( str, "Sn", 2 ) == 0 ) *type = SN;
      else if ( strncmp ( str, "Sr", 2 ) == 0 ) *type = SR;
      else if ( *str == 'S' ) *type = S;
      else if ( strncmp ( str, "Te", 2 ) == 0 ) *type = TE;
      else if ( *str == 'T' ) *type = TL;
      else if ( *str == 'U' ) *type = U;
      else if ( *str == 'V' ) *type = V;
      else if ( *str == 'Y' ) *type = Y_;
      else if ( strncmp ( str, "Zn", 2 ) == 0 ) *type = ZN;
      else if ( strncmp ( str, "Zr", 2 ) == 0 ) *type = ZR;
      else {
	  if ( strncmp ( str, "Du", 2 ) != 0 )
	    /* don't panic if it is just a dummy */
	    {
	      *type = UNKNOWN;
	      fprintf ( stderr,
			"WARNING: atom definition: %s\n",
			str );
	      printf (	"WARNING: atom definition: %s\n",
			str );
	      sprintf ( err_msg, "WARNING: atom definition: %s\n",
			str);
	      err_print (err_msg);
	      err_warning2 ( "assign_type_and_orbit", "unknown atom type");
	    }
	  return FAILURE;
	}
    }
  /* no orbit type specified */
  if ( orbit_str == NULL ) *orbit = ANY;
  /* determine orbit type */
  else{
    orbit_str++;   /* skip '.' */
    /* this is for rules like "( . )", i.e. something has to be
       bound to this, but neither atom type nor orbit are specfied,
       and for rules like "N.", where any orbit is ok */
    if ( *orbit_str == ' ' || *orbit_str == '\0' ) *orbit = ANY;
      else if ( *orbit_str == '1' ) *orbit = SP1;
      else if ( *orbit_str == '2' ) *orbit = SP2;
      else if ( *orbit_str == '3' ) *orbit = SP3;
      else if ( *orbit_str == '4' ) *orbit = SP4;
      else if ( strncmp ( orbit_str, "ar", 2 ) == 0 ) *orbit = AR;
      else if ( strncmp ( orbit_str, "cat", 3 ) == 0 ) *orbit = CAT;
      else if ( strncmp ( orbit_str, "co2", 3 ) == 0 ) *orbit = CO2;
      else if ( strncmp ( orbit_str, "oh", 2 ) == 0 ) *orbit = OH;
      else if ( strncmp ( orbit_str, "o2", 2 ) == 0 ) *orbit = O2;
      else if ( strncmp ( orbit_str, "am", 2 ) == 0 ) *orbit = AM;
      else if ( strncmp ( orbit_str, "pl3", 3 ) == 0 ) *orbit = PL3;
      else if ( strncmp ( orbit_str, "th", 2 ) == 0 ) *orbit = TH;
      else if ( *orbit_str == 'o' ) *orbit = O;
      else {
        sprintf ( err_msg, "atom definition: %s\n", str );
        fprintf ( stderr, "atom definition: %s\n", str );
        printf (  "atom definition: %s\n", str );
        err_print (err_msg);
        err_warning2 ( "assign_type_and_orbit", "unknown orbit type");
        return FAILURE;
      }
    }
  return SUCCESS;
}

