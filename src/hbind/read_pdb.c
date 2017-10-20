#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "defs.h"
#include "types.h" 
#include "hydro.h"
#include "err_handle.h"
#include "read_mol2.h"
#include "vdwrad.h"

#define column(n)       (linebuf+(n)-1)
#define strkeq(s1,s2,k) ( k<=0 ? 1 : s1[0]==s2[0] ? !strncmp(s1,s2,(k)) : 0)
#define TRAC

float radius (int  type, int  orbit)
{
  float rad = -1.0;

  switch ( type ) {
    case AL:
      rad = MRAD_AL;
      break;
    case AG:
      rad = MRAD_AG;
      break;
    case AS:
      rad = MRAD_AS;
      break;
    case BR:
      rad = MRAD_BR;
      break;
    case C:
     rad = MRAD_CDEF;
     break;
    case CA:
      rad = MRAD_CA;
      break;
    case CD:
      rad = MRAD_CD;
      break;
    case CE:
      rad = MRAD_CE;
      break;
    case CL:
      rad = MRAD_CL;
      break;
    case CO:
      rad = MRAD_CO;
      break;
    case CR:
      rad = MRAD_CR;
      break;
    case CS:
      rad = MRAD_CS;
      break;
    case CU:
      rad = MRAD_CU;
      break;
    case ER:
      rad = MRAD_ER;
      break;
    case F: 
      rad = MRAD_F;
      break;
    case FE:
      rad = MRAD_FE;
      break;
    case GA:
      rad = MRAD_GA;
      break;
    case GE:
      rad = MRAD_GE;
      break;
    case H: 
      rad = MRAD_H;
      break;
    case I: 
      rad = MRAD_I;
      break;
    case K: 
      rad = MRAD_K;
      break;
    case LA:
      rad = MRAD_LA;
      break;
    case LI:
      rad = MRAD_LI;
      break;
    case MG:
      rad = MRAD_MG;
      break;
    case MN:
      rad = MRAD_MN;
      break;
    case MO:
      rad = MRAD_MO;
      break;
    case N:
      rad = MRAD_NDEF;
      break;
    case NA:
      rad = MRAD_NA;
      break;
    case NI:
      rad = MRAD_NI;
      break;
    case O:
      rad = MRAD_ODEF;
      break;
    case OS:
      rad = MRAD_OS;
      break;
    case P:
      rad = MRAD_P;
      break;
    case PD:
      rad = MRAD_PD;
      break;
    case RB:
      rad = MRAD_RB;
      break;
    case RE:
      rad = MRAD_RE;
      break;
    case RH:
      rad = MRAD_RH;
      break;
    case RU:
      rad = MRAD_RU; 
      break;
    case S:
      rad = MRAD_S;
      break; 
    case SI:
      rad = MRAD_SI;
      break;
    case TI:
      rad = MRAD_TI;
      break;
    case TL:
      rad = MRAD_TL;;
      break;
    case U:
      rad = MRAD_U;
      break;
    case V:
      rad = MRAD_V;
      break; 
    case Y_:
      rad = MRAD_Y_;
      break;
    case ZN:
      rad = MRAD_ZN;
      break;
    case ZR:
      rad = MRAD_ZR;
      break;
    case UNKNOWN:
      rad = MRAD_UNKNOWN;
      break;
    default:
      break;
    }
  return rad;
}

void  assign_residue_type ( atom_pt atm)
{
  switch ( atm->residue[0] )
    {
    case 'A':
      switch ( atm->residue[1] )
	{
	case 'C':
	  if ( strcmp ( atm->residue, "ACE" ) == 0 )
	    atm->res = ACE; 
	  break;
	case 'L':
	  atm->res = ALA;
	  break;
	case 'R':
	  atm->res = ARG;
	  break;
	case 'S':
	  switch ( atm->residue[2] )
	    {
	    case 'N':
	      atm->res = ASN;
	      break;
	    case 'P':
	      atm->res = ASP;
	      break;
	    default:
	      err_unknown_residue2 ( atm->residue, atm->residue_num);
	      break;
	    }
	  break;
	default:
	  err_unknown_residue2 ( atm->residue, atm->residue_num);
	  break;
	}
      break;
    case 'C':
      if ( strncmp ( atm->residue, "CYS", 3 ) == 0 )
	atm->res = CYS;
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'G':
      if ( atm->residue[1] == 'L' )
	switch ( atm->residue[2] )
	  {
	  case 'N':
	      atm->res = GLN;
	      break;
	  case 'U':
	    atm->res = GLU;
	      break;
	  case 'Y':
	    atm->res = GLY;
	      break;
	  default:
	    err_unknown_residue2 ( atm->residue, atm->residue_num);
	    break;
	  }
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'H':
      if ( strncmp ( atm->residue, "HIS", 3 ) == 0 )
	atm->res = HIS;
      else
	if ( strncmp ( atm->residue, "HOH", 3 ) == 0 )
	  atm->res = HOH;
	else
	  err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'I':
      if ( strncmp ( atm->residue, "ILE", 3 ) == 0 )
	atm->res = ILE;
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'L':
      switch ( atm->residue[1] )
	{
	case 'E':
	  atm->res = LEU;
	  break;
	case 'Y':
	  atm->res = LYS;
	  break;
	default:
	  err_unknown_residue2 ( atm->residue, atm->residue_num);
	  break;
	}
      break;
    case 'M':
      if ( strncmp ( atm->residue, "MET", 3 ) == 0 )
	atm->res = MET;
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'P':
      switch ( atm->residue[1] )
	{
	case 'H':
	  atm->res = PHE;
	  break;
	case 'R':
	  atm->res = PRO;
	  break;
	case 'C':           /* treat PCA (pyrrolidone carboxylic acid) like Proline */
	  atm->res = PCA;
	  break;
	default:
	  err_unknown_residue2 ( atm->residue, atm->residue_num);
	  break;
	}
      break;
    case 'S':
      if ( strncmp ( atm->residue, "SER", 3 ) == 0 )
	atm->res = SER;
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    case 'T':
      switch ( atm->residue[1] )
	{
	case 'H':
	  atm->res = THR;
	  break;
	case 'R':
	  atm->res = TRP;
	  break;
	case 'Y':
	  atm->res = TYR;
	  break;
	default:
	  err_unknown_residue2 ( atm->residue, atm->residue_num);
	  break;
	}
      break;
    case 'V':
      if ( strncmp ( atm->residue, "VAL", 3 ) == 0 )
	atm->res = VAL;
      else
	err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    default:
      err_unknown_residue2 ( atm->residue, atm->residue_num);
      break;
    }
}
	  
void  assign_type_and_radius ( atom_pt atm)
{
  atm->level = ALPHA;     /* default level */
  atm->rad = -99.9;
  atm->act = NOTHING;
  if ( ( strncmp ( atm->name, "N", 1 ) == 0 
	 && atm->res != PRO && atm->res != PCA )
       || strncmp ( atm->name, "O", 1 ) == 0 ) 
    { 
      if ( strncmp ( atm->name, "N", 1 ) == 0 ) 
	atm->act = DONOR;
      if ( strncmp ( atm->name, "O", 1 ) == 0 ) 
	atm->act = ACCEPTOR;
      if ( strncmp ( atm->name, "OG", 2 ) == 0 
	   || strncmp ( atm->name, "OH", 2 ) == 0 
	   || strncmp ( atm->name, "AD", 2 ) == 0
	   || strncmp ( atm->name, "AE", 2 ) == 0 ) 
	atm->act = DONEPTOR;
    }

#ifdef TRACE
  printf ("atom residue %d\n", atm->res );
#endif
  switch (atm->res) 
    {
    case ALA:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_ALAN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ALACA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ALAC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ALAO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_ALAOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_ALACB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case ACE:   
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ACEO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ACEC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"CH3") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ACECA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;     
    case ARG:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_ARGN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ARGCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ARGC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ARGO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_ARGOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_ARGCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_ARGCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_ARGCD;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"NE") == 0) {
	atm->level = EPSILON;
	atm->type = NE;
	atm->rad = RAD_ARGNE;
	atm->orbit = PL3;
	strcpy(atm->type_str,  "N.pl3");
      }
      if (strcmp(atm->name,"CZ") == 0) {
	atm->level = ZETA;
	atm->type = CZ;
	atm->rad = RAD_ARGCZ;
	atm->orbit = CAT;
	strcpy(atm->type_str,  "C.cat");
      }
      if (strcmp(atm->name,"NH1") == 0) {
	atm->level = ETA;
	atm->type = NH1;
	atm->rad = RAD_ARGNH1;
	atm->orbit = PL3;
	strcpy(atm->type_str,  "N.pl3");
      }
      if (strcmp(atm->name,"NH2") == 0) {
	atm->level = ETA;
	atm->type = NH2;
	atm->rad = RAD_ARGNH2;
	atm->orbit = PL3;
	strcpy(atm->type_str,  "N.pl3");
      }
      break;
    case ASN:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_ASNN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ASNCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ASNC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ASNO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_ASNOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_ASNCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_ASNCG;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"OD1") == 0) {
	atm->level = DELTA;
	atm->type = OD1;
	atm->rad = RAD_ASNOD1;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"ND2") == 0) {
	atm->level = DELTA;
	atm->type = ND2;
	atm->rad = RAD_ASNND2;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"AD1") == 0) {
	atm->level = DELTA;
	atm->type = AD1;
	atm->rad = RAD_ASNAD1;
	atm->orbit = AMBIG;
	strcpy(atm->type_str,  "Ambig");
      }
      if (strcmp(atm->name,"AD2") == 0) {
	atm->level = DELTA;
	atm->type = AD2;
	atm->rad = RAD_ASNAD2;
	atm->orbit = AMBIG;
	strcpy(atm->type_str,  "Ambig");
      }
      break;
    case ASP:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_ASPN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ASPCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ASPC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ASPO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_ASPOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_ASPCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_ASPCG;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"OD1") == 0) {
	atm->level = DELTA;
	atm->type = OD1;
	atm->rad = RAD_ASPOD1;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OD2") == 0) {
	atm->level = DELTA;
	atm->type = OD2;
	atm->rad = RAD_ASPOD2;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }     
      break;
    case CYS:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_CYSN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_CYSCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_CYSC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_CYSO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_CYSOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_CYSCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"SG") == 0) {
	atm->level = GAMMA;
	atm->type = SG;
	atm->rad = RAD_CYSSG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "S.3");
      }
      break;
    case GLN:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_GLNN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_GLNCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_GLNC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_GLNO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_GLNOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_GLNCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_GLNCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_GLNCD;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"OE1") == 0) {
	atm->level = EPSILON;
	atm->type = OE1;
	atm->rad = RAD_GLNOE1;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"NE2") == 0) {
	atm->level = EPSILON;
	atm->type = NE2;
	atm->rad = RAD_GLNNE2;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"AE1") == 0) {
	atm->level = EPSILON;
	atm->type = AE1;
	atm->rad = RAD_GLNAE2;
	atm->orbit = AMBIG;
	strcpy(atm->type_str,  "Ambig");
      }
      if (strcmp(atm->name,"AE2") == 0) {
	atm->level = EPSILON;
	atm->type = AE2;
	atm->rad = RAD_GLNAE2;
	atm->orbit = AMBIG;
	strcpy(atm->type_str,  "Ambig");
      }
      break;
    case GLU:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_GLUN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_GLUCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_GLUC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_GLUO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_GLUOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_GLUCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_GLUCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_GLUCD;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"OE1") == 0) {
	atm->level = EPSILON;
	atm->type = OE1;
	atm->rad = RAD_GLUOE1;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"OE2") == 0) {
	atm->level = EPSILON;
	atm->type = OE2;
	atm->rad = RAD_GLUOE2;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }      
      break;
    case GLY:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_GLYN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_GLYCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_GLYC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_GLYO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
       if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_GLYOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
     break;
    case HIS:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_HISN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_HISCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_HISC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_HISO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_HISOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_HISCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_HISCG;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"ND1") == 0) {
	atm->level = DELTA;
	atm->type = ND1;
	atm->rad = RAD_HISND1;
	atm->orbit = HISN;
	strcpy(atm->type_str,  "N.pl3");
      }
      if (strcmp(atm->name,"CD2") == 0) {
	atm->level = DELTA;
	atm->type = CD2;
	atm->rad = RAD_HISCD2;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"CE1") == 0) {
	atm->level = EPSILON;
	atm->type = CE1;
	atm->rad = RAD_HISCE1;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"NE2") == 0) {
	atm->level = EPSILON;
	atm->type = NE2;
	atm->rad = RAD_HISNE2;
	atm->orbit = HISN;
	strcpy(atm->type_str,  "N.pl3");
      }
      break;
    case ILE:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_ILEN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_ILECA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_ILEC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_ILEO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_ILEOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_ILECB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG1") == 0) {
	atm->level = GAMMA;
	atm->type = CG1;
	atm->rad = RAD_ILECG1;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG2") == 0) {
	atm->level = GAMMA;
	atm->type = CG2;
	atm->rad = RAD_ILECG2;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strncmp(atm->name,"CD",2) == 0) {
	atm->level = DELTA;
	atm->type = CD1;
	atm->rad = RAD_ILECD1;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case LEU:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_LEUN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_LEUCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_LEUC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_LEUO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_LEUOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_LEUCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_LEUCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD1") == 0) {
	atm->level = DELTA;
	atm->type = CD1;
	atm->rad = RAD_LEUCD1;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD2") == 0) {
	atm->level = DELTA;
	atm->type = CD2;
	atm->rad = RAD_LEUCD2;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case LYS:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_LYSN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_LYSCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_LYSC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_LYSO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_LYSOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
     if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_LYSCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_LYSCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_LYSCD;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CE") == 0) {
	atm->level = EPSILON;
	atm->type = CE;
	atm->rad = RAD_LYSCE;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"NZ") == 0) {
	atm->level = ZETA;
	atm->type = NZ;
	atm->rad = RAD_LYSNZ;
	atm->orbit = SP4;
	strcpy(atm->type_str,  "N.4");
      }
      break;
    case MET:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_METN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_METCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_METC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_METO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_METOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_METCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_METCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"SD") == 0) {
	atm->level = DELTA;
	atm->type = SD;
	atm->rad = RAD_METSD;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "S.3");
      }
      if (strcmp(atm->name,"CE") == 0) {
	atm->level = EPSILON;
	atm->type = CE;
	atm->rad = RAD_METCE;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case PCA:  
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_PRON;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_PROCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_PROC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_PROO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strncmp(atm->name,"OE", 2) == 0) {
	atm->level = EPSILON;
	atm->type = OE1;
	atm->rad = RAD_GLNOE1;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
     if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_PROCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_PROCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_PROCD;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case PHE:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_PHEN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_PHECA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_PHEC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_PHEO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_PHEOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_PHECB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_PHECG;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CD1") == 0) {
	atm->level = DELTA;
	atm->type = CD1;
	atm->rad = RAD_PHECD1;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CD2") == 0) {
	atm->level = DELTA;
	atm->type = CD2;
	atm->rad = RAD_PHECD2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CE1") == 0) {
	atm->level = EPSILON;
	atm->type = CE1;
	atm->rad = RAD_PHECE1;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CE2") == 0) {
	atm->level = EPSILON;
	atm->type = CE2;
	atm->rad = RAD_PHECE2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CZ") == 0) {
	atm->level = ZETA;
	atm->type = CZ;
	atm->rad = RAD_PHECZ;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      break;
    case PRO:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_PRON;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_PROCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_PROC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_PROO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_PROOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
     if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_PROCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_PROCG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CD") == 0) {
	atm->level = DELTA;
	atm->type = CD;
	atm->rad = RAD_PROCD;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case SER:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_SERN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_SERCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_SERC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_SERO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_SEROXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_SERCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"OG") == 0) {
	atm->level = GAMMA;
	atm->type = OG;
	atm->rad = RAD_SEROG;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "O.3");
      }
      break;
    case THR:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_THRN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_THRCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_THRC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_THRO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_THROXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_THRCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"OG1") == 0) {
	atm->level = GAMMA;
	atm->type = OG1;
	atm->rad = RAD_THROG1;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "O.3");
      }
      if (strcmp(atm->name,"CG2") == 0) {
	atm->level = GAMMA;
	atm->type = CG2;
	atm->rad = RAD_THRCG2;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      break;
    case TRP:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_TRPN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_TRPCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_TRPC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_TRPO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_TRPOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_TRPCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_TRPCG;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }      
      if (strcmp(atm->name,"CD1") == 0) {
	atm->level = DELTA;
	atm->type = CD1;
	atm->rad = RAD_TRPCD1;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }      
      if (strcmp(atm->name,"CD2") == 0) {
	atm->level = DELTA;
	atm->type = CD2;
	atm->rad = RAD_TRPCD2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      if (strcmp(atm->name,"NE1") == 0) {
	atm->level = EPSILON;
	atm->type = NE1;
	atm->rad = RAD_TRPNE1;
	atm->orbit = PL3;
	strcpy(atm->type_str,  "N.pl3");
      }      
      if (strcmp(atm->name,"CE2") == 0) {
	atm->level = EPSILON;
	atm->type = CE2;
	atm->rad = RAD_TRPCE2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      if (strcmp(atm->name,"CE3") == 0) {
	atm->level = EPSILON;
	atm->type = CE3;
	atm->rad = RAD_TRPCE3;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      if (strcmp(atm->name,"CZ2") == 0) {
	atm->level = ZETA;
	atm->type = CZ2;
	atm->rad = RAD_TRPCZ2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      if (strcmp(atm->name,"CZ3") == 0) {
	atm->level = ZETA;
	atm->type = CZ3;
	atm->rad = RAD_TRPCZ3;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      if (strcmp(atm->name,"CH2") == 0) {
	atm->level = ETA;
	atm->type = CH2;
	atm->rad = RAD_TRPCH2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      }      
      break;
    case TYR:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_TYRN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_TYRCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_TYRC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_TYRO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_TYROXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
      if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_TYRCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG") == 0) {
	atm->level = GAMMA;
	atm->type = CG;
	atm->rad = RAD_TYRCG;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CD1") == 0) {
	atm->level = DELTA;
	atm->type = CD1;
	atm->rad = RAD_TYRCD1;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CD2") == 0) {
	atm->level = DELTA;
	atm->type = CD2;
	atm->rad = RAD_TYRCD2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CE1") == 0) {
	atm->level = EPSILON;
	atm->type = CE1;
	atm->rad = RAD_TYRCE1;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CE2") == 0) {
	atm->level = EPSILON;
	atm->type = CE2;
	atm->rad = RAD_TYRCE2;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"CZ") == 0) {
	atm->level = ZETA;
	atm->type = CZ;
	atm->rad = RAD_TYRCZ;
	atm->orbit = AR;
	strcpy(atm->type_str,  "C.ar");
      } 
      if (strcmp(atm->name,"OH") == 0) {
	atm->level = ETA;
	atm->type = OH;
	atm->rad = RAD_TYROH;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "O.3");
      } 
      break;
    case VAL:
      if (strcmp(atm->name,"N") == 0) {
	atm->level = ALPHA;
	atm->type = N;
	atm->rad = RAD_VALN;
	atm->orbit = AM;
	strcpy(atm->type_str,  "N.am");
      }
      if (strcmp(atm->name,"CA") == 0) {
	atm->level = ALPHA;
	atm->type = CA;
	atm->rad = RAD_VALCA;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"C") == 0) {
	atm->level = ALPHA;
	atm->type = C;
	atm->rad = RAD_VALC;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "C.2");
      }
      if (strcmp(atm->name,"O") == 0) {
	atm->level = ALPHA;
	atm->type = O;
	atm->rad = RAD_VALO;
	atm->orbit = SP2;
	strcpy(atm->type_str,  "O.2");
      }
      if (strcmp(atm->name,"OXT") == 0) {
	atm->level = ALPHA;
	atm->type = OXT;
	atm->rad = RAD_VALOXT;
	atm->orbit = CO2;
	strcpy(atm->type_str,  "O.co2");
      }
     if (strcmp(atm->name,"CB") == 0) {
	atm->level = BETA;
	atm->type = CB;
	atm->rad = RAD_VALCB;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }
      if (strcmp(atm->name,"CG1") == 0) {
	atm->level = GAMMA;
	atm->type = CG1;
	atm->rad = RAD_VALCG1;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }      
      if (strcmp(atm->name,"CG2") == 0) {
	atm->level = GAMMA;
	atm->type = CG2;
	atm->rad = RAD_VALCG2;
	atm->orbit = SP3;
	strcpy(atm->type_str,  "C.3");
      }      
      break;
    default:
      atm->type = UNKNOWN;
      err_unknown_atom2 ( atm->name, atm->residue, atm->residue_num);
      break;
    }
  
  if (strncmp(atm->name,"H",1) == 0) {
      atm->rad = 1.0;
      atm->type = H;
      switch ( atm->name[1] )
	{
	case '\0':
	case 'A':
	  atm->level = ALPHA;
	  break;
	case 'B':
	  atm->level = BETA;
	  break;
	case 'G':
	  atm->level = GAMMA;
	  break;
	case 'D':
	  atm->level = DELTA;
	  break;
	case 'E':
	  atm->level = EPSILON;
	  break;
	case 'H':
	  atm->level = ETA;
	  break;
	case 'Z':
	  atm->level = ZETA;
	  break;
	default:
	  atm->type = UNKNOWN;
	  err_unknown_level2 ( atm->name, atm->residue, atm->residue_num);
	  break;
	}
  }
#ifdef TRACE
  printf ("atom name = %s [%7.3f,%7.3f,%7.3f]\n", atm->name,atm->pos[X],atm->pos[Y],atm->pos[Z]);
#endif
 if ( atm->rad < 0 )
      err_unknown_atom2 ( atm->name, atm->residue, atm->residue_num);
}


/* assign charges to certain 'O' and 'N' */
void  assign_charge ( atom_pt atm )
{
  if ( strcmp(atm->name, "OT1") == 0 || 
       strcmp(atm->name, "OT2") == 0 ||
       strcmp(atm->name, "OT") == 0 ||
       strcmp(atm->name, "OXT") == 0 ||
       ( atm->res == GLU && 
	 ( strcmp(atm->name, "OE1") == 0 || strcmp(atm->name, "OE2") == 0 ) 
	 ) ||
       ( atm->res == ASP && 
	 ( strcmp(atm->name, "OD1") == 0 ||strcmp(atm->name, "OD2" ) == 0 ) 
	 ) )
    {
      atm->charge = NEGATIVE_CHARGE;
      /*printf (" charge target atom %s\n", atm->name );*/
    }
  else
    if ( 
	 ( atm->res == ARG && 
	   ( strcmp(atm->name, "NE") == 0 || strcmp(atm->name, "NH1") == 0
	     || strcmp(atm->name, "NH2") == 0 ) 
	   ) ||
	 ( atm->res == LYS && strcmp(atm->name, "NZ") == 0 
	   ) ||
	 ( atm->res == HIS && 
	   ( strcmp(atm->name, "ND1") == 0 || strcmp(atm->name, "NE2") == 0 ) 
	   ) )
      {
	atm->charge = POSITIVE_CHARGE;
      }
}


void read_pdb(char *fname, atom_pt atoms, residue_pt residues,
              int what_to_read, int *number_of_atoms, int *number_of_residues)
{
  FILE           *in;
  atom_pt        atm;
  residue_pt     res;
  char           linebuf[MAX_PDB_LINELENGTH];
  char           prev_res[6];
  int            a, i, r; 
  int            is_atom, is_hydro;
  char tmp[6];

  in = fopen ( fname, "r" );
  if ( in == NULL ) {
      perror ( "ERROR  " );
      sprintf ( linebuf, "unable to open file %s", fname );
      err_panic2 ( "read_pdb", linebuf);
  }
  sprintf ( prev_res, "XXX" );
  a = i = r = 0;        /* counters on the number of entries for the vectors */
  atm = atoms;          /* local pointers on current vector elements */
  res = residues;
  res->start_interaction = UNKNOWN;
  res->start_atom = 0;

  /* read all the lines in the pdb-file */
  while ( fgets ( linebuf, sizeof linebuf, in ) != NULL
	  && ! strkeq ( ".", linebuf, 1 ) 
	  && ! strkeq ( "END\n", linebuf, 4 ) 
	  && ! strkeq ( "end\n", linebuf, 4 ) 
	  && ! strkeq ( "END ", linebuf, 4 ) 
	  /* dont trigger on ENDx extensions */
	  && ! strkeq ( "end ", linebuf, 4 ) )
    {
      /* grep only the 'ATOM'-lines */
      if ( what_to_read == ATOMS_ONLY )
	/* this is the 'normal' case, so get only ATOMs */
	is_atom = strkeq ( "ATOM", linebuf, 4 )
	  || strkeq ( "atom", linebuf, 4 );
      else
	/* in this case also read the HETATMs, but no waters */
	{
	  is_atom = strkeq ( "ATOM", linebuf, 4 )
	    || strkeq ( "atom", linebuf, 4 )
	    || ( ( strkeq ( "HETATM", linebuf, 6 )
		   || strkeq ( "hetatm", linebuf, 6 ) )
		 && strncmp ( "HOH", (char *) linebuf + 17, 3 ) != 0 );	  
	}

      is_hydro = ( toupper ( *column ( 14 ) ) == 68       /* 'D' or 'd' */
		   || toupper ( *column ( 14 ) ) == 72    /* 'H' or 'h' */
		   || toupper ( *column ( 14 ) ) == 81 ); /* 'Q' or 'q' */

      /* take care of only ATOM that no hydrogens and other junk is read */
      if ( is_atom 
	   && !is_hydro 
	   && 1 == sscanf ( column ( 14 ), " %3s", atm->name ) )
	{
	  atm->pos[X] = (float) atof ( column ( 31 ) );
	  atm->pos[Y] = (float) atof ( column ( 39 ) );
	  atm->pos[Z] = (float) atof ( column ( 47 ) );
#ifdef TRACE
	  if ( a > 2780 && a < 2800 )
	  printf ("[%7.3f,%7.3f,%7.3f]\n",atm->pos[X], atm->pos[Y], atm->pos[Z] );
#endif
	  /* get residue type, number, chainID, altloc, and insetion code */
	  strncpy ( atm->residue, column ( 18 ), 3 );
	  atm->residue[3] = '\0';
	  strncpy ( &atm->alt_location, column ( 17 ), 1);
	  strncpy ( &atm->chain_id, column ( 22 ), 1);
	  strncpy ( &atm->insertion_code, column( 27 ), 1);
	  /* toneroma 20FEB07  
           * changes made so that insertion codes work properly*/
	  strncpy ( atm->residue_num, column( 23 ), 5);
          atm->residue_num[5] = 0;
          strncpy(tmp, column(7), 5);
          tmp[5] = 0;
          atm->atom_number = strtol(tmp, (char**) NULL, 10);

	  if(strkeq ( "ATOM", linebuf, 4 ) || strkeq ( "atom", linebuf, 4 )){ 
            /* assign_residue_type to current atom */
	    assign_residue_type ( atm);
	    /* assign types and radii to current atom */
	    assign_type_and_radius ( atm);
	    atm->hydro = hydro[atm->res][atm->type];

            if(atm->type == OXT && atm->hydro == 0) atm->hydro = 530;

	  }else{
	      /* this marks that this atom is a HETATM */
	      atm->level = UNKNOWN;
	      atm->res = HETATM;
	      atm->act = NOTHING;
      
	      switch ( (int) atm->name[0] )
		{
		case 67:
		  atm->type = C;
		  atm->hydro = 80;
		  strcpy(atm->type_str,  "C.het");
		  break;
		case 69:
		  atm->type = FE;
		  atm->hydro = 317;
		  strcpy(atm->type_str,  "Fe");
		  break;
		case 78:
		  atm->type = N;
		  atm->hydro = 350;
		  strcpy(atm->type_str,  "N.het");
		  break;
		case 79:
		  atm->type = O;
		  atm->hydro = 530;
		  strcpy(atm->type_str,  "O.het");
		  break;
		case 80:
		  atm->type = P;
		  atm->hydro = 317;
		  strcpy(atm->type_str,  "P.het");
		  break;
		case 83:
		  atm->type = S;
		  atm->hydro = 80;
		  strcpy(atm->type_str,  "S.het");
		  break;
		default:
		  atm->type = UNKNOWN;
		  atm->hydro = 317;
		  strcpy(atm->type_str,  "HET");
		  break;
		}

	      /*****************************************************************
		handle the metal HETATM to count the metal hbond in scoring

                   metal_class  metal_name                    metal_hbond_distance
                    1            Ca, Na, K                       2.9 angstroms
                    2            Co, Cu, Fe, Mg, Mn, Ni, Zn      2.6 angstroms

                		       -- added by Litian He  12/02/2003  
	      ******************************************************************/

	      /* "CA" is a special case, both atom name and residue name should be "CA"  */
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'A' && 
		   *column ( 15 ) == ' ' && *column ( 18 ) == ' ' && 
		   *column ( 19 ) == 'C' && *column ( 20 ) == 'A' )
		 {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = CA_M;
		      atm->act = METAL_1;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Ca");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atm->name, linebuf);
#endif
		  }
	      /* we only need to check the atom name for all other metals */
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'O' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = CO;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Co");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'C' && *column ( 14 ) == 'U' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = CU;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Cu");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'F' && *column ( 14 ) == 'E' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = FE;
		      atm->act = METAL_2;
		      atm->hydro = 317;  
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Fe");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'K' && *column ( 14 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = K;
		      atm->act = METAL_1;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "K");
		      
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'M' && *column ( 14 ) == 'G' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = MG;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Mg");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'M' && *column ( 14 ) == 'N' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = MN;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Mn");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'N' && *column ( 14 ) == 'A' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = NA;
		      atm->act = METAL_1;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Na");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_1 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'N' && *column ( 14 ) == 'I' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = NI;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Ni");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }
	      if ( *column ( 13 ) == 'Z' && *column ( 14 ) == 'N' && *column ( 15 ) == ' ' )
		  {
		      sscanf ( column ( 13 ), " %3s", atm->name );
		      atm->type = ZN;
		      atm->act = METAL_2;
		      atm->hydro = 635;
		      atm->level = HETATM;
		      strcpy(atm->type_str,  "Zn");
#ifdef TRACE_METAL
		      printf("\"%s\" - METAL_2 : %s", atm->name, linebuf);
#endif
		  }

	      /******************* end of metal check ********************/

	      if ( atm->name[0] == 'O' || atm->name[0] == 'N' )
		{
		  if ( strncmp ( atm->name + 1, "DD", 2 ) == 0 )
		    atm->act = DONOR;
		  if ( strncmp ( atm->name + 1, "1D", 2 ) == 0 )
		    atm->act = DONOR;
		  if ( strncmp ( atm->name + 1, "2D", 2 ) == 0 )
		    atm->act = DONOR;
		  if ( strncmp ( atm->name + 1, "AA", 2 ) == 0 )
		    atm->act = ACCEPTOR;
		  if ( strncmp ( atm->name + 1, "NN", 2 ) == 0 )
		    atm->act = DONEPTOR;
		}
	      atm->rad = radius ( atm->type, UNKNOWN );
	    }

	  atm->charge = 0.0;
	  
	  assign_charge ( atm );

	  /* is this a new residue? (0 = No) */
	  if ( strncmp ( atm->residue_num, prev_res, 5 ) != 0 )
	    {
	      if ( prev_res[0] != 'X' )
		{
		  res[r].number_of_atoms = a - res[r].start_atom;
#ifdef TRACE
		  printf("finish a residue %d: %d - %d, total %d \n", 
			 r, res[r].start_atom, a-1, res[r].number_of_atoms);
		  /*printf("%s", linebuf);*/
		  printf("a=%d\n", a);
#endif
		  r++;
		  res[r].start_interaction = UNKNOWN;
		}
#ifdef TRACE
	      if ( r < 2 ) 
		printf ("residue %d: atom %d [%7.3f, %7.3f, %7.3f]\n",
			r, a, atm->pos[X], atm->pos[Y], atm->pos[Z] );
#endif
	      res[r].start_atom = a;
	      res[r].start_interaction = UNKNOWN;
	      res[r].type = atm->res;
	      strcpy ( res->name, atm->residue );
	      strcpy ( res->num, atm->residue_num );
	      strcpy ( prev_res, atm->residue_num );
	    }
	  if ( atm->level == UNKNOWN && atm->act == DONOR && 
	       atm->res != HETATM )
	    /* sorry, but we cannot assign hydrogens to this one later */
	    atm->act = NOTHING;
	  atm->residue_index = r;

	  if ( atm->type != UNKNOWN )
	    {
	      a++;
	      atm++;
	      /*printf("a=%d\n", a);*/
	    }
	  
	  if ( a >= MAX_PDB_ATOMS ) 
	    err_panic2 ( "read_pdb", "more than MAX_PDB_ATOMS atoms");
	  if ( r >= MAX_PDB_RESIDUES ) 
	    err_panic2 ( "read_pdb", "more than MAX_PDB_RESIDUES residues");
	}
    }
  fclose ( in );
  res[r].number_of_atoms = a - res[r].start_atom;
  if ( a == 0 ) a = 1;
  *number_of_atoms = a;
  *number_of_residues = r + 1;
  if ( res[r].start_interaction == UNKNOWN )
    {
      /* the last residue doesn't have an interaction point, thus skip it
	 during the screening, since this causes problems otherwise  e.g.
	 for 1lkk_a */
      (*number_of_residues)--;
      r--;
    }
  r ++;
  res[r].start_interaction = i;  /* for the last peptide in the chain, needed
				  for the global loop */


#ifdef TRACE

  printf ( "CHECK: rad file read\n");
  for ( i = 0; i < 200; i++)
    if(atoms[i].insertion_code != ' ')
      {
	printf ( "%4s %c\n", atoms[i].residue_num, atoms[i].insertion_code);
	printf ( "%d\n", strlen(atoms[i].residue_num));
      }
  for ( i = 0; i < 3; i++)
    printf ( "   %3s %3s %c %4s%c - %3d - %8.3f, %8.3f, %8.3f\n", 
	     atoms[i].residue, atoms[i].name, atoms[i].chain_id,
	     atoms[i].residue_num, atoms[i].insertion_code,
	     atoms[i].hydro,
	     atoms[i].pos[X],atoms[i].pos[Y],atoms[i].pos[Z]);
/*************************************/
#endif

#ifdef TRACE
  for ( i = 0; i < *number_of_atoms; i++ )
    printf ( "%4d: %3s-%3s-%4s-%3d-%3d-%c-%c-%c-%8.3f-%8.3f-%8.3f\n", 
	     i,atoms[i].residue, atoms[i].name, atoms[i].residue_num,
	     atoms[i].hydro, atoms[i].type, atoms[i].alt_location,
	     atoms[i].chain_id, atoms[i].insertion_code,
	     atoms[i].pos[X],atoms[i].pos[Y],atoms[i].pos[Z]);
	       
#endif

}



