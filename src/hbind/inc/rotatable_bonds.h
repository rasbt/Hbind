static int  number_of_rotatable_bonds[20] = 
{  0,  4,  2,  2,  1, 
   3,  3,  0,  2,  2, 
   2,  4,  3,  2,  0,  /* proline (last entry in this line) had originally 2 */
   1,  1,  2,  2,  1  };

static int  rotatable_bonds[20][6] = 
{
  /* ALA */
  { CA },
  /* ARG */
  { CA, CB, CG, CD, NE, CZ },
  /* ASN */
  { CA, CB, CG, OD1 },
  /* ASP */
  { CA, CB, CG, OD1 },
  /* CYS */
  { CA, CB, SG },
  /* GLN */
  { CA, CB, CG, CD, OE1 },
  /* GLU */
  { CA, CB, CG, CD, OE1 },
  /* GLY */
  { CA },
  /* HIS */
  { CA, CB, CG, ND1 },
  /* ILE */
  { CA, CB, CG1, CD1 }, 
  /* LEU */
  { CA, CB, CG, CD1 },
  /* LYS */
  { CA, CB, CG, CD, CE, NZ },
  /* MET */
  { CA, CB, CG, SD, CE },
  /* PHE */
  { CA, CB, CG, CD1 },
  /* PRO */
  { CA }, /* , CB, CG, CD }, in the rotamer library rotamers for proline
	     have been given, but simple rotating those bonds does not work */
  /* SER */
  { CA, CB, OG },
  /* THR */
  { CA, CB, OG1 },
  /* TRP */
  { CA, CB, CG, CD1 },
  /* TYR */
  { CA, CB, CG, CD1 },
  /* VAL */
  { CA, CB, CG1 }
};
