#include <stdio.h>
#include <string.h>
#include "types.h"
#include "defs.h"
#include "count_flexible_bonds.h"
#include "find_cycles.h"
#include "err_handle.h"


int  count_flexible_ligand_bonds ( molecule_pt     molecule,
                                   short           ligand_flag[MAX_LIGAND_ATOMS] )
{
    atom_pt    atoms;
    int number_of_flexible_interfacial_ligand_bonds;
    int i, j,
        atom1, atom2,
        only_hydrogen_neighbors;
    
    
    atoms = molecule->atoms;
    number_of_flexible_interfacial_ligand_bonds = 0;
    
    for ( i = 0; i < molecule->number_of_bonds; i++ )
    {
        atom1 = molecule->bonds[i].atom1;
        atom2 = molecule->bonds[i].atom2;
        

        if ( molecule->bonds[i].type == FLEX &&
            ( ligand_flag[atom1] != INITIAL || ligand_flag[atom2] != INITIAL ) &&
            ( molecule->number_of_neighbors[atom1] > 1 &&
                molecule->number_of_neighbors[atom2] > 1 ) )
        {

            only_hydrogen_neighbors = 0;
            
            if ( atoms[atom1].name[0] == 'C' )
            {
                only_hydrogen_neighbors = 1;
                for ( j = 0; j < molecule->number_of_neighbors[atom1]; j ++ )
                {
                    
                    if ( molecule->neighbors[atom1][j] != atom2 &&
                        atoms[molecule->neighbors[atom1][j]].type != H )
                    {
                        only_hydrogen_neighbors = 0;  
                        break;
                    }
                }
            }

            if ( only_hydrogen_neighbors == 0 &&
                atoms[atom2].name[0] == 'C' )
            {
                only_hydrogen_neighbors = 1;
                for ( j = 0; j < molecule->number_of_neighbors[atom2]; j ++ )
                {
                    
                    if ( molecule->neighbors[atom2][j] != atom1 &&
                        atoms[molecule->neighbors[atom2][j]].type != H )
                    {
                        only_hydrogen_neighbors = 0;  
                        break;
                    }
                }
            }
            
            if ( only_hydrogen_neighbors == 0 )
            {   
                number_of_flexible_interfacial_ligand_bonds ++;
       
            }
        }
#ifdef TRACE3
        else if ( molecule->bonds[i].type == FLEX &&
            ( ligand_flag[atom1] != INITIAL || ligand_flag[atom2] != INITIAL ) )
        {
            printf("%d has %d neighbors, %d has %d neighbors  ",
                atom1, molecule->number_of_neighbors[atom1],
                atom2, molecule->number_of_neighbors[atom2] );
            printf("non-interfacial bond %d: %d %s [%7.2f, %7.2f, %7.2f] - %d %s [%7.2f, %7.2f, %7.2f]\n", i,
                molecule->bonds[i].atom1,
                atoms[molecule->bonds[i].atom1].name,
                atoms[molecule->bonds[i].atom1].pos[X],
                atoms[molecule->bonds[i].atom1].pos[Y],
                atoms[molecule->bonds[i].atom1].pos[Z],
                molecule->bonds[i].atom2,
                atoms[molecule->bonds[i].atom2].name,
                atoms[molecule->bonds[i].atom2].pos[X],
                atoms[molecule->bonds[i].atom2].pos[Y],
                atoms[molecule->bonds[i].atom2].pos[Z]);
        }
#endif        
    }
#ifdef TRACE
    printf("number_of_flexible_interfacial_ligand_bonds = %d\n",
        number_of_flexible_interfacial_ligand_bonds );
#endif    
    return (number_of_flexible_interfacial_ligand_bonds);
}

int  count_flexible_target_bonds (
    global_data_pt       global,
    short           target_flag[MAX_PDB_ATOMS] )
{
    int number_of_flexible_interfacial_target_bonds;
    int i,
      current_residue;
    short atom_CA, atom_CB, atom_CG, atom_CG1, atom_CD, atom_CE,
        atom_CZ, atom_CH, atom_NE, atom_NZ, atom_NH1, atom_NH2,
        atom_ND2, atom_NE2, atom_SD, atom_OG, atom_OH;
    atom_pt    atoms;
    
    
#ifdef TRACE2
    printf("count_flexible_target_bonds:\n");
#endif

    atoms = global->target_atoms;
    number_of_flexible_interfacial_target_bonds = 0;
    current_residue = 0;
    atom_CA = 0;
    atom_CB = 0;
    atom_CG = 0;
    atom_CG1 = 0;
    atom_CD = 0;
    atom_CE = 0;
    atom_CZ = 0;
    atom_CH = 0;
    atom_NE = 0;
    atom_NZ = 0;
    atom_NH1 = 0;
    atom_NH2 = 0;
    atom_ND2 = 0;
    atom_NE2 = 0;
    atom_SD = 0;
    atom_OG = 0;
    atom_OH = 0;
    
    for ( i = 0; i < global->number_of_target_atoms; i++ )
    {
#ifdef TRACE2
        printf("residue %d - atom %d: %s [%7.2f, %7.2f, %7.2f] \n",
            atoms[i].res, i, atoms[i].name, 
            atoms[i].pos[X], atoms[i].pos[Y], atoms[i].pos[Z] );
#endif
        if ( current_residue != atoms[i].res )
        {
#ifdef TRACE2
            printf(" current_residue = %d, atoms[i].res = %d\n",
                current_residue, atoms[i].res );
#endif

              switch (current_residue) 
              {
                case ALA:                  
                  break;
                case ARG:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_CD == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CD == 1 || atom_NE == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_NE == 1 || atom_CZ == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CZ == 1 || atom_NH1 == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CZ == 1 || atom_NH2 == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case ASN:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_ND2 == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case ASP:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case CYS:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case GLU:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_CD == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case GLN:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_CD == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CD == 1 || atom_NE2 == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case GLY:
                  break;
                case HIS:
#ifdef TRACE2
            printf("residue %d: atom_CA = %d, atom_CB = %d, atom_CG = %d, atom_CG1 = %d, atom_CD = %d, atom_CE = %d, \
        atom_CZ = %d, atom_CH = %d, atom_NE = %d, atom_NZ = %d, atom_NH1 = %d, atom_NH2 = %d, \
        atom_ND2 = %d, atom_NE2 = %d, atom_SD = %d, atom_OG = %d, atom_OH = %d \n",
                current_residue, atom_CA, atom_CB, atom_CG, atom_CG1, atom_CD, atom_CE,
                atom_CZ, atom_CH, atom_NE, atom_NZ, atom_NH1,atom_NH2,
                atom_ND2, atom_NE2, atom_SD, atom_OG, atom_OH );
#endif
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif
                  break;
                case ILE:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG1 == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case LEU:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case LYS:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_CD == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CD == 1 || atom_CE == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CE == 1 || atom_NZ == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif
                  break;
                case MET:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CG == 1 || atom_SD == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case PHE:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case PRO:
                  break;
                case SER:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_OG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case THR:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_OG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
                case TRP:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;                  
                case TYR:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CB == 1 || atom_CG == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
                  if ( atom_CZ == 1 || atom_OH == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif
                  break;
                case VAL:
                  if ( atom_CA == 1 || atom_CB == 1 )
                      number_of_flexible_interfacial_target_bonds ++;
#ifdef TRACE2
                  printf("residue %d: number_of_flexible_interfacial_target_bonds = %d\n",
                      current_residue, number_of_flexible_interfacial_target_bonds);
#endif                  
                  break;
              }
              
            /* re-initialization */
            current_residue = atoms[i].res;
            atom_CA = 0;
            atom_CB = 0;
            atom_CG = 0;
            atom_CG1 = 0;
            atom_CD = 0;
            atom_CE = 0;
            atom_CZ = 0;
            atom_CH = 0;
            atom_NE = 0;
            atom_NZ = 0;
            atom_NH1 = 0;
            atom_NH2 = 0;
            atom_ND2 = 0;
            atom_NE2 = 0;
            atom_SD = 0;
            atom_OG = 0;
            atom_OH = 0;
        }
        else
        {
#ifdef TRACE2          
            if ( target_flag[i] != INITIAL )
                printf("residue %d - interfacial atom %d: %s [%7.2f, %7.2f, %7.2f] \n",
                    atoms[i].res, i, atoms[i].name, 
                    atoms[i].pos[X], atoms[i].pos[Y], atoms[i].pos[Z] );
#endif      
              switch (atoms[i].res) 
              {
                case ALA:
                  break;
                case ARG:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "CD" ) == 0 && target_flag[i] != INITIAL )
                      atom_CD = 1;
                  if ( strcmp ( atoms[i].name, "NE" ) == 0 && target_flag[i] != INITIAL )
                      atom_NE = 1;
                  if ( strcmp ( atoms[i].name, "CZ" ) == 0 && target_flag[i] != INITIAL )
                      atom_CZ = 1;
                  if ( strcmp ( atoms[i].name, "NH1" ) == 0 && target_flag[i] != INITIAL )
                      atom_NH1 = 1;
                  if ( strcmp ( atoms[i].name, "NH2" ) == 0 && target_flag[i] != INITIAL )
                      atom_NH2 = 1;                  
                  break;
                case ASN:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "ND2" ) == 0 && target_flag[i] != INITIAL )
                      atom_ND2 = 1;                 
                  break;
                case ASP:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  break;
                case CYS:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  break;
                case GLU:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "CD" ) == 0 && target_flag[i] != INITIAL )
                      atom_CD = 1;
                  break;
                case GLN:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "CD" ) == 0 && target_flag[i] != INITIAL )
                      atom_CD = 1;
                  if ( strcmp ( atoms[i].name, "NE2" ) == 0 && target_flag[i] != INITIAL )
                      atom_NE2 = 1;
                  break;
                case GLY:
                  break;
                case HIS:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;                  
                  break;
                case ILE:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG1" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG1 = 1; 
                  break;
                case LEU:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1; 
                  break;
                case LYS:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "CD" ) == 0 && target_flag[i] != INITIAL )
                      atom_CD = 1;
                  if ( strcmp ( atoms[i].name, "CE" ) == 0 && target_flag[i] != INITIAL )
                      atom_CE = 1;
                  if ( strcmp ( atoms[i].name, "NZ" ) == 0 && target_flag[i] != INITIAL )
                      atom_NZ = 1;
                  break;
                case MET:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "SD" ) == 0 && target_flag[i] != INITIAL )
                      atom_SD = 1;
                  break;
                case PHE:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  break;
                case PRO:
                  break;
                case SER:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "OG" ) == 0 && target_flag[i] != INITIAL )
                      atom_OG = 1;
                  break;
                case THR:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "OG" ) == 0 && target_flag[i] != INITIAL )
                      atom_OG = 1;
                  break;
                case TRP:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  break;                  
                case TYR:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  if ( strcmp ( atoms[i].name, "CG" ) == 0 && target_flag[i] != INITIAL )
                      atom_CG = 1;
                  if ( strcmp ( atoms[i].name, "CZ" ) == 0 && target_flag[i] != INITIAL )
                      atom_CZ = 1;
                  if ( strcmp ( atoms[i].name, "OH" ) == 0 && target_flag[i] != INITIAL )
                      atom_OH = 1;
                  break;
                case VAL:
                  if ( strcmp ( atoms[i].name, "CA" ) == 0 && target_flag[i] != INITIAL )
                      atom_CA = 1;
                  if ( strcmp ( atoms[i].name, "CB" ) == 0 && target_flag[i] != INITIAL )
                      atom_CB = 1;
                  break;
              }
        } 
    }

    switch (current_residue) 
    {
      case ALA:
        break;
      case ARG:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_CD == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CD == 1 || atom_NE == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_NE == 1 || atom_CZ == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CZ == 1 || atom_NH1 == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CZ == 1 || atom_NH2 == 1 )
            number_of_flexible_interfacial_target_bonds ++;                  
        break;
      case ASN:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_ND2 == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case ASP:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case CYS:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case GLU:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_CD == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case GLN:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_CD == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CD == 1 || atom_NE2 == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case GLY:
        break;
      case HIS:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;                  
        break;
      case ILE:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG1 == 1 )
            number_of_flexible_interfacial_target_bonds ++; 
        break;
      case LEU:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case LYS:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_CD == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CD == 1 || atom_CE == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CE == 1 || atom_NZ == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case MET:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CG == 1 || atom_SD == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case PHE:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case PRO:
        break;
      case SER:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_OG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case THR:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_OG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case TRP:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;                  
      case TYR:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CB == 1 || atom_CG == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        if ( atom_CZ == 1 || atom_OH == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
      case VAL:
        if ( atom_CA == 1 || atom_CB == 1 )
            number_of_flexible_interfacial_target_bonds ++;
        break;
    }

#ifdef TRACE2
    printf("number_of_flexible_interfacial_target_bonds = %d\n",
        number_of_flexible_interfacial_target_bonds );
#endif
    
    return (number_of_flexible_interfacial_target_bonds);

}
