#include "defs.h"
#include "types.h"
#include "dist_fun.h"
#include "hbond_check.h"

int  intra_target_hbonds ( global_data_pt  global )
{
  atom_pt target;
  float   angle;
  int     donor,
          acceptor,
          sum,
          number_of_hbonds;
  int     i, j;
  float   hydrogen_position[3];

  target = global->target_atoms;
  number_of_hbonds = 0;
  for ( i = 0; i < global->number_of_target_atoms; i++ )
    for ( j = i + 1; j < global->number_of_target_atoms; j++ )
      if ( target[i].act != NOTHING 
	   && target[j].act != NOTHING
	   && dist_fun ( target[i].pos, target[j].pos ) <= MAX_HBOND_LENGTH )

	{
	  sum = target[i].act + target[j].act;
	  if ( sum == ACCEPTOR_DONOR
	       || sum == DONOR_DONEPTOR
	       || sum == ACCEPTOR_DONEPTOR
	       || sum == DONEPTOR_DONEPTOR ) 
	    /* this is an intramolecular hydrogen bond in the target protein */
	    {
	      if ( target[i].act == DONOR 
		   || ( target[i].act == DONEPTOR
			&& target[j].act != DONOR ) )
		{
		  donor = i;
		  acceptor = j;
		}
	      else
		{
		  donor = j;
		  acceptor = i;
		}		 
	      /* now compute the donor-H-acceptor angle */
	      angle = 
		compute_target_hydrogen_angle ( target,
						global->target_residues,
						donor,
						target[acceptor].pos,
						hydrogen_position,
						global);
	      if ( angle >= MIN_HYDROGEN_ANGLE
		   && angle <= MAX_HYDROGEN_ANGLE )
		number_of_hbonds++;
	    }
	}
  return number_of_hbonds;
}
      
int  target_to_water_hbonds ( global_data_pt  global )
{ 
  atom_pt target,
          waters;
  float   angle;
  float   number_of_hbonds;
  int     i, j;
  float   hydrogen_position[3];

  if ( global->number_of_waters == 0 )
    /* no binding-site waters */
    return 0;
  target = global->target_atoms;
  waters = global->waters;
  number_of_hbonds = 0.0;
  for ( i = 0; i < global->number_of_target_atoms; i++ )
    if ( target[i].act != NOTHING )
      for ( j = 0; j < global->number_of_waters; j++ )
	if ( waters[j].state == CONSERVED
	     && dist_fun ( target[i].pos, waters[j].pos ) < MAX_HBOND_LENGTH )
	  {
	    if ( target[i].act == DONOR )
	      /* compute the donor-H-acceptor angle */
	      {
		angle = 
		  compute_target_hydrogen_angle ( target,
						  global->target_residues,
						  i,
						  waters[j].pos,
						  hydrogen_position,
						  global);
		if ( angle >= MIN_HYDROGEN_ANGLE
		     && angle <= MAX_HYDROGEN_ANGLE )
		  number_of_hbonds += waters[j].prediction;
	      }
	    else
	      /* target is either DONEPTOR or ACCEPTOR, in both cases
		 we already have a hydrogen bond, since we only consider
		 the distance between this atom and the water */
	      number_of_hbonds += waters[j].prediction;
	  }
  return (int) number_of_hbonds;
}
 
