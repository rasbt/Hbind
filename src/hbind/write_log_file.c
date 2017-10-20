#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/times.h>
#include <limits.h>
#include <unistd.h>
#include <netdb.h>
#include <sys/param.h>
#include "defs.h"
#include "types.h"
#include "err_handle.h"
#include "hbind_itimer.h"

void write_log_file(global_data_pt global, char *filename )
{
  FILE       *fp;
  char hname[MAXHOSTNAMELEN];
  char tbuf[80];
  double real, virt, prof;

  fp = fopen ( filename, "w");
  if ( fp == NULL ) {
      perror ( "ERROR  ");
      err_panic2 ( "write_log_file", "file open failed");
  }
  fprintf(fp, "\n********************************\n");
  fprintf(fp, "   HBIND screening statistics\n");
  fprintf(fp, "********************************\n\n");
  fprintf(fp, "target                   : %s\n", global->protein);
  fprintf(fp, "template                 : %s\n", global->template);
  fprintf(fp, "database                 : %s\n", global->database);
  fprintf(fp, "compounds                : %d\n\n", 
	    global->number_of_screened_compounds);
  fprintf(fp, "parameters\n");
  fprintf(fp, "----------\n");
  fprintf(fp, "DME threshold            : %3.1f\n", global->dmetol);
  fprintf(fp, "RMS threshold            : %3.1f\n", global->rmstol);
  fprintf(fp, "max anchor translation   : %3.1f\n", global->anchor_translation);
  fprintf(fp, "anchor overlap           : %3.1f\n", global->anchor_overlap);
  fprintf(fp, "side-chain overlap       : %3.1f\n", global->side_chain_overlap);
  fprintf(fp, "intra overlap            : %3.1f\n", global->intra_overlap);
  fprintf(fp, "intermediate overlap     : %3.1f\n",
	    global->intermediate_overlap);
  fprintf(fp, "tolerated overlap        : %3.1f\n",
	    global->finally_tolerated_overlap);
  fprintf(fp, "tolerated max bump       : %3.1f\n",
	    global->finally_tolerated_max_bump);
  fprintf(fp, "score cutoff             : %3.1f\n\n", global->score_cutoff);
  fprintf(fp, "hashing\n");
  fprintf(fp, "-------\n");
  fprintf(fp, "total template triangles : %d\n",
	    global->number_of_template_triangles);
  fprintf(fp, "non-empty buckets        : %d\n",
	    global->number_of_hash_buckets);
  fprintf(fp, "total triangle pointers  : %d\n",
	    global->number_of_triangle_pointers);
  fprintf(fp, "average triangle pointers: %d\n",
	  global->number_of_triangle_pointers / global->number_of_hash_buckets);
  fprintf(fp, "maximal triangle pointers: %d\n\n",
	    global->max_number_of_triangles);
  fprintf(fp, "filter effectiveness\n");
  fprintf(fp, "--------------------\n");
  fprintf(fp, "triangle match           : %d\n", 
	    global->filter_counter[TRIANGLE_MATCH]);
  fprintf(fp, "triangle DME             : %d\n", global->filter_counter[DME]);
  fprintf(fp, "RMS deviation            : %d\n", 
	    global->filter_counter[RMS_DEVI]);
  fprintf(fp, "anchor unbump            : %d\n", 
	    global->filter_counter[BUMP_ANCHOR]);
  fprintf(fp, "side-chain unbump        : %d\n", 
	    global->filter_counter[BUMP_SIDE_CHAIN]);
  fprintf(fp, "scoring                  : %d\n\n", 
	    global->filter_counter[SCORING]);
  fprintf(fp, "potential ligands        : %d\n\n", 
	    global->number_of_potential_ligands);
  fprintf(fp, "best affiscore           : %0.3f (name: %s, significance: %0.3f)\n",
	  global->best_affiscore_so_far, global->best_affi_name, global->affiscore_significance);
  fprintf(fp, "affiscore mean           : %0.3f\n",
	  global->affiscore_mean);
  fprintf(fp, "affiscore stdd           : %0.3f\n\n",
	  global->affiscore_stdd);
  fprintf(fp, "timing\n------\n");	    
  gethostname(hname, MAXHOSTNAMELEN);
  fprintf(fp, "host                     : %s\n", hname);
  hbind_get_local_time(tbuf, 80);
  fprintf(fp, "finished                 : %s", tbuf);
  hbind_itimer_get(&real, &virt, &prof);
  fprintf(fp, "Wall Clock Time (seconds): %8.2f", real);
  fclose ( fp);
}

