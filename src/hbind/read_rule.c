#include <stdio.h>
#include "read_rule.h"
#include "err_handle.h"
#include "assign_type.h"
#include "mymalloc.h"

char *get_rule_item(char *line, rule_list_pt rule)
{
  char  *rulestart;
  char   err_msg[FILENAME_MAX];

  rule->required = rule->prohibited = NULL;
  rulestart = line;
  /* move to the first item */
  while ( *line == ' ' || *line == '[' ||  *line == '(' ) line++;
  /* number of required/prohibited atoms */
  if ( *line >= '1' && *line <= '9' ) {
    rule->number = (int) *line - 48;
    line++;
  }
  /* no number, so at least one atom of the corresponding type */
  else rule->number = 1;
  while ( *line == ' ' || *line == '\t'  ) line++;
  /* line points to the definition of the atom bonded to the current bond */
  assign_type_and_orbit ( line, &rule->type, &rule->orbit);
  /* skip atom field and appended spaces */
  while ( *line != ' ' && *line != '\t' ) line++;
  while ( *line == ' ' || *line == '\t' ) line++;
  /* rule for a prohibited atom connected to the current atom */
  if ( *line == '[' ) {
    rule->prohibited = (rule_list_pt) mymalloc ( sizeof (rule_list_t) );
    line = get_rule_item ( line, rule->prohibited);
    if ( *line != ']' ) {
      fprintf ( stderr, "something is wrong with rule %s\n", rulestart );
      sprintf ( err_msg, "something is wrong with rule %s\n", rulestart );
      err_print (err_msg);
    }
  /* rule for a required atom, else, because only one rule grepped 
     per time */
  }else if ( *line == '(' ){
    rule->required = (rule_list_pt) mymalloc ( sizeof (rule_list_t) );
    line = get_rule_item ( line, rule->required);
    if ( *line != ')' ) {
      sprintf ( err_msg, "something is wrong with rule %s\n", rulestart );
      fprintf ( stderr, err_msg);
      err_print(err_msg);
    }
  }
  /* return the rest of the line */
  return line;
}
  
  
void  parse_definition_line ( char    *line, rule_pt rule)
{
  char  *field;
  int   p_counter, 
        r_counter;
  int   i;

  p_counter = r_counter = 0;
  /* jump to the beginning of the second entry in the line */
  while ( *line != ' ' && *line != '\t' )
    line++;
  /* extract first field (atom definition) in line */
  while ( *line == ' ' || *line == '\t' )
    line++;
  field = line;
  while ( *line != ' ' && *line != '\0' )
    line++;
  *line++ = '\0'; 
  /* get definition of the specified atom */
  assign_type_and_orbit ( field, &rule->type, &rule->orbit);
  for ( i = 0; i < MAX_RULES_PER_BONDED_ATOM; i++ )
    rule->prohibited[i] = rule->required[i] = NULL;
  /* grep all rules for this atom */
  while ( *line != '\0' )
    {
      if ( *line == '[' )
	{
	  if ( p_counter >= MAX_RULES_PER_BONDED_ATOM )
	    err_panic2 ( "parse_definition_line",
			"too many prohibitive rules for a single atom");
	  rule->prohibited[p_counter] = 
	    (rule_list_pt) mymalloc ( sizeof (rule_list_t) );
	  line = get_rule_item ( line, rule->prohibited[p_counter]);
	  p_counter++;
	}
      if ( *line == '(' )
	{
	  if ( r_counter >= MAX_RULES_PER_BONDED_ATOM )
	    err_panic2 ( "parse_definition_line",
			"too many required rules for a single atom");
	  rule->required[r_counter] = 
	    (rule_list_pt) mymalloc ( sizeof (rule_list_t) );
	  line = get_rule_item ( line, rule->required[r_counter]);
	  r_counter++;
	}
      line++;
    }  
}

void  print_rule_item ( rule_list_pt rule )
{
  printf ( "%d %s %s ",
	   rule->number,
	   rule->type == C ? "C"
	   : rule->type == O ? "O"
	   : rule->type == F ? "F"
	   : rule->type == CL ? "Cl"
	   : rule->type == N ? "N"
	   : rule->type == S ? "S"
	   : rule->type == P ? "P"
	   : rule->type == H ? "H"
	   : rule->type == SI ? "SI" : "ANY",
	   rule->orbit == SP1 ? "SP1"
	   : rule->orbit == SP2 ? "SP2"
	   : rule->orbit == SP3 ? "SP3"
	   : rule->orbit == SP4 ? "SP4"
	   : rule->orbit == AR ? "AR"
	   : rule->orbit == CO2 ? "CO2"
	   : rule->orbit == AM ? "AM"
	   : rule->orbit == PL3 ? "PL3"
	   : rule->orbit == CAT ? "CAT" : "ANY" );
  if ( rule->required != NULL )
    {
      printf ( "required: " );
      print_rule_item ( rule->required );
    }
  if ( rule->prohibited != NULL )
    {
      printf ( "prohibited: " );
      print_rule_item ( rule->prohibited );
    }
}  

void  print_rule ( rule_pt  rule )
{
  int  i;

  printf ( "   atom : %s %s\n", 
	   rule->type == C ? "C"
	   : rule->type == O ? "O"
	   : rule->type == N ? "N"
	   : rule->type == S ? "S"
	   : rule->type == F ? "F"
	   : rule->type == CL ? "Cl"
	   : rule->type == P ? "P"
	   : rule->type == H ? "H"
	   : rule->type == SI ? "SI" : "ANY",
	   rule->orbit == SP1 ? "SP1"
	   : rule->orbit == SP2 ? "SP2"
	   : rule->orbit == SP3 ? "SP3"
	   : rule->orbit == SP4 ? "SP4"
	   : rule->orbit == AR ? "AR"
	   : rule->orbit == CO2 ? "CO2"
	   : rule->orbit == AM ? "AM"
	   : rule->orbit == PL3 ? "PL3"
	   : rule->orbit == CAT ? "CAT" : "ANY" );
  for ( i = 0; i < MAX_RULES_PER_BONDED_ATOM; i++ )
    {
      if ( rule->required[i] != NULL )
	{
	  printf ( "   required: " );
	  print_rule_item ( rule->required[i] );
	  printf ( "\n" );
	}
      if ( rule->prohibited[i] != NULL )
	{
	  printf ( "   prohibited: " );
	  print_rule_item ( rule->prohibited[i] );
	  printf ( "\n" );
	}
    }
}

