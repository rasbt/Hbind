#ifndef _READ_RULE_H
#define _READ_RULE_H
#include "types.h"

char  *get_rule_item ( char         *line, rule_list_pt rule);
  
void  parse_definition_line ( char    *line, rule_pt rule);

void  print_rule_item ( rule_list_pt rule );

void  print_rule ( rule_pt  rule );

#endif
