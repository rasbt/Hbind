#ifndef _ERR_HANDLE_H
#define _ERR_HANDLE_H

void set_err_filename(char *fname);
void  err_panic ( char  *function, char  *message );

void  err_error ( char  *function, char  *message );

void  err_warning ( char  *function, char  *message );

void  err_log ( char  *message );

void  err_usage ( char  *message );

void  err_unknown_atom ( char  *name, char  *residue, char  *residue_num );

void  err_unknown_residue ( char  *residue, char  *residue_num );

void  err_unknown_level ( char  *name, char  *residue, char  *residue_num );


void  err_print ( char  *message);

void  err_panic2 ( char  *function, char  *message);

void  err_error2 ( char  *function, char  *message);

void  err_warning2 ( char  *function, char  *message);

void  err_log2 ( char  *message );

void  err_usage2 ( char  *message );

void  err_unknown_atom2 ( char  *name, char  *residue, char  *residue_num);

void  err_unknown_residue2 ( char  *residue, char  *residue_num);

void  err_unknown_level2 ( char  *name, char  *residue, char  *residue_num);

#endif
 
