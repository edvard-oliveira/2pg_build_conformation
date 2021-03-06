#ifndef OLD_PARAMETERS_TYPE_H
#define OLD_PARAMETERS_TYPE_H

#include "enums.h"

typedef struct sinput_parameters {
    int size_population;
    char *seq_protein_file_name;
    char *z_matrix_file;
    char *top_file;
    char *initial_pop_file_name;
    type_energy_minimization_t gromacs_energy_min; //Set which energy minimization: none, implicit or explicit. none is default.
    char *file_final_pdb;
    char *path_database;
    char *path_local_execute;
    char *path_gromacs_programs;
    type_rotamer_library_t rotamer_library; // shows what kind of rotamer library
    int processor_number;
    char *mdp_file_min; // Indicates the name of mdp file
    char *force_field; // Indicates what is the force field
    /*set who can be neutralize N-Terminal*/
    type_terminal_charge_t n_terminal_charge;    
    /*set who can be neutralize C-Terminal*/    
    type_terminal_charge_t c_terminal_charge;    

 }input_parameters_t;

#endif
