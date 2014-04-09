#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "load_parameters.h"
#include "LoadConfig.h"
#include "messages.h"
#include "functions.h"
#include "string_owner.h"

static void initialize_parameters(input_parameters_t *param){
	param->seq_protein_file_name = Malloc(char, MAX_PATH_FILE_NAME );
	param->file_final_pdb = Malloc(char, MAX_FILE_NAME );
	param->path_database =  Malloc(char, MAX_PATH );
	param->path_local_execute = Malloc(char, MAX_PATH );	
	param->initial_pop_file_name = Malloc(char, MAX_FILE_NAME );
	param->z_matrix_file = Malloc(char, MAX_FILE_NAME );
	param->path_gromacs_programs = Malloc(char, MAX_PATH );
	param->gromacs_energy_min = ener_min_none;
	param->rotamer_library = rotamer_library_none;
}


static void set_parameter_gromacs_minimization(input_parameters_t *param,
		char *param_value){
	if (is_equal(param_value,"none") ){
		param->gromacs_energy_min = ener_min_none;
	}else if (is_equal(param_value,"ener_implicit") ){
		param->gromacs_energy_min = ener_min_implicit;
	}else if (is_equal(param_value,"ener_explicit") ){
		param->gromacs_energy_min = ener_min_explicit;
	}else{
		char msg[300];
		sprintf(msg, "gromacs_energy_min parameter was typed %s. But, this value must be none, ener_implicit or ener_explicit \n",param_value);
		fatal_error(msg);
	}
}


static void set_parameter_rotamer_library(input_parameters_t *param,
		char *param_value){
	if (is_equal(param_value,"none") ){
		param->rotamer_library = rotamer_library_none;
	}else if (is_equal(param_value,"cad_tuffery") ){
		param->rotamer_library = rotamer_library_cad_tuffery;
	}else{
		char msg[300];
		sprintf(msg, "rotamer_library parameter was typed %s. But, this value must be either none or cad_tuffery \n",param_value);
		fatal_error(msg);
	}
}


void deAllocateload_parameters(input_parameters_t *param){
  
    free(param->seq_protein_file_name );
    free(param->file_final_pdb );
	free(param->path_database );
	free(param->path_local_execute);	
	free(param->initial_pop_file_name );
	free(param->z_matrix_file);	
	free(param->path_gromacs_programs);

}

void load_parameters_from_file(input_parameters_t *param,
		const char *conf_file_name){
	/*Loading the configuration from file*/

	initialize_parameters(param);

	LoadConfig conf(conf_file_name);

	param->size_population = atoi(conf.getParameter("SizePopulation").c_str());
	strcpy(param->seq_protein_file_name, conf.getParameterChar("SequenceAminoAcidsPathFileName"));
	strcpy(param->path_database, conf.getParameterChar("Database"));
	strcpy(param->path_local_execute, conf.getParameterChar("Local_Execute"));
	strcpy(param->initial_pop_file_name, conf.getParameterChar("IniPopFileName"));
	strcpy(param->z_matrix_file, conf.getParameterChar("z_matrix_fileName"));
	strcpy(param->path_gromacs_programs, conf.getParameterChar("Path_Gromacs_Programs"));
	set_parameter_gromacs_minimization(param,conf.getParameterChar("gromacs_energy_min"));		
	set_parameter_rotamer_library(param,conf.getParameterChar("rotamer_library"));
	
}