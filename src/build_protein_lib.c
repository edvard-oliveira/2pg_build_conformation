#include <stdio.h>

#include "build_initial_population.h"
#include "build_protein_database_types.h"
#include "build_protein_lib.h"
#include "populationio.h"
#include "futil.h"
#include "enums.h"
#include "messages.h"
#include "defines.h"
#include "randomlib.h"
#include "math_owner.h"

static amino_database_parameters_t amino_database_parameters [] = {
		{aARG,18283,"argenina.txt",80,4,"arg_lat.txt","%f %f %f %f %f %f"},
		{aALA,27265,"alanina.txt",1,1,"ala_lat.txt","%f %f %f %f"},
		{aGLY,25067,"glicina.txt",0,0,NULL,NULL},
		{aASN,16157,"asparagina.txt",7,2,"asn_lat.txt","%f %f %f %f"},
		{aASP,21265,"acido_aspartico.txt",8,2,"asp_lat.txt","%f %f %f %f"},
		{aCYS,5996,"cisteina.txt",3,1,"cys_lat.txt","%f %f %f"},
		{aGLN,14107,"glutanina.txt",32,3,"gln_lat.txt","%f %f %f %f %f"},
		{aGLU,24325,"acido_glutamico.txt",28,3,"glu_lat.txt","%f %f %f %f %f"},
		{aILE,21344,"isoleucina.txt",9,2,"ile_lat.txt","%f %f %f %f"},
		{aLEU,32483,"leucina.txt",10,2,"leu_lat.txt","%f %f %f %f"},
		{aLYS,21851,"lisina.txt",78,4,"lys_lat.txt","%f %f %f %f %f %f"},
		{aSER,22025,"serina.txt",3,1,"ser_lat.txt","%f %f %f"},
		{aMET,6880,"metionina.txt",29,3,"met_lat.txt","%f %f %f %f %f"},
		{aPHE,15476,"fenilalanina.txt",8,2,"phe_lat.txt","%f %f %f %f"},
		{aPRO,16321,"prolina.txt",12,2,"pro_lat.txt","%f %f %f %f"},
		{aTHR,20181,"treonina.txt",3,1,"thr_lat.txt","%f %f %f"},
		{aTRP,5493,"triptofano.txt",8,2,"trp_lat.txt","%f %f %f %f"},
		{aTYR,13042,"tirosina.txt",7,2,"tyr_lat.txt","%f %f %f %f"},
		{aVAL,25616,"valina.txt",3,1,"val_lat.txt","%f %f %f"},
		{aHIS,8788,"histidina.txt",8,2,"his_lat.txt","%f %f %f %f"},
		{aX,1,"",0,0,"",""}


};

/*This variable is applied to check if database is already loaded or not.
 * This concept is same singleton. However, it shows an error message instead of
 * controls the allocation.
 */
static database_initialized_t database_started;

void set_database_started(boolean_t b){
	/*Set value for database_started*/
	database_started.initialized = 	b;
}

database_initialized_t * get_database_started(){
	return &database_started;
}

library_dihedral_info_t* allocate_library_dihedral_info(int num_res){
	/*Receives the number of protein residues*/
	library_dihedral_info_t* aux;
	aux = Malloc(library_dihedral_info_t,num_res);
	for (int i=0; i < num_res;i++){
		aux[i].num_side_chains = 0;
		aux[i].num_torsional = 0;
		aux[i].num_angles_angles = 0;
		aux[i].torsional = NULL;
		aux[i].side_chains = NULL;
	}
    return aux;
}

library_dihedral_info_tors_t * allocate_library_dihedral_info_tors(int num_angles){
	library_dihedral_info_tors_t * aux;
	aux = Malloc(library_dihedral_info_tors_t,num_angles);
	/*
	for (int i = 0;i < num_angles;i++){
		aux[i].id = 0;
		aux[i].phi = 0;
		aux[i].psi = 0;
	}
	*/
	return aux;
}

library_dihedral_info_side_chains_t* allocate_library_dihedral_info_side_chains(int num_angles){
	library_dihedral_info_side_chains_t *aux;
	int i = 0;
	aux = Malloc(library_dihedral_info_side_chains_t,num_angles);
	for (i = 0; i < num_angles;i++){
		aux[i].side_chains_angles = allocate_side_chains(1);
	}
	return aux;

}

static void desAllocate_library_dihedral_info_tors (library_dihedral_info_tors_t *aux){
        free(aux);
	/*for (int r = 0; r< *nr_res; r++){
		for (int a = 0; a < lib_dihe_info[r].num_torsional;a++){
			free(&lib_dihe_info[r].torsional[a]);
		}
		free(&lib_dihe_info[r]);
	}*/
}

static void desAllocate_library_dihedral_info_side_chains(library_dihedral_info_side_chains_t *aux, int num_angles){
	for(int i =0; i< num_angles; i++)
	  free(aux[i].side_chains_angles);
	free(aux);
}

void desAllocate_library_dihedral_info(library_dihedral_info_t * lib_dihe_info, const int *nr_res){
	for (int r = 0; r< *nr_res; r++){
	        desAllocate_library_dihedral_info_tors(lib_dihe_info[r].torsional);
		desAllocate_library_dihedral_info_side_chains(lib_dihe_info[r].side_chains,lib_dihe_info[r].num_side_chains);
	}
	free(lib_dihe_info);	

	set_database_started(bfalse);
	//desAllocate_library_dihedral_info_tors(lib_dihe_info, nr_res);
	//free(lib_dihe_info);
}

static void _check_max_number(const long int *max_param,const long int *max_file){
	if (*max_param != *max_file){
		fatal_error("Please, check the number in database file and amino_database_parameters struct. They are different.");

	}
}

static void check_amino_type(type_aminos_t *amino){
	/* This function checks type of amino.
	 * It is necessary to save memory. HSE and HSD have the same internal
	 * coordinaties than HIS
	 */
	if ( (*amino == aHSE) ||  (*amino == aHSD)){
		*amino = aHIS;
	}
}

static int get_index_amino_database_parameters(type_aminos_t amino){
	/*Receives an amino.
	 * Returns its index in amino_database_parameters
	 */
	int aux = -1;
	int TAM = asize(amino_database_parameters);
	check_amino_type(&amino);
	for (int i = 0;i< TAM;i++){
		if (amino_database_parameters[i].aminoid == amino ){
			aux = i;
			break;
		}
	}
	if ( aux == -1){
		fatal_error("index not found for amino_database_parameters. See get_index_amino_database_parameters functions!!");
	}
	return aux;
}

int _get_index_library_dihedral_info(const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res){
	/*Receives amino, library of dihedral angles and the number of kind
	 * of residues in primary sequence
	 * Returns the index of amino in library*/
	int aux = -1;
	for (int i = 0;i < *nr_kind_res;i++){
		if (lib_dihe_info[i].aminoid == *amino){
			aux = i;
			break;
		}
	}
	if (aux == -1){
		fatal_error("Index not found for amino in _get_index_library_dihedral_info function. Please, check it.");
	}
	return aux;
}

float _get_side_chain_angle_from_lib(const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino){
	_check_side_chain_angle_library(index_amino, lib_dihe_info);
	if (*index_lib == 0){
		return lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if (*index_lib == 1){
		return lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*index_lib == 2){
		return lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
	}else if (*index_lib == 3){
		return lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi4;
	}else{
		fatal_error("In _get_side_chain_angle_from_lib the value of index_lib is more than number of side chains");
	}
}

void _load_amino_database(library_dihedral_info_t *lib_dihe, type_aminos_t amino,
		const int *index_res, const char *database){
/*Receives a pointer of library_dihedral_info_t, an amino and the database path
 * Read amino files which contain all informations about dihedrals angles (torsional and
 * side chains) from database.
 */
    float phi, psi;
    float chi1, chi2, chi3, chi4, chi5;
    float freq, despad;
    int fscanfError;
    FILE *tor_file;
    FILE *side_chain_file;
    char *pathFileName;
    char *pathFileName_side_chain;
    long int max_value;
    int index_amino;

    if (amino != aX){
	    //Obtain the index of amino_database_parameters for amino
	    index_amino = get_index_amino_database_parameters(amino);

	    //Opening torsional file database    
	    pathFileName = path_join_file(database,amino_database_parameters[index_amino].file_torsion);
	    tor_file = open_file(pathFileName,fREAD);
	    //Get the number of torsional angle from database file
		fscanfError = fscanf(tor_file,"%ld",&max_value);
		_check_max_number(&amino_database_parameters[index_amino].max_torsional_angles, &max_value);
		//Loading database file torsional angles
		lib_dihe[*index_res].aminoid = amino;
		lib_dihe[*index_res].torsional = allocate_library_dihedral_info_tors(max_value);
		lib_dihe[*index_res].num_torsional = max_value;
		for (int l = 0;l <max_value;l++ ){
			fscanfError = fscanf(tor_file,"%f %f",&phi, &psi);
			lib_dihe[*index_res].torsional[l].id = l;
			lib_dihe[*index_res].torsional[l].phi = degree2radians(&phi);
			lib_dihe[*index_res].torsional[l].psi = degree2radians(&psi);
		}

		lib_dihe[*index_res].side_chains = NULL;
		lib_dihe[*index_res].num_side_chains = amino_database_parameters[index_amino].max_side_chains_angles;
		lib_dihe[*index_res].num_angles_angles = amino_database_parameters[index_amino].num_side_angles;

		if (amino_database_parameters[index_amino].max_side_chains_angles > 0){
			pathFileName_side_chain = path_join_file(database,amino_database_parameters[index_amino].file_side_chains);
		    side_chain_file = open_file(pathFileName_side_chain,fREAD);
		    //Get the number of side_chains angle from database file
			fscanfError = fscanf(side_chain_file,"%ld",&max_value);
			_check_max_number(&amino_database_parameters[index_amino].max_side_chains_angles, &max_value);
			//Loading database file side chains angles
			lib_dihe[*index_res].side_chains = allocate_library_dihedral_info_side_chains(max_value);
			lib_dihe[*index_res].num_side_chains = max_value;
			lib_dihe[*index_res].num_angles_angles = amino_database_parameters[index_amino].num_side_angles;
			for (int l = 0;l <max_value;l++ ){
				_get_line_values(side_chain_file,&amino,amino_database_parameters[index_amino].format_file,&chi1,
						&chi2,&chi3,&chi4,&chi5,&freq,&despad,&fscanfError);
				_side_chain_database_line2lib_dieh_info(lib_dihe,index_res, &amino,
						&chi1,&chi2,&chi3,&chi4,&chi5,&freq,&despad,&l);
			}
			free(pathFileName_side_chain);
			fclose(side_chain_file);
		}
		free(pathFileName);
		fclose(tor_file);
    }else{
    	//Here residue is aX. So it is either ACE or NME.
		lib_dihe[*index_res].aminoid = amino;
		lib_dihe[*index_res].torsional = allocate_library_dihedral_info_tors(1);
		lib_dihe[*index_res].num_torsional = 1;
		lib_dihe[*index_res].torsional[0].phi = 0;
		lib_dihe[*index_res].torsional[0].psi = 0;

    	lib_dihe[*index_res].num_side_chains = 0;
    	lib_dihe[*index_res].side_chains = NULL;
    }
}

static void _get_line_values(FILE *side_chain_file, const type_aminos_t *amino, const char *format,
		float *chi1, float *chi2, float *chi3, float *chi4, float *chi5,
		float *freq, float *despad, int *fscanfError){
	/*Receives the variables which are necessary to obtain information from database file line*/
	if (*amino == aARG){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,chi3,chi4,freq,despad);
	}else if (*amino == aALA){
		*fscanfError = fscanf(side_chain_file,format,chi1,freq,despad);
	}else if (*amino == aASN){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aASP){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aCYS){
		*fscanfError = fscanf(side_chain_file,format,chi1,freq,despad);
	}else if (*amino == aGLN){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,chi3,freq,despad);
	}else if (*amino == aGLU){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,chi3,freq,despad);
	}else if (*amino == aILE){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aLEU){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aLYS){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,chi3,chi4,freq,despad);
	}else if (*amino == aSER){
		*fscanfError = fscanf(side_chain_file,format,chi1,freq,despad);
	}else if (*amino == aMET){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,chi3,freq,despad);
	}else if (*amino == aPHE){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aPRO){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aTHR){
		*fscanfError = fscanf(side_chain_file,format,chi1,freq,despad);
	}else if (*amino == aTRP){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aTYR){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}else if (*amino == aVAL){
		*fscanfError = fscanf(side_chain_file,format,chi1,freq,despad);
	}else if ( (*amino == aHIS) || (*amino == aHSE) || (*amino == aHSD)){
		*fscanfError = fscanf(side_chain_file,format,chi1,chi2,freq,despad);
	}
	else{
		fatal_error("Amino is not found.Please see it in _get_line_values \n");
	}
}

static void _side_chain_database_line2lib_dieh_info(library_dihedral_info_t *lib_dihe, const int *index_res,
		const type_aminos_t* amino,const float *chi1,const float *chi2,const float *chi3,const float *chi4,
		const float *chi5, const float *freq,const float *despad, const int *index_angles){
	if (*amino == aARG){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi3 = degree2radians(chi3);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi4 = degree2radians(chi4);
	}else if (*amino == aALA){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
	}else if (*amino == aASN){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aASP){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aCYS){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
	}else if (*amino == aGLN){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi3 = degree2radians(chi3);
	}else if (*amino == aGLU){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi3 = degree2radians(chi3);
	}else if (*amino == aILE){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aLEU){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aLYS){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi3 = degree2radians(chi3);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi4 = degree2radians(chi4);
	}else if (*amino == aSER){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
	}else if (*amino == aMET){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi3 = degree2radians(chi3);
	}else if (*amino == aPHE){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aPRO){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aTHR){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
	}else if (*amino == aTRP){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aTYR){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}else if (*amino == aVAL){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
	}else if ( (*amino == aHIS) || (*amino == aHSE) || (*amino == aHSD) ){
		lib_dihe[*index_res].side_chains[*index_angles].id = *index_angles;
		lib_dihe[*index_res].side_chains[*index_angles].despad = *despad;
		lib_dihe[*index_res].side_chains[*index_angles].freq = *freq;
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi1 = degree2radians(chi1);
		lib_dihe[*index_res].side_chains[*index_angles].side_chains_angles->chi2 = degree2radians(chi2);
	}
	else{
		fatal_error("Amino is not found.Please see it in _side_chain_database_line2lib_dieh_info \n");
	}
}

void _build_random_amino_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res){
	/*Receives an amino_t, amino, the library that contains all possible
	 * dihedrals angles and the number of kind of residues
	 * amino_t receives random values of dihdrals angles based on database
	 */
	_build_random_dihedral_angles_rotamer_library(amino_aux,amino,lib_dihe_info,nr_kind_res);
}

void _build_random_amino_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res){
	/*Receives an amino_t, amino, the library that contains all possible
	 * dihedrals angles and the number of kind of residues
	 * amino_t receives random values of dihdrals angles based on database
	 */
	_build_random_dihedral_angles_gsl(amino_aux,amino,lib_dihe_info,nr_kind_res);
}

static void _build_random_dihedral_angles_rotamer_library(amino_t *amino_aux,
		const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res){
	int index_amino;
	int rand_number;

	//Get index of amino in lib_dihe_info
	index_amino = _get_index_library_dihedral_info(amino,lib_dihe_info,nr_kind_res);

	//Get a random number for obtaining the torsional angles
	rand_number = _get_int_random_number(&lib_dihe_info[index_amino].num_torsional);
	//Set phi and psi
	amino_aux->phi = lib_dihe_info[index_amino].torsional[rand_number].phi;
	amino_aux->psi = lib_dihe_info[index_amino].torsional[rand_number].psi;
	amino_aux->number_late = 0;
	/*When side_chains is not null means this amino has side chains	*/
	if (lib_dihe_info[index_amino].side_chains != NULL){
		_set_random_side_chains_rotamer_library(amino_aux,amino,lib_dihe_info,&rand_number,&index_amino);
	}
}

static void _build_random_dihedral_angles_gsl(amino_t *amino_aux,
		const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res){
	int index_amino;
	int rand_number;

	//Get index of amino in lib_dihe_info
	index_amino = _get_index_library_dihedral_info(amino,lib_dihe_info,nr_kind_res);

	//Get a random number for obtaining the torsional angles
	rand_number = _get_int_random_number(&lib_dihe_info[index_amino].num_torsional);

	//Set phi and psi
	amino_aux->phi = _get_random_angle_gsl();
	amino_aux->psi = _get_random_angle_gsl();
	amino_aux->number_late = 0;
	/*When side_chains is not null means this amino has side chains	*/
	if (lib_dihe_info[index_amino].side_chains != NULL){
		_set_random_side_chains_gsl(amino_aux,amino,lib_dihe_info,&rand_number,&index_amino);
	}

}

static double _get_random_angle_gsl(){
	/* This function returns a random angle value based on GSL*/
	return _get_double_random_number()*3.14;
}

static void _check_side_chain_angle_library(const int *index_amino,
		const library_dihedral_info_t* lib_dihe_info){
	/* Checks side chains in index_amino. If not, a fatal error is called
	*/
	if (lib_dihe_info[*index_amino].side_chains == NULL){
		fatal_error("The side chains of amino is not allocated or this amino has not side chains.");
	}
}

static void _set_random_side_chains_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino){
	//If below guarantees this function is executed for amino which has side chains.
	_check_side_chain_angle_library(index_amino,lib_dihe_info);
	*index_lib = _get_int_random_number(&lib_dihe_info[*index_amino].num_side_chains);
	amino_aux->number_late = lib_dihe_info[*index_amino].num_angles_angles;
	if (amino_aux->late == NULL){
		amino_aux->late = Malloc(float,amino_aux->number_late);
	}else{
	        free(amino_aux->late);
		amino_aux->late = Malloc(float,amino_aux->number_late);
		//fatal_error("The side chains of amino_aux are not null which means that they have already allocated.");
	}
	if (*amino == aARG){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
		amino_aux->late[2] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
		amino_aux->late[3] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi4;
	}else if (*amino == aALA){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if (*amino == aASN){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aASP){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aCYS){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if (*amino == aGLN){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
		amino_aux->late[2] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
	}else if (*amino == aGLU){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
		amino_aux->late[2] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
	}else if (*amino == aILE){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aLEU){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aLYS){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
		amino_aux->late[2] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
		amino_aux->late[3] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi4;
	}else if (*amino == aSER){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if (*amino == aMET){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
		amino_aux->late[2] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi3;
	}else if (*amino == aPHE){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aPRO){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aTHR){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if (*amino == aTRP){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aTYR){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}else if (*amino == aVAL){
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
	}else if ( (*amino == aHIS) || (*amino == aHSE) || (*amino == aHSD) ) {
		amino_aux->late[0] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi1;
		amino_aux->late[1] = lib_dihe_info[*index_amino].side_chains[*index_lib].side_chains_angles->chi2;
	}
	else {
		fatal_error("In _set_random_side_chains_rotamer_library have to implemented side chains for new residue.");
	}

}


static void _set_random_side_chains_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino){
	//If below guarantees this function is executed for amino which has side chains.
	_check_side_chain_angle_library(index_amino,lib_dihe_info);
	*index_lib = _get_int_random_number(&lib_dihe_info[*index_amino].num_side_chains);
	amino_aux->number_late = lib_dihe_info[*index_amino].num_angles_angles;
	if (amino_aux->late == NULL){
		amino_aux->late = Malloc(float,amino_aux->number_late);
	}else{
	    free(amino_aux->late);
		amino_aux->late = Malloc(float,amino_aux->number_late);
		//fatal_error("The side chains of amino_aux are not null which means that they have already allocated.");
	}
	if (*amino == aARG){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
		amino_aux->late[2] =  _get_random_angle_gsl();
		amino_aux->late[3] =  _get_random_angle_gsl();
	}else if (*amino == aALA){
		amino_aux->late[0] =  _get_random_angle_gsl();
	}else if (*amino == aASN){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aASP){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aCYS){
		amino_aux->late[0] =  _get_random_angle_gsl();
	}else if (*amino == aGLN){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
		amino_aux->late[2] =  _get_random_angle_gsl();
	}else if (*amino == aGLU){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
		amino_aux->late[2] =  _get_random_angle_gsl();
	}else if (*amino == aILE){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aLEU){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aLYS){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
		amino_aux->late[2] =  _get_random_angle_gsl();
		amino_aux->late[3] =  _get_random_angle_gsl();
	}else if (*amino == aSER){
		amino_aux->late[0] =  _get_random_angle_gsl();
	}else if (*amino == aMET){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
		amino_aux->late[2] =  _get_random_angle_gsl();
	}else if (*amino == aPHE){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aPRO){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aTHR){
		amino_aux->late[0] =  _get_random_angle_gsl();
	}else if (*amino == aTRP){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aTYR){
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}else if (*amino == aVAL){
		amino_aux->late[0] =  _get_random_angle_gsl();
	}else if ( (*amino == aHIS) || (*amino == aHSE) || (*amino == aHSD) ) {
		amino_aux->late[0] =  _get_random_angle_gsl();
		amino_aux->late[1] =  _get_random_angle_gsl();
	}
	else {
		fatal_error("In _set_random_side_chains_gsl have to implemented side chains for new residue.");
	}

}

amino_t * _get_unique_res(int *nr_kind_res,const amino_t *primary_sequence,
		const int *nrresidues){
	/*Receives the primary sequence and its residue number.
	 * Returns the number of unique residue frm primary sequence and that sequece.
	 * Example: input ARGERA, output AREG
	 */
	amino_t * aux_seq;
	amino_t *amino_default = allocateAmino(21);

    //Obtain the number of kind residues. This number is used for
    //allocating aux_seq
    *nr_kind_res = 0;
    for (int r = 0; r < *nrresidues; r++){
    	if (_search_amino(amino_default, &primary_sequence[r].id,nr_kind_res) == bfalse){
    		copy_amino_allocating_only(&amino_default[*nr_kind_res],&primary_sequence[r]);
    		*nr_kind_res = *nr_kind_res +1;
    	}
    }
    aux_seq = allocateAmino(*nr_kind_res);
    for (int r = 0; r < *nr_kind_res; r++){
    	aux_seq[r].late = NULL; // It is forced to be NULL, because there is a check at copy_amino_allocating
    	copy_amino_allocating_only(&aux_seq[r],&amino_default[r]);
    }
    deAllocateAmino(amino_default,21);
	return aux_seq;
}

boolean_t _search_amino(const amino_t *aminos, const type_aminos_t* amino,
		const int *len_aminos){
	/* Receives a array of amino_t and an amino.
	 * Returns True when is found when occurence that amino in array.
	 * Otherwise, returns False.
	 */
	boolean_t aux = bfalse;
	for (int i = 0;i<*len_aminos;i++){
		if (aminos[i].id == *amino){
			aux = btrue;
			break;
		}
	}
	return aux;
}

