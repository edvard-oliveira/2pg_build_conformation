#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "build_initial_population.h"
#include "populationio.h"
#include "build_protein.h"
#include "build_protein_lib.h"
#include "build_protein_database_types.h"
#include "parameters_type.h"
#include "messages.h"
#include "protein.h"
#include "topology.h"
#include "defines.h"
#include "functions.h"
#include "string_owner.h"
#include "pdbatom.h"
#include "pdbio.h"
#include "gromacs.h"


static void build_name_pdb_file(char *file_name,
    const int *ind){
        char *ind_s;
    ind_s = Malloc(char, 10);    
    int2str(ind_s, ind);
    strcpy(file_name, PREFIX_PDB_FILE_NAME_AUX);
    strcat(file_name, ind_s);
    strcat(file_name, ".pdb");
    free(ind_s);        
}

static void save_pdb_files(pdb_atom_t** pdb_pop, const int *numatom, 
    const int *size, const char *path ){

    const pdb_atom_t* pdb_ind;
    const pdb_seqres_t *seqres = NULL;
    char *file_name;    

    file_name = Malloc(char, MAX_FILE_NAME);
    for (int ind = 0; ind < *size; ind++){
        build_name_pdb_file(file_name, &ind);
        pdb_ind = pdb_pop[ind];
        save_pdb_file(path, file_name, numatom, pdb_ind, seqres);
    }
    free(file_name);
}

static void save_pop_file(const char *path, const char *file_name,
		protein ** pop, int pop_size){
	_save_population_file(path, file_name, pop, &pop_size);
}

static void create_message(char *msg,const int *ind){
	char aux[5];
	sprintf(aux, "%d", *ind+1);
	strcpy(msg,"Individual - ");
	strcat(msg,aux);
	strcat(msg,"\n");
}

void build_initial_population(const input_parameters_t *in_para){
	int nresiduos;
	int numatom;
	int num_bond_angle;
	int nr_kind_res;
	int num_protein_side_chains;
	int bond_angles;
    int numatom_after_min;
    char *pdbfile;
    int size_to_min;

	amino_t *primary_sequence;
	library_dihedral_info_t* lib_dihe_info;
	protein ** population_p;
	char message[200];
    int num_dihedral_angles_type;

    top_global_t *top_global=NULL;
    z_matrix_global_t *z_matrix=NULL;


    display_msg("Loading primary sequence \n");
    primary_sequence = load_protein2(in_para->seq_protein_file_name, &nresiduos,
    		&numatom,&bond_angles,&num_protein_side_chains,
    		&num_dihedral_angles_type, in_para);

    display_msg("Loading Database \n");
    lib_dihe_info = allocate_library_dihedral_info(nresiduos);
    initialize_build_protein(lib_dihe_info,primary_sequence,&nresiduos,
    		in_para->path_database,&nr_kind_res);

    display_msg("Build Topology \n");
	top_global = allocateTop_Global(&numatom, &nresiduos,&bond_angles,
			&num_protein_side_chains,&num_dihedral_angles_type);
    //build_top_global(primary_sequence,top_global);
    build_top_global_without_dihedral(primary_sequence,top_global);
    save_topology(in_para->path_local_execute, in_para->top_file,top_global);

    /* Here we have must buid the z matrix outside the protein struct because,
       here there is not known the population yet.
    */
    display_msg("Build z_matrix \n");
    z_matrix = allocateZ_matrix(&top_global->numatom);
    build_z_matrix(z_matrix,top_global);
    save_z_matrix_file(in_para->path_local_execute, in_para->z_matrix_file,
            z_matrix);

    display_msg("Building Population \n");
    population_p = allocatePopulation(in_para->size_population,nresiduos, 
        1, top_global->numatom);    
    for (int p = 0; p < in_para->size_population;p++){
    	create_message(message,&p);
    	display_msg(message);
    	if (in_para->rotamer_library == rotamer_library_cad_tuffery ){
        	build_random_protein_rotamer_library(population_p,&p,primary_sequence,
        			&nr_kind_res,lib_dihe_info,&nresiduos);
    	}else if (in_para->rotamer_library == rotamer_library_none){
        	build_random_protein_gsl(population_p,&p,primary_sequence,
        			&nr_kind_res,lib_dihe_info,&nresiduos);
    	}else{
    		fatal_error("Build population falied. Parameter not found \n");
    	}
        copy_z_matrix(population_p[p]->z_matrix, z_matrix);        
    }    
    //Saving Diehdral population to Cartesian Population
    pdb_atom_t** pdb_pop = allocate_Population_pdb(&in_para->size_population, &top_global->numatom);
    pdb_atom_t** pdb_pop_after_min = NULL;
    pdb_atom_t** pdb_pop_input_min = NULL;
    population2pdb_atom(pdb_pop, population_p, &in_para->size_population, top_global);
    display_msg("Saving file of before Minimization or Adding Hydrogen Atoms \n");
    save_pdb_files(pdb_pop, &top_global->numatom,  &in_para->size_population, 
        in_para->path_local_execute);    
    display_msg("Minimization or Adding Hydrogen Atoms\n");        
    init_gromacs_execution();
    pdbfile = Malloc(char, MAX_FILE_NAME);
    //Allocating newm population of atoms. It will be saved.
    int model;
    size_to_min = 1;
    pdb_pop_after_min = Malloc(pdb_atom_t*,size_to_min);
    pdb_pop_input_min = Malloc(pdb_atom_t*,size_to_min);
    pdb_pop_input_min[0] = allocate_pdbatom(&top_global->numatom);
    model = 1;
    for (int p = 0; p < in_para->size_population;p++){
        create_message(message,&p);
        display_msg(message);
        /* Apply minimization process or adding Hydrogen atoms in 
         * population that was built by Dihedral representation.
        */
        model = model + 1;
        build_name_pdb_file(pdbfile,&p); 
        load_pdb_file(pdb_pop_input_min[0], NULL, in_para->path_local_execute, pdbfile, 
            &top_global->numatom);        
        minimization_gromacs(pdb_pop_input_min[0], pdbfile, &numatom_after_min, in_para, 
            &top_global->numatom );
        /* Loading new conformation. This conformation either miminizated or adding 
         * Hydrogen atom by pdb2gmx. The number of atoms of this
         * conformation should changed because the number of Hydrogen atoms is based on
         * protonation state of residues.
        */
        pdb_pop_after_min[0] = allocate_pdbatom(&numatom_after_min);
        load_pdb_file(pdb_pop_after_min[0], NULL, in_para->path_local_execute, pdbfile, 
            &numatom_after_min);
        save_adding_model_pdb_file(in_para->path_local_execute, 
            in_para->initial_pop_file_name, 
        &model, &numatom_after_min, pdb_pop_after_min[0], NULL );            
        desAllocate_pdbatom(pdb_pop_after_min[0]);        
    }
    finish_gromacs_execution();
    free(pdbfile);
    display_msg("Creating the Population file \n");
    //save_model_pdb_file(in_para->path_local_execute,in_para->initial_pop_file_name, 
    //    &in_para->size_population, &numatom_after_min, pdb_pop_after_min, 
    //    NULL );    
}
