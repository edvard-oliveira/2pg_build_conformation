#include <string.h>

#include"build_protein.h"
#include "build_protein_lib.h"
#include "messages.h"
#include "pdbatom.h"
#include "pdbio.h"
#include "topology.h"
#include "string_owner.h"
#include "protein.h"
#include "defines.h"
#include "functions.h"
#include "osutil.h"

static void load_database_diehdral_angles(library_dihedral_info_t* lib_dihe_info, const amino_t *primary_sequence,
		const int *nrresidues, const char *path_database,int *nr_kind_res){
/*This function loads all informations of dihedral angles from database
 * Those information are called default informations and they are based on
 * library for dihedral angles.
 */
	amino_t *sequence_unique;
	sequence_unique = _get_unique_res(nr_kind_res,primary_sequence,nrresidues);

	for (int i=0; i < *nr_kind_res; i++) {
		_load_amino_database(lib_dihe_info, sequence_unique[i].id, &i,
				path_database);
        // free this memory
    }
	deAllocateAmino(sequence_unique,*nr_kind_res);

}

void initialize_build_protein(library_dihedral_info_t* lib_dihe_info, const amino_t *primary_sequence,
		const int *nrresidues, const char *path_database,int *nr_kind_res){
	if (get_database_started()->initialized == btrue){
		fatal_error("The database is already initialized. Check your code because you are calling twice\n");
	}
	load_database_diehdral_angles(lib_dihe_info, primary_sequence, nrresidues,
			path_database,nr_kind_res);
	set_database_started(btrue);
}

void build_random_protein_rotamer_library(protein ** pop,const int *index_pop,
		const amino_t *primary_sequence,const int *nr_kind_res,
		const library_dihedral_info_t* lib_dihe_info,
		const int *num_res){
    /*Build a random protein using a rotamer library such as CAD and Tuffery*/
	amino_t *amino_aux;
	for (int r  = 0; r < *num_res;r++){
		amino_aux = allocateAmino(1);
		amino_aux->id = primary_sequence[r].id;
		strcpy(amino_aux->aminoacido, primary_sequence[r].aminoacido);
		strcpy(amino_aux->idl3, primary_sequence[r].idl3);
		build_random_amino_rotamer_library(amino_aux,&primary_sequence[r].id,lib_dihe_info,
				nr_kind_res);
		set_amino_protein(pop[*index_pop] ,amino_aux, &r);
		deAllocateAmino(amino_aux,1);
	}
}

void build_random_protein_gsl(protein ** pop,const int *index_pop,
		const amino_t *primary_sequence,const int *nr_kind_res,
		const library_dihedral_info_t* lib_dihe_info,
		const int *num_res){
    /*Build a random protein whithout rotamer library*/
	amino_t *amino_aux;
	for (int r  = 0; r < *num_res;r++){
		amino_aux = allocateAmino(1);
		amino_aux->id = primary_sequence[r].id;
		strcpy(amino_aux->aminoacido, primary_sequence[r].aminoacido);
		strcpy(amino_aux->idl3, primary_sequence[r].idl3);
		build_random_amino_gsl(amino_aux,&primary_sequence[r].id,lib_dihe_info,
				nr_kind_res);
		set_amino_protein(pop[*index_pop] ,amino_aux, &r);
		deAllocateAmino(amino_aux,1);
	}
}

void build_random_amino_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res ){
	_build_random_amino_rotamer_library(amino_aux,amino,lib_dihe_info,
			nr_kind_res);
}

void build_random_amino_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res ){
	_build_random_amino_gsl(amino_aux,amino,lib_dihe_info,
			nr_kind_res);
}

