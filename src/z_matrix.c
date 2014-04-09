#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "z_matrix_types.h"
#include "z_matrix_lib.h"
#include "z_matrix_io.h"
#include "topology_types.h"
#include "z_matrix.h"
#include "messages.h"

z_matrix_global_t * allocateZ_matrix(const int *num_atoms){
	z_matrix_global_t *aux;
	aux = Malloc(z_matrix_global_t, 1);
	aux->z_matrix_info = Malloc(z_matrix_information_t, *num_atoms);
	//aux->prot_back = prot_back;
	aux->num_elements = *num_atoms;//-1;
	initialize_z_matrix_info(aux, num_atoms);
	return aux;
}

void initialize_z_matrix_info(z_matrix_global_t * aux, const int *num_atoms){
	int i;
	for (i = 0; i < *num_atoms; i++){
		aux->z_matrix_info[i].atom_angle = -1;
		aux->z_matrix_info[i].atom_connected = -1;
		aux->z_matrix_info[i].atom_reference = -1;
		aux->z_matrix_info[i].atomid = atmNR;
		strcpy(aux->z_matrix_info[i].atomname,"");
		aux->z_matrix_info[i].bond_angle = 0.00;
		aux->z_matrix_info[i].bond_len = 0.00;
		aux->z_matrix_info[i].bond_len_2 = 0.00;
		aux->z_matrix_info[i].dihedral_angle = 0.00;
		aux->z_matrix_info[i].dihedral_connect = -1;
		aux->z_matrix_info[i].index_top = -1;
		aux->z_matrix_info[i].tpAngle = angl_typ_dieh_0;
	}
}

void desAllocateZ_matrix(z_matrix_global_t * z_matrix){
        free(z_matrix->z_matrix_info);
	free(z_matrix);
}

void build_z_matrix(z_matrix_global_t *z_matrix, const top_global_t *top){
	_build_z_matrix(z_matrix, top);
}

void build_z_matrix_DATABASE(z_matrix_global_t *z_matrix,
		const top_global_t *top, const char *path){
	_build_z_matrix_database(z_matrix, top, path);
}

void show_z_matrix(const z_matrix_global_t *z_matrix){
	_show_z_matrix(z_matrix);
}

void save_z_matrix_file(const char *path, const char *file_name,
			const z_matrix_global_t *z_matrix){
	_save_z_matrix_file(path, file_name, z_matrix);
}

static void set_values_z_matriz_information(z_matrix_information_t *dest, 
	const z_matrix_information_t *source){
    dest->atom_angle       = source->atom_angle;
	dest->atom_connected   = source->atom_connected;
	dest->atom_reference   = source->atom_reference;
	dest->atomid           = source->atomid;		
	dest->bond_angle       = source->bond_angle;
	dest->bond_len         = source->bond_len;
	dest->bond_len_2       = source->bond_len_2;
	dest->dihedral_angle   = source->dihedral_angle;
	dest->dihedral_connect = source->dihedral_connect;
	dest->index_top        = source->index_top;
	dest->tpAngle          = source->tpAngle;
	strcpy(dest->atomname, source->atomname);
}
void copy_z_matrix(z_matrix_global_t *dest, const z_matrix_global_t *source){
	/* This function copies from source to dest. Both are z_matrix_global_t
	*/
	if (source == NULL){
		if (dest != NULL){
			fatal_error("In copy_z_matrix source is null, but dest is not null. \n");
		}
	}else if (dest == NULL){
		fatal_error("In copy_z_matrix dest is null. Please, allocate it. \n");
	}else{ //source and dest are not null
    	int i;
	    for (i = 0; i < source->num_elements; i++){
	    	set_values_z_matriz_information(&dest->z_matrix_info[i], &source->z_matrix_info[i]);
		}
	}
}

void copy_z_matrix_residue(z_matrix_global_t *dest, const z_matrix_global_t *source, 
	const int *first_atom_residue, const int *last_atom_residue){
	/* This function copies from source to dest. Both are z_matrix_global_t. 
	 * first_atom_residue means the first atom of residue. 
	 * last_atom_residue means the last atom of residue. 
	 * These values are obtained from topology.
	*/
	if (source == NULL){
		if (dest != NULL){
			fatal_error("In copy_z_matrix source is null, but dest is not null. \n");
		}
	}else if (dest == NULL){
		fatal_error("In copy_z_matrix dest is null. Please, allocate it. \n");
	}else{ //source and dest are not null
    	int i;
	    for (i = *first_atom_residue-1; i < *last_atom_residue; i++){
	    	set_values_z_matriz_information(&dest->z_matrix_info[i], &source->z_matrix_info[i]);
		}
	}
}
