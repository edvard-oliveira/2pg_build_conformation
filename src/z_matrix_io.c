#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "defines.h"
#include "z_matrix_types.h"
#include "z_matrix_io.h"
#include "z_matrix_lib.h"
#include "messages.h"
#include "futil.h"
#include "math_owner.h"

static char *stype;

void _show_z_matrix(const z_matrix_global_t *z_matrix){
	char msg[100];
	build_header_z_matrix(msg);
	display_msg(msg);
	for (int i = 0; i < z_matrix->num_elements;i++){
		build_line_output_z_matrix(msg,&i,z_matrix);
    	display_msg(msg);
	}
}

static void write_z_matrix_file(FILE *z_matrix_file,
		const z_matrix_global_t *z_matrix){
	char msg[500];
	build_header_z_matrix(msg);
	fprintf(z_matrix_file, "%s", msg);
	for (int i = 0; i < z_matrix->num_elements; i++){
		build_line_output_z_matrix(msg,&i,z_matrix);
		fprintf(z_matrix_file,"%s",msg);
	}
}

static void build_line_output_z_matrix(char *line, const int *i,
		const z_matrix_global_t *z_matrix){
	int index_atom;
	index_atom = *i;
	sprintf(line,"");
	float aux_bond_angle, aux_dihedral_angle, aux_rad;
	if (index_atom == 0){
		sprintf(line,"%i %i %i %s\n",
				z_matrix->z_matrix_info[index_atom].index_top,
				index_atom+1,
				z_matrix->z_matrix_info[index_atom].atom_reference,
				z_matrix->z_matrix_info[index_atom].atomname);

	}else if (index_atom == 1){
		sprintf(line,"%i %i %i %s %i %f \n", //%f
				z_matrix->z_matrix_info[index_atom].index_top,
				index_atom+1,
				z_matrix->z_matrix_info[index_atom].atom_reference,
				z_matrix->z_matrix_info[index_atom].atomname,
				z_matrix->z_matrix_info[index_atom].atom_connected,
		        z_matrix->z_matrix_info[index_atom].bond_len
		        //,z_matrix->z_matrix_info[index_atom].bond_len_2
		        );
	}else if (index_atom == 2){
		aux_rad = radians2degree_double(&z_matrix->z_matrix_info[index_atom].bond_angle);
		sprintf(line,"%i %i %i %s %i %f %f %f\n",
				z_matrix->z_matrix_info[index_atom].index_top,
				index_atom+1,
				z_matrix->z_matrix_info[index_atom].atom_reference,
				z_matrix->z_matrix_info[index_atom].atomname,
				z_matrix->z_matrix_info[index_atom].atom_connected,
		        z_matrix->z_matrix_info[index_atom].bond_len,
		        aux_rad,
		        z_matrix->z_matrix_info[index_atom].bond_len_2);
	}else{
		stype = Malloc(char,11);
		_get_type_of_dihedral(stype,z_matrix,&index_atom);
		aux_rad = radians2degree_double(&z_matrix->z_matrix_info[index_atom].bond_angle);
		sprintf(line,"%i %i %i %s %i %f %i %lf %i %s %f\n", //%f
				z_matrix->z_matrix_info[index_atom].index_top,
				index_atom+1,
				z_matrix->z_matrix_info[index_atom].atom_reference,
				z_matrix->z_matrix_info[index_atom].atomname,
				z_matrix->z_matrix_info[index_atom].atom_connected,
		        z_matrix->z_matrix_info[index_atom].bond_len,
		        z_matrix->z_matrix_info[index_atom].atom_angle,
		        aux_rad,
				z_matrix->z_matrix_info[index_atom].dihedral_connect,
				//radians2degree(&z_matrix->z_matrix_info[index_atom].dihedral_angle),
				stype,
				z_matrix->z_matrix_info[index_atom].bond_len_2);
		free(stype);
	}
}

static void build_header_z_matrix(char *line){
	sprintf(line,"Index_Top Seq Atom_Reference Atom_Name Atom_Connected Bond_Len Atom_Angle Bond_Angle Atom_Dihedral Dihedral_Angle type_diedhral bond_len_2\n");
}

void _save_z_matrix_file(const char *path, const char *file_name,
		const z_matrix_global_t *z_matrix){
	FILE *z_matrix_file=NULL;
	char *fname = path_join_file(path,file_name);
	z_matrix_file = open_file(fname, fWRITE);
	write_z_matrix_file(z_matrix_file, z_matrix);
	free(fname);
	fclose(z_matrix_file);
}


void read_file_distance_parameters_z_matrix(z_matrix_distance_parameters_t *dist_para,
		const char *path, const int *res_num){
	FILE *file=NULL;
	int r;
	char *fname;
	char header[TAM_BLOCO];
	fname = path_join_file(path,"1FUG_distancias.dat");
	file = open_file(fname, fREAD);
	r = 0;
	while(file != NULL){
		    if (r == *res_num){
		    	break;
		    }
			fscanf(file,"%i %lf %lf %lf ",&dist_para[r].res_num, &dist_para[r].N_CA,
	    			&dist_para[r].CA_C, &dist_para[r].C_Nplus);
			r++;
	}
	free(fname);
	fclose(file);
}

void read_file_angle_parameters_z_matrix(z_matrix_angle_parameters_t *ange_para,
		const char *path, const int *res_num){
	FILE *file=NULL;
	int r;
	char *fname;
	char chBloco[TAM_BLOCO];
	fname = path_join_file(path,"1FUG_angulos.dat");
	file = open_file(fname, fREAD);
	r = 0;
	while(file != NULL){
		    if (r == *res_num){
		    	break;
		    }
			fscanf(file,"%i %lf %lf %lf ",&ange_para[r].res_num, &ange_para[r].N_CA_C,
    			&ange_para[r].CA_C_Nplus, &ange_para[r].C_Nplus_CAplus);
			r++;
	}
	free(fname);
	fclose(file);
}


