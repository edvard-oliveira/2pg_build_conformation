#include <stdio.h>

#include"nerfio.h"
#include"nerf.h"
#include"futil.h"
#include"enums.h"
#include <stdlib.h>

void save_nerf_file(const char *path, const char *file_name,
		const int *atom_index, const int *atm_num_topol_A, const char *atm_name_A,
		__const char *atm_name_B,
		__const char *atm_name_C, __const char *atm_name_D,
		const own_vector_t *D, const double *bond_len_BC,
		const double *bond_len_CD, const double *bond_angle_BCD,
		const double *torsion_BC, const own_vector_t *A,
		const own_vector_t *B, const own_vector_t *C){
	FILE *nerf_file=NULL;
	char *fname = path_join_file(path,file_name);
	nerf_file = open_file(fname, fAPPEND);
	if (file_is_empty(nerf_file) == btrue){
		write_header_nerf(nerf_file);
	}
	write_line_nerf(nerf_file, atom_index, atm_num_topol_A, atm_name_A,
			atm_name_B, atm_name_C,
			atm_name_D, D,	bond_len_BC, bond_len_CD,
			bond_angle_BCD, torsion_BC, A, B, C);
	free(fname);
	fclose(nerf_file);
}

static void write_header_nerf(FILE *nerf_file){
	fprintf(nerf_file,"Atom_Index \t Atom_Num_Topol_A \t Atom_Name_A \t Atom_Name_B \t Atom_Name_C \t Atom_Name_D \t bond_len_BC \t bond_len_CD \t bond_angle_BCD \t torsion_BC \t Ax \t Ay \t Az \t Bx \t By \t Bz \t Cx \t Cy \t Cz \t Dx \t Dy \t Dz \n");
}

static void write_line_nerf(FILE *nerf_file, const int *atom_index,
		const int *atm_num_topol_A,
		const char *atm_name_A, __const char *atm_name_B, __const char *atm_name_C,
		__const char *atm_name_D, 	const own_vector_t *D,
		const double *bond_len_BC,	const double *bond_len_CD,
		const double *bond_angle_BCD, 	const double *torsion_BC,
		const own_vector_t *A, 	const own_vector_t *B, const own_vector_t *C){
	fprintf(nerf_file,"%i \t %i \t %s \t %s \t %s \t %s \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \t %8.5f \n",
			*atom_index, *atm_num_topol_A, atm_name_A, atm_name_B,
			atm_name_C, atm_name_D,
			*bond_len_BC, *bond_len_CD, *bond_angle_BCD, *torsion_BC,
			A->x, A->y, A->z,
			B->x, B->y, B->z,
			C->x, C->y, C->z,
			D->x, D->y, D->z);
}
