#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "z_matrix.h"
#include "z_matrix_lib.h"
#include "z_matrix_types.h"
#include "topology_types.h"
#include "topology.h"
#include "messages.h"
#include "z_matrix_io.h"

static z_matrix_distance_parameters_t *dist_para = NULL;
static z_matrix_angle_parameters_t *angle_para = NULL;

static void _set_type_dihedral_angle(z_matrix_global_t *z_matrix,
		const int *index_z_matrix, const type_dihedral_angles_t *type_dihedral){
/* This function informs the type of dihedral angles for z matrix.
 * Basically, this information is already stored in topology structure. However,
 * you can build a z matrix information to maintain a better control for the
 * algorithm which converts dihedral space to cartezian space such as nerf
 * algorithm.
 * *type_dihedral is obtained at get_four_atom_number_for_bond_angle function
 */
	z_matrix->z_matrix_info[*index_z_matrix].tpAngle = *type_dihedral;
}

void _build_z_matrix(z_matrix_global_t *z_matrix, const top_global_t *top){
	type_aminos_t amino_id;
	int index_z_matrix = -1;
	z_matrix->num_elements = 0;

	for (int r = 1; r <= top->numres;r++){
		amino_id = get_amino_id_from_res_id(&r,top);
		if (r == 1){
			_add_atoms_backbone_N_Terminal(z_matrix,&index_z_matrix,&r,top);
		}else if ( (r > 1) && (r < top->numres) ){
			_add_atoms_backbone(z_matrix,&index_z_matrix,&r,top);
		}else if (r == top->numres){
			_add_atoms_backbone_C_terminal(z_matrix,&index_z_matrix,&r,top);
		}
		_add_hydrogen_atoms_backbone(z_matrix,&index_z_matrix,&r,&amino_id,
				top);
		_add_atoms_side_chains(z_matrix,&index_z_matrix,&r,&amino_id,top);
		_add_hydrogen_atoms_side_chains(z_matrix,&index_z_matrix,&r,&amino_id,
				top);

	}
	z_matrix->num_elements = index_z_matrix + 1;
}

static z_matrix_distance_parameters_t* allocate_z_matrix_distance_parameters(
		const int *numres){
	int i;
	z_matrix_distance_parameters_t *aux;
	aux = Malloc(z_matrix_distance_parameters_t, *numres);
	//Initialize values
	for (i = 0; i < *numres; i++){
		aux->res_num = -1;
		aux->CA_C = 0;
		aux->C_Nplus = 0;
		aux->N_CA = 0;
	}
	return aux;
}

static z_matrix_angle_parameters_t* allocate_z_matrix_angle_parameters(
		const int *numres){
	int i;
	z_matrix_angle_parameters_t *aux;
	aux = Malloc(z_matrix_angle_parameters_t, *numres);
	//Initialize values
	for (i = 0; i < *numres; i++){
		aux->res_num = -1;
		aux->CA_C_Nplus = 0;
		aux->C_Nplus_CAplus = 0;
		aux->N_CA_C = 0;
	}
	return aux;
}

void _build_z_matrix_database(z_matrix_global_t *z_matrix,
		const top_global_t *top, const char *path){
	/* Builds a Z matrix that has its parameters obtained from files.
	 * These parameters are distances and angles.
	 * _build_z_matrix function is a general parameters for Z matrix
	*/

	type_aminos_t amino_id;
	int index_z_matrix = -1;
	z_matrix->num_elements = 0;

	dist_para = allocate_z_matrix_distance_parameters(&top->numres);
	angle_para = allocate_z_matrix_angle_parameters(&top->numres);

	read_file_distance_parameters_z_matrix(dist_para, path, &top->numres);
	read_file_angle_parameters_z_matrix(angle_para, path, &top->numres);

	for (int r = 1; r <= top->numres;r++){
		amino_id = get_amino_id_from_res_id(&r,top);
		if (r == 1){
			_add_atoms_backbone_N_Terminal_DATABASE(z_matrix,&index_z_matrix,&r,top);
		}else if ( (r > 1) && (r < top->numres) ){
			_add_atoms_backbone_DATABASE(z_matrix,&index_z_matrix,&r,top);
		}else if (r == top->numres){
			_add_atoms_backbone_C_terminal_DATABASE(z_matrix,&index_z_matrix,&r,top);
		}
		_add_hydrogen_atoms_backbone(z_matrix,&index_z_matrix,&r,&amino_id,
				top);
		_add_atoms_side_chains(z_matrix,&index_z_matrix,&r,&amino_id,top);
		_add_hydrogen_atoms_side_chains(z_matrix,&index_z_matrix,&r,&amino_id,
				top);
	}
	z_matrix->num_elements = index_z_matrix + 1;

	free(dist_para);
}


static int get_index_atom_reference(const int *i, const top_global_t *top){
	/* the index of atom reference is the atom number of i less 1, always.
	 * atom_number[i] - 1
	 */
	int index = top->top_global_atom[*i].atom_number-1;
	return index;
}

void _get_type_of_dihedral(char stype[], const z_matrix_global_t *z_matrix,
		const int *i){
	/* This function returns string of type_of_diedhral
	 * Basically, it converts from enum to string representation
	 */
	if (z_matrix->z_matrix_info[*i].tpAngle == angl_phi){
		strcpy(stype,"PHI");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_psi){
		strcpy(stype,"PSI");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_psi_){
		strcpy(stype,"PSI_");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_omega){
		strcpy(stype,"OMEGA");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_chi1){
		strcpy(stype,"CHI1");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_chi2){
		strcpy(stype,"CHI2");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_chi3){
		strcpy(stype,"CHI3");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_chi4){
		strcpy(stype,"CHI4");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_1){
		strcpy(stype,"dieh_1");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_2){
		strcpy(stype,"dieh_2");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_3){
		strcpy(stype,"dieh_3");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_180){
		strcpy(stype,"dieh_180");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_90){
		strcpy(stype,"dieh_90");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_117_){
		strcpy(stype,"dieh_117_");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_0){
		strcpy(stype,"dieh_0");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_trans_120){
		strcpy(stype,"trans_120");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_trans_240){
		strcpy(stype,"trans_240");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_trans_123_){
		strcpy(stype,"trans_123_");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_4){
		strcpy(stype,"dieh_4");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_5){
		strcpy(stype,"dieh_5");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_6){
		strcpy(stype,"dieh_6");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_7){
		strcpy(stype,"dieh_7");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_29_6){
		strcpy(stype,"dieh_29_6");
	}else if (z_matrix->z_matrix_info[*i].tpAngle == angl_typ_dieh_37_4_){
		strcpy(stype,"dieh_37_4_");
	}
	else{
		fatal_error("Type of dihedral not found to display it\n");
	}
}

static void _set_information_atom_z_matrix(z_matrix_global_t *z_matrix,
		int *index_atom_reference,	int *index_z_matrix,
		const type_atoms_t *atom_reference, const int *res_id,
		const top_global_t *top){
/* Set generic information about atom reference at Z matrix. These information
 * are used for all case which an atom is added at Z matrix.
 */
	/*Obtain index of atom reference. It is (atom number at topology) -1 */
	*index_atom_reference = (get_num_atom_from_topol(res_id, *atom_reference,
			top) - 1);
	z_matrix->z_matrix_info[*index_z_matrix].atomid = top->top_global_atom[*atom_reference].atom_id;
	z_matrix->z_matrix_info[*index_z_matrix].atom_reference = top->top_global_atom[*index_atom_reference].atom_number;
	strcpy(z_matrix->z_matrix_info[*index_z_matrix].atomname,
			top->top_global_atom[*index_atom_reference].atom_name);
	z_matrix->z_matrix_info[*index_z_matrix].index_top = *index_atom_reference;
}


static void _add_atom_z_matriz(z_matrix_global_t *z_matrix, int *index_z_matrix,
		const int *res_id,
		type_atoms_t atom_reference, type_atoms_t atom_conected,
		type_atoms_t atom_angle,type_atoms_t atom_dihedral,
		type_dihedral_angles_t type_dihedral,
		float bond_len, float bond_angle, float bond_len_2,
		const top_global_t *top){
	/*bond_len_2 represents bond length BC (atom_angle bond with atom_conected*/
	type_atoms_t atom_aux;
	int res_id_aux, index_atom_reference;

	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_atom_reference,
			index_z_matrix, &atom_reference,res_id, top );
	set_residue_atom_values(&atom_aux,&res_id_aux,atom_conected,res_id);
	z_matrix->z_matrix_info[*index_z_matrix].atom_connected = get_num_atom_from_topol(&res_id_aux,
			atom_aux, top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_len = bond_len;
	z_matrix->z_matrix_info[*index_z_matrix].bond_len_2 = bond_len_2;
	set_residue_atom_values(&atom_aux,&res_id_aux,atom_angle,res_id);
	z_matrix->z_matrix_info[*index_z_matrix].atom_angle = get_num_atom_from_topol(&res_id_aux,
			atom_aux, top);
	/* PI - bond_angle means that we are working with complementary angle */
	z_matrix->z_matrix_info[*index_z_matrix].bond_angle = PI - bond_angle;
	set_residue_atom_values(&atom_aux,&res_id_aux,atom_dihedral,res_id);
	z_matrix->z_matrix_info[*index_z_matrix].dihedral_connect = get_num_atom_from_topol(&res_id_aux,
			atom_aux, top);
	_set_type_dihedral_angle(z_matrix, index_z_matrix, &type_dihedral);
}

static void _add_atoms_backbone(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	type_aminos_t amino_id;
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmN, atmC_,atmCA_,
			atmN_, angl_psi_, 1.30, 2.06, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmCA, atmN,atmC_,atmCA_,
			angl_typ_dieh_180, 1.49, 2.12, 1.30, top);
	_add_atom_z_matriz(z_matrix,index_z_matrix, res_id,atmC, atmCA, atmN, atmC_,
			angl_phi, 1.53, 1.9, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmO, atmC, atmCA,atmN,
			angl_typ_dieh_1, 1.21, 2.1145, 1.49, top);
}

static void _add_atoms_backbone_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/* res_id is started with 1. However, index is started with 0
	 * Therefore, res_id-1 means the reference residue and res_id-2 means
	 * last residue
	 */
	type_aminos_t amino_id;
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmN, atmC_,atmCA_,
			atmN_, angl_psi_, dist_para[*res_id-2].C_Nplus,
			angle_para[*res_id-2].CA_C_Nplus,
			dist_para[*res_id-2].CA_C, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmCA, atmN,atmC_,atmCA_,
			angl_typ_dieh_omega, dist_para[*res_id-1].N_CA,
			angle_para[*res_id-2].C_Nplus_CAplus,
			dist_para[*res_id-2].C_Nplus, top);
	_add_atom_z_matriz(z_matrix,index_z_matrix, res_id,atmC, atmCA, atmN, atmC_,
			angl_phi, dist_para[*res_id-1].CA_C,
			angle_para[*res_id-1].N_CA_C,
			dist_para[*res_id-1].N_CA, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmO, atmC, atmCA,atmN,
			angl_typ_dieh_1, 1.21, 2.1145, 1.49, top);
}

static void _add_atoms_backbone_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	type_atoms_t atom_reference;
	int index_reference;
	//Adding first atom - Nitrogen
	atom_reference = atmN;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	//Adding second atom - Alpha Carbon
	atom_reference = atmCA;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	z_matrix->z_matrix_info[*index_z_matrix].atom_connected = get_num_atom_from_topol(res_id,atmN,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_len = 1.49;
	//Adding third atom - Carbon
	atom_reference = atmC;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	z_matrix->z_matrix_info[*index_z_matrix].atom_connected = get_num_atom_from_topol(res_id,
			atmCA,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_len = 1.53;
	z_matrix->z_matrix_info[*index_z_matrix].bond_len_2 = 1.49;
	z_matrix->z_matrix_info[*index_z_matrix].atom_angle = get_num_atom_from_topol(res_id,
			atmN,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_angle = PI - 1.9042;
	//Adding four atom - Oxygen
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmO, atmC, atmCA,atmN,
			angl_typ_dieh_1, 1.21, 2.1145, 1.49, top);
}

static void _add_atoms_backbone_N_Terminal_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	type_atoms_t atom_reference;
	int index_reference;
	//Adding first atom - Nitrogen
	atom_reference = atmN;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	//Adding second atom - Alpha Carbon
	atom_reference = atmCA;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	z_matrix->z_matrix_info[*index_z_matrix].atom_connected = get_num_atom_from_topol(res_id,atmN,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_len = dist_para[0].N_CA;
	//Adding third atom - Carbon
	atom_reference = atmC;
	*index_z_matrix = *index_z_matrix + 1;
	_set_information_atom_z_matrix(z_matrix,&index_reference,index_z_matrix,
			&atom_reference, res_id,top);
	z_matrix->z_matrix_info[*index_z_matrix].atom_connected = get_num_atom_from_topol(res_id,
			atmCA,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_len = dist_para[0].CA_C;
	z_matrix->z_matrix_info[*index_z_matrix].bond_len_2 = 1.49;
	z_matrix->z_matrix_info[*index_z_matrix].atom_angle = get_num_atom_from_topol(res_id,
			atmN,top);
	z_matrix->z_matrix_info[*index_z_matrix].bond_angle = angle_para[0].N_CA_C;
	//Adding four atom - Oxygen
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmO, atmC, atmCA,atmN,
			angl_typ_dieh_1, 1.21, 2.1145, 1.49, top);
}

static void _add_atoms_backbone_C_terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmN, atmC_,atmCA_,
			atmN_, angl_psi_, 1.53, 2.0488, 1.30, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmCA, atmN,atmC_,atmCA_,
			angl_typ_dieh_180, 1.30, 2.0933, 1.49, top);
	_add_atom_z_matriz(z_matrix,index_z_matrix, res_id,atmC, atmCA, atmN, atmC_,
			angl_phi, 1.49, 1.9042, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmOT1, atmC, atmCA,atmN,
			angl_typ_dieh_3, 1.53, 2.1145, 1.21, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmOT2, atmC, atmCA,atmN,
			angl_typ_dieh_3, 1.53, 2.1145, 1.21, top);
}

static void _add_atoms_backbone_C_terminal_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmN, atmC_,atmCA_,
			atmN_, angl_psi_, dist_para[*res_id-2].C_Nplus,
			angle_para[*res_id-2].CA_C_Nplus,
			dist_para[*res_id-1].CA_C, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmCA, atmN,atmC_,atmCA_,
			angl_typ_dieh_180,
			dist_para[*res_id-1].N_CA,
			angle_para[*res_id-2].C_Nplus_CAplus ,
			dist_para[*res_id-2].C_Nplus, top);
	_add_atom_z_matriz(z_matrix,index_z_matrix, res_id,atmC, atmCA, atmN, atmC_,
			angl_phi,
			dist_para[*res_id-1].CA_C,
			angle_para[*res_id-1].N_CA_C,
			dist_para[*res_id-2].N_CA, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmOT1, atmC, atmCA,atmN,
			angl_typ_dieh_3, 1.53, 2.1145, 1.21, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id,atmOT2, atmC, atmCA,atmN,
			angl_typ_dieh_3, 1.53, 2.1145, 1.21, top);
}

static void _add_atoms_side_chains(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
const top_global_t *top){
	if (*amino_id == aALA){
		_add_atoms_side_chains_ALA(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aARG){
		_add_atoms_side_chains_ARG(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aASN){
		_add_atoms_side_chains_ASN(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aASP){
		_add_atoms_side_chains_ASP(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aCYS){
		_add_atoms_side_chains_CYS(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aGLN){
		_add_atoms_side_chains_GLN(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aGLU){
		_add_atoms_side_chains_GLU(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aILE){
		_add_atoms_side_chains_ILE(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aLEU){
		_add_atoms_side_chains_LEU(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aLYS){
		_add_atoms_side_chains_LYS(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aSER){
		_add_atoms_side_chains_SER(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aMET){
		_add_atoms_side_chains_MET(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aPHE){
		_add_atoms_side_chains_PHE(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aPRO){
		_add_atoms_side_chains_PRO(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aTHR){
		_add_atoms_side_chains_THR(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aTRP){
		_add_atoms_side_chains_TRP(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aTYR){
		_add_atoms_side_chains_TYR(z_matrix, index_z_matrix, res_id, top);
	}else if (*amino_id == aVAL){
		_add_atoms_side_chains_VAL(z_matrix, index_z_matrix, res_id, top);
	}else if ( (*amino_id == aHIS) || (*amino_id == aHSD) || (*amino_id == aHSE) ){
		_add_atoms_side_chains_HIS(z_matrix, index_z_matrix, res_id, top);
	}
}

static void _add_atoms_side_chains_HIS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){

	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top); // 110.5 = 1.9286
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.549, 1.9862, 1.53, top); // 113.8 = 1.9862
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmND1, atmCG, atmCB,
			atmCA, angl_chi2, 1.4, 2.1974, 1.53, top); // 125.9 = 2.1974
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE1, atmND1, atmCG,
			atmCB, angl_typ_dieh_180, 1.321, 1.7453, 1.4, top); // 100.0 = 1.7453
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNE2, atmCE1, atmND1,
			atmCG, angl_typ_dieh_180, 1.4, 4.222, 1.321, top); // 118.1 = 2.061
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD2, atmNE2, atmCE1,
			atmND1, angl_typ_dieh_0, 1.4, 1.7907, 1.4, top); // 102.6 = 1.7907
}


static void _add_atoms_side_chains_VAL(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.54, 1.9460, 1.49, top); // 111.5 = 1.9460
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG1, atmCB, atmCA,
			atmN, angl_chi1, 1.521, 1.9286, 1.54, top); // 110.5 = 1.9286
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG2, atmCB, atmCA,
			atmN, angl_typ_dieh_6, 1.521, 1.9286, 1.53, top); // 110.5 = 1.9286
}

static void _add_atoms_side_chains_TYR(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top); // 110.5 = 1.9286
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.512, 1.9879, 1.53, top); // 113.9 = 1.9879
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD1, atmCG, atmCB,
			atmCA, angl_chi2, 1.389, 2.1084, 1.512, top); // 120.8 = 2.1084
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD2, atmCG, atmCB,
			atmCA, angl_typ_dieh_4, 1.389, 2.1084, 1.512, top); // 120.8 = 2.1084
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE1, atmCD1, atmCG,
			atmCB, angl_typ_dieh_180, 1.382, 2.1153, 1.389, top); // 121.2 = 2.1153
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE2, atmCD2, atmCG,
			atmCB, angl_typ_dieh_180, 1.382, 2.1153, 1.389, top); // 121.2 = 2.1153
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCZ, atmCE1, atmCD1,
			atmCG, angl_typ_dieh_0, 1.378, 2.0874, 1.382, top); // 119.6 = 2.0874
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOH, atmCZ, atmCE1,
			atmCE2, angl_typ_dieh_180, 1.376, 2.0926, 1.378, top); // 119.9 = 2.0926
}

static void _add_atoms_side_chains_TRP(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top); // 110.5 = 1.9286
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.498, 1.9827, 1.53, top); // 113.6 = 1.9827
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD2, atmCG, atmCB,
			atmCA, angl_chi2, 1.4, 2.239, 1.53, top); // 128.3 = 2.239
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE2, atmCD2, atmCG,
			atmCB, angl_typ_dieh_180, 1.4, 1.8815, 1.4, top); // 107.8 = 1.8815
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNE1, atmCE2, atmCD2,
			atmCG, angl_typ_dieh_0, 1.4, 1.9338, 1.4, top); // 110.8 = 1.9338
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD1, atmNE1, atmCE2,
			atmCD2, angl_typ_dieh_0, 1.4, 1.8134, 1.4, top); // 103.9 = 1.8134
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCZ2, atmCE2, atmCD2,
			atmCG, angl_typ_dieh_180, 1.4, 2.094, 1.4, top); // 120.0 = 2.094
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCH2, atmCZ2, atmCE2,
			atmCD2, angl_typ_dieh_0, 1.4, 2.094, 1.4, top); // 120.0 = 2.094
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCZ3, atmCH2, atmCZ2,
			atmCE2, angl_typ_dieh_0, 1.4, 2.094, 1.4, top); // 120.0 = 2.094
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE3, atmCZ3, atmCH2,
			atmCZ2, angl_typ_dieh_0, 1.4, 2.094, 1.4, top); // 120.0 = 2.094

}

static void _add_atoms_side_chains_THR(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9460, 1.49, top); // 111.5 = 1.9460
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOG1, atmCB, atmCA,
			atmN, angl_chi1, 1.433, 1.9129, 1.53, top); // 109.6 = 1.9129
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG2, atmCB, atmCA,
			atmN, angl_typ_dieh_6, 1.521, 1.9286, 1.53, top); // 110.5 = 1.9286
}


static void _add_atoms_side_chains_PRO(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_117_, 1.53, 1.7977, 1.49, top); // 103 = 1.7977
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.492, 1.8239, 1.53, top); // 104.5 = 1.8239
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG, atmCB,
			atmCA, angl_chi2, 1.503, 1.8518, 1.492, top); // 106.1 = 1.8518

}


static void _add_atoms_side_chains_PHE(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.502, 1.9862, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD1, atmCG, atmCB,
			atmCA, angl_chi2, 1.384, 2.1066, 1.502, top); // 120.7 = 2.1066
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD2, atmCG, atmCB,
			atmCA, angl_typ_dieh_4 , 1.384, 2.1066, 1.502, top); // 120.7 = 2.1066
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE1, atmCD1, atmCG,
			atmCB, angl_typ_dieh_180 , 1.382, 2.1066, 1.384, top); // 120.7 = 2.1066
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE2, atmCD2, atmCG,
			atmCB, angl_typ_dieh_180 , 1.382, 2.1066, 1.384, top); // 120.7 = 2.1066
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCZ, atmCE1, atmCD1,
			atmCB, angl_typ_dieh_0 , 1.382, 2.0944, 1.382, top); // 120.0 = 2.0944
}


static void _add_atoms_side_chains_MET(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.52, 1.9914, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmSD, atmCG, atmCB,
			atmCA, angl_chi2, 1.803, 1.9670, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE, atmSD, atmCG,
			atmCB, angl_chi3, 1.791, 1.9356, 1.803, top);
}

static void _add_atoms_side_chains_SER(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOG, atmCB, atmCA,
			atmN, angl_chi1, 1.417, 1.9391, 1.53, top);
}

static void _add_atoms_side_chains_LYS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.52, 1.9914, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG, atmCB,
			atmCA, angl_chi2, 1.52, 1.9426, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCE, atmCD, atmCG,
			atmCB, angl_chi3, 1.52, 1.9426, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNZ, atmCE, atmCD,
			atmCG, angl_chi4, 1.489, 1.9530, 1.52, top);
}

static void _add_atoms_side_chains_LEU(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.54, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.53, 2.0298, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD1, atmCG, atmCB,
			atmCA, angl_chi2, 1.521, 1.9321, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD2, atmCG, atmCB,
			atmCA, angl_typ_dieh_7, 1.521, 1.9862, 1.53, top);
}

static void _add_atoms_side_chains_ILE(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.54, 1.9460, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG1, atmCB, atmCA,
			atmN, angl_chi1, 1.53, 1.9268, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG2, atmCB, atmCA,
			atmN, angl_typ_dieh_6, 1.521, 1.9286, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG1, atmCB,
			atmCA, angl_chi2, 1.513, 1.9862, 1.53, top);

}

static void _add_atoms_side_chains_GLU(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9120, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.52, 1.9914, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG, atmCB,
			atmCA, angl_chi2, 1.516, 1.9652, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOE1, atmCD, atmCG,
			atmCB, angl_chi3, 1.249, 2.0665, 1.516, top); //118.4 = 2.0665
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOE2, atmCD, atmCG,
			atmCB, angl_typ_dieh_5, 1.249, 2.0665, 1.516, top); //118.4 = 2.0665
}

static void _add_atoms_side_chains_GLN(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9120, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.52, 1.9914, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG, atmCB,
			atmCA, angl_chi2, 1.516, 1.9652, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOE1, atmCD, atmCG,
			atmCB, angl_chi3, 1.231, 2.1084, 1.516, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNE2, atmCD, atmCG,
			atmCB, angl_typ_dieh_5, 1.328, 2.2061, 1.516, top);
}

static void _add_atoms_side_chains_CYS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmSG, atmCB, atmCA,
			atmN, angl_chi1, 1.822, 1.9967, 1.53, top);
}


static void _add_atoms_side_chains_ASP(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.516, 1.9652, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOD1, atmCG, atmCB,
			atmCA, angl_chi2, 1.249, 2.0665, 1.516, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOD2, atmCG, atmCB,
			atmCA, angl_typ_dieh_4, 1.249, 2.0665, 1.516, top);
}

static void _add_atoms_side_chains_ASN(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmN,
			atmC, angl_typ_dieh_trans_123_, 1.53, 1.9286, 1.49, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.516, 1.9652, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmOD1, atmCG, atmCB,
			atmCA, angl_chi2, 1.231, 2.1084, 1.516, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmND2, atmCG, atmCB,
			atmCA, angl_typ_dieh_4, 1.328, 2.0316, 1.516, top);
}

static void _add_atoms_side_chains_ALA(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmC,
			atmN, angl_typ_dieh_2, 1.53, 1.9120, 1.53, top);
}

static void _add_atoms_side_chains_ARG(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top){
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCB, atmCA, atmC,
			atmN, angl_typ_dieh_2, 1.53, 1.9120, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCG, atmCB, atmCA,
			atmN, angl_chi1, 1.53, 1.9120, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCD, atmCG, atmCB,
			atmCA, angl_chi2, 1.52, 1.9426, 1.53, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNE, atmCD, atmCG,
			atmCB, angl_chi3, 1.46, 1.9548, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmCZ, atmNE, atmCD,
			atmCG, angl_chi4, 1.33, 2.1677, 1.52, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNH1, atmCZ, atmNE,
			atmCD, angl_typ_dieh_0, 1.33, 2.0944, 1.30, top);
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmNH2, atmCZ, atmNE,
			atmCD, angl_typ_dieh_180, 1.33, 2.0944, 1.30, top);
}


static void _add_hydrogen_atoms_backbone(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top){
	/* Adds all Hydrogen atoms for residue backbone.
	 * The parameters were based on aminoacids.hdb file from charmm force field
	 * implemented at Gromacs 4.5.4. More details about these parameters can be
	 * found in gromacs manual section 5.6.4 Hydrogen Database. In Gromacs
	 * manual 4.5.4 the page is 118.
	 *
	 * Firstly is checked kind of amino. GLY has HA2 atom which is not present
	 * in other amino. Because of this, GLY has specific check.
	 */
	if (*amino_id == aGLY){
		if (*res_id == 1){//N-Terminal
			_add_hydrogen_atoms_backbone_GLY_N_Terminal(z_matrix,
					index_z_matrix, res_id, top);
		}else{
			_add_hydrogen_atoms_backbone_GLY(z_matrix, index_z_matrix, res_id,
					top);
		}
	}else{
		if (*res_id == 1){//N-Terminal
			if ( (*amino_id != aPRO) && (*amino_id != aHIS) )
				_add_hydrogen_atoms_backbone_N_Terminal(z_matrix, index_z_matrix,
						res_id, top);
			else{
				if (*amino_id == aPRO){
					_add_hydrogen_atoms_backbone_N_Terminal_PRO(z_matrix, index_z_matrix,
							res_id, top);
				}
				if (*amino_id == aHIS){
					_add_hydrogen_atoms_backbone_N_Terminal_HIS(z_matrix, index_z_matrix,
							res_id, top);
				}
			}
		}else  {
			if (*amino_id == aARG){
				_add_hydrogen_atoms_backbone_ARG(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aALA){
				_add_hydrogen_atoms_backbone_ALA(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aASN){
				_add_hydrogen_atoms_backbone_ASN(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aASP){
				_add_hydrogen_atoms_backbone_ASP(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aCYS){
				_add_hydrogen_atoms_backbone_CYS(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aGLN){
				_add_hydrogen_atoms_backbone_GLN(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aGLU){
				_add_hydrogen_atoms_backbone_GLU(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aILE){
				_add_hydrogen_atoms_backbone_ILE(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aLEU){
				_add_hydrogen_atoms_backbone_LEU(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aLYS){
				_add_hydrogen_atoms_backbone_LYS(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aSER){
				_add_hydrogen_atoms_backbone_SER(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aMET){
				_add_hydrogen_atoms_backbone_MET(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aPHE){
				_add_hydrogen_atoms_backbone_PHE(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aPRO){
				_add_hydrogen_atoms_backbone_PRO(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aTHR){
				_add_hydrogen_atoms_backbone_THR(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aTRP){
				_add_hydrogen_atoms_backbone_TRP(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aTYR){
				_add_hydrogen_atoms_backbone_TYR(z_matrix, index_z_matrix, res_id,top);
			}else if (*amino_id == aVAL){
				_add_hydrogen_atoms_backbone_VAL(z_matrix, index_z_matrix, res_id,top);
			}else if ( (*amino_id == aHIS) || (*amino_id == aHSD) || (*amino_id == aHSE)){
				_add_hydrogen_atoms_backbone_HIS(z_matrix, index_z_matrix, res_id,top);
			}
		}
	}

}

static void _add_hydrogen_atoms_backbone_N_Terminal_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PRO when it is a N-Terminal
	 * According to Gromacs topology  this residue misses H3 atom.
	 */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH1, atmN, atmCA,
			atmC, angl_typ_dieh_180, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH2, atmN, atmCA,
			atmC, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH3, atmN, atmCA,
			atmC, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_N_Terminal_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PRO when it is a N-Terminal
	 * According to Gromacs topology  this residue misses H3 atom.
	 */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH1, atmN, atmCA,
			atmC, angl_typ_dieh_180, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH2, atmN, atmCA,
			atmC, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	//_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH3, atmN, atmCA,
//			atmC, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for N-Terminal*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH1, atmN, atmCA,
			atmC, angl_typ_dieh_180, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH2, atmN, atmCA,
			atmC, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH3, atmN, atmCA,
			atmC, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for HIS, HSD and HSE residue in backbone
	 * These conformations are equal. HSD conformation uses HD1 atom instead of
	 * HE2 atom. But these atoms are their side-chains. Therefore, their backbone
	 * are equal.
	 */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_VAL(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for VAL residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}


static void _add_hydrogen_atoms_backbone_TYR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for TYR residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_TRP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for TRP residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_THR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for THR residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PRO residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_PHE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PHE residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_MET(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for MET residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_SER(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for SER residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_LYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for LYS residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_LEU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for LEU residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_ILE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ILE residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_GLU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLU residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_GLN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLN residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_CYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for CYS residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_ASP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ASP residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_ASN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ASN residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_GLY(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLY residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA1, atmCA, atmC,
			atmN, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA2, atmCA, atmC,
			atmN, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_GLY_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLY residue in backbone
	 * GLY has HA2 atom
	 */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH1, atmN, atmCA,
			atmC, angl_typ_dieh_180, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH2, atmN, atmCA,
			atmC, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH3, atmN, atmCA,
			atmC, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA1, atmCA, atmC,
			atmN, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA2, atmCA, atmC,
			atmN, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
}
static void _add_hydrogen_atoms_backbone_ARG(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ARG residue in backbone */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_backbone_ALA(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ALA residue in backbone*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHN, atmN, atmC_,
			atmCA_, angl_typ_dieh_0, 1.0, 2.093, 1.30, top); // 2.093 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHA, atmCA, atmN,
			atmC, angl_typ_dieh_0, 1.0, 1.9106, 1.49, top); // 1.9106 = 109.47

}

static void _add_hydrogen_atoms_side_chains(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top){
	/* Adds all Hydrogen atoms for residue side-chains.
	 * The parameters were based on aminoacids.hdb file from charmm force field
	 * implemented at Gromacs 4.5.4. More details about these parameters can be
	 * found in gromacs manual section 5.6.4 Hydrogen Database. In Gromacs
	 * manual 4.5.4 the page is 118.
	 */
	 if (*amino_id == aARG){
		 _add_hydrogen_atoms_side_chains_ARG(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aALA){
		 _add_hydrogen_atoms_side_chains_ALA(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aASN){
		 _add_hydrogen_atoms_side_chains_ASN(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aASP){
		 _add_hydrogen_atoms_side_chains_ASP(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aCYS){
		  if ( (*res_id == 1) || (*res_id == top->numres) )
		    _add_hydrogen_atoms_side_chains_CYS(z_matrix, index_z_matrix, res_id,top);
		  else{
			  _add_hydrogen_atoms_side_chains_CYS_Internal(z_matrix,
					  index_z_matrix, res_id,top);
		  }
	 }else if (*amino_id == aGLN){
		 _add_hydrogen_atoms_side_chains_GLN(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aGLU){
		 _add_hydrogen_atoms_side_chains_GLU(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aILE){
		 _add_hydrogen_atoms_side_chains_ILE(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aLEU){
		 _add_hydrogen_atoms_side_chains_LEU(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aLYS){
		 _add_hydrogen_atoms_side_chains_LYS(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aSER){
		 _add_hydrogen_atoms_side_chains_SER(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aMET){
		 _add_hydrogen_atoms_side_chains_MET(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aPHE){
		 _add_hydrogen_atoms_side_chains_PHE(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aPRO){
		 _add_hydrogen_atoms_side_chains_PRO(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aTHR){
		 _add_hydrogen_atoms_side_chains_THR(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aTRP){
		 _add_hydrogen_atoms_side_chains_TRP(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aTYR){
		 _add_hydrogen_atoms_side_chains_TYR(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aVAL){
		 _add_hydrogen_atoms_side_chains_VAL(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aHSE){
		_add_hydrogen_atoms_side_chains_HSE(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aHSD){
			_add_hydrogen_atoms_side_chains_HSD(z_matrix, index_z_matrix, res_id,top);
	 }else if (*amino_id == aHIS){
			_add_hydrogen_atoms_side_chains_HIS(z_matrix, index_z_matrix, res_id,top);
	 }
}

static void _add_hydrogen_atoms_side_chains_VAL(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for VAL residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB, atmCB, atmCG1,
			atmCG2, angl_typ_dieh_0, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG11, atmCG1, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG12, atmCG1, atmCB,
			atmCA, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG13, atmCG1, atmCB,
			atmCA, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG21, atmCG2, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG22, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG23, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47

}

static void _add_hydrogen_atoms_side_chains_HSE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for HSE residue in side chains  */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
//	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmND1, atmCG,
//			atmCE1, angl_typ_dieh_0, 1.0, 1.9120, 1.378, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD2, atmCG,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.354, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE1, atmND1,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE2, atmNE2, atmCE1,
				atmCD2, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50

}

static void _add_hydrogen_atoms_side_chains_HSD(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for HSD residue in side chains  */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD2, atmCG,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.354, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE1, atmND1,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmND1, atmCG,
			atmCE1, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50
}

static void _add_hydrogen_atoms_side_chains_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for HIS residue in side chains  */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.549, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD2, atmCG,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.354, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE1, atmND1,
			atmNE2, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmH, atmND1, atmCG,
			atmCE1, angl_typ_dieh_0, 1.0, 1.9120, 1.321, top); // // 1.9120 = 109.50
}

static void _add_hydrogen_atoms_side_chains_TYR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for TYR residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.512, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.512, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD1, atmCG,
			atmCE1, angl_typ_dieh_0, 1.0, 1.9120, 1.389, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD2, atmCG,
			atmCE2, angl_typ_dieh_0, 1.0, 1.9120, 1.389, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE1, atmCD1,
			atmCZ, angl_typ_dieh_0, 1.0, 1.9120, 1.382, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE2, atmCE2, atmCD2,
			atmCZ, angl_typ_dieh_0, 1.0, 1.9120, 1.382, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH, atmOH, atmCZ,
			atmCE1, angl_typ_dieh_180, 1.0, 1.9120, 1.376, top); // 1.9120 = 109.50
}

static void _add_hydrogen_atoms_side_chains_TRP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for TRP residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.498, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.498, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD1, atmNE1,
			atmCG, angl_typ_dieh_0, 1.0, 1.9120, 1.374, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmNE1, atmCD1,
			atmCE2, angl_typ_dieh_0, 1.0, 1.9120, 1.374, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE3, atmCE3, atmCD2,
			atmCZ3, angl_typ_dieh_0, 1.0, 1.9120, 1.398, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ3, atmCZ3, atmCE3,
			atmCH2, angl_typ_dieh_0, 1.0, 1.9120, 1.382, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH2, atmCH2, atmCZ3,
			atmCZ2, angl_typ_dieh_0, 1.0, 1.9120, 1.4, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ2, atmCZ2, atmCE2,
			atmCH2, angl_typ_dieh_0, 1.0, 1.9120, 1.394, top); // 1.9120 = 109.50
}

static void _add_hydrogen_atoms_side_chains_THR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for THR residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB, atmCB, atmCA,
			atmOG1, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmOG1, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9111, 1.433, top); // 1.9111 = 109.5
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG21, atmCG2, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG22, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG23, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47

}


static void _add_hydrogen_atoms_side_chains_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PRO residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.492, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.492, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD, atmN,
			atmCG, angl_typ_dieh_0, 1.0, 1.9106, 1.473, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD, atmN,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.473, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.503, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmCD,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.503, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_PHE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for PHE residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD1, atmCG,
			atmCE1, angl_typ_dieh_0, 1.0, 1.5708, 1.384, top); // 90 = 1.5708
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD2, atmCG,
			atmCE2, angl_typ_dieh_0, 1.0, 1.5708, 1.384, top); // 90 = 1.5708
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE1, atmCD1,
				atmCZ, angl_typ_dieh_0, 1.0, 1.5708, 1.384, top); // 90 = 1.5708
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE2, atmCE2, atmCD2,
				atmCZ, angl_typ_dieh_0, 1.0, 1.5708, 1.384, top); // 90 = 1.5708
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ, atmCZ, atmCE1,
				atmCE2, angl_typ_dieh_0, 1.0, 1.5708, 1.384, top); // 90 = 1.5708
}

static void _add_hydrogen_atoms_side_chains_MET(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for MET residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmSD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.803, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmSD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.803, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE, atmSD,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.791, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE2, atmCE, atmSD,
			atmCG, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.791, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE3, atmCE, atmSD,
			atmCG, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.791, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_SER(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for SER residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmOG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmOG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmOG, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9111, 1.417, top); // 1.9111 = 109.50
}

static void _add_hydrogen_atoms_side_chains_LYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for LYS residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCA,
			atmCG, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCA,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmCB,
			atmCD, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmCB,
			atmCD, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD, atmCE,
			atmCG, angl_typ_dieh_0, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD, atmCE,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE1, atmCE, atmNZ,
			atmCD, angl_typ_dieh_0, 1.0, 1.9106, 1.489, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE2, atmCE, atmNZ,
			atmCD, angl_typ_dieh_180, 1.0, 1.9106, 1.489, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ1, atmNZ, atmCE,
			atmCD, angl_typ_dieh_180, 1.0, 1.9106, 1.489, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ2, atmNZ, atmCE,
			atmCD, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.489, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHZ2, atmNZ, atmCE,
			atmCD, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.489, top); // 1.9106 = 109.47
}


static void _add_hydrogen_atoms_side_chains_LEU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for LEU residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCA,
			atmCG, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCA,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG, atmCG, atmCB,
			atmCD1, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD11, atmCD1, atmCG,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD12, atmCD1, atmCG,
			atmCB, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD13, atmCD1, atmCG,
			atmCB, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD21, atmCD2, atmCG,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD22, atmCD2, atmCG,
			atmCB, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD23, atmCD2, atmCG,
			atmCB, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_ILE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ILE residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB, atmCB, atmCA,
			atmCG1, angl_typ_dieh_0, 1.0, 1.9106, 1.54, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG11, atmCG1, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.513, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG12, atmCG1, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.513, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG21, atmCG2, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG22, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG23, atmCG2, atmCB,
			atmCA, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.521, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD, atmCG1,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.513, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD, atmCG1,
			atmCB, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.513, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD3, atmCD, atmCG1,
			atmCB, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.513, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_GLU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLU residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.516, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmCD,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.516, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_GLN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for GLN residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.52, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.516, top); //1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmCD,
			atmCB, angl_typ_dieh_180, 1.0, 1.9106, 1.516, top); //1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE21, atmNE2, atmCD,
			atmCG, angl_typ_dieh_180, 1.0, 2.0944, 1.328, top); // 2.0944 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE22, atmNE2, atmCD,
			atmCG, angl_typ_dieh_0, 1.0, 2.0944, 1.328, top); // 2.0944 = 120.0
}

static void _add_hydrogen_atoms_side_chains_CYS_Internal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for CYS residue in side chains
	 * The word internal means CYS is NOT N or C terminal. According to
	 * aminoacids.rtp file this residue does not have HG1 atom.
	 */
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmSG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.822, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmSG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.822, top); // 1.9106 = 109.47
	//_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmSG, atmCB,
//			atmCA, angl_typ_dieh_180, 1.0, 1.9111, 1.822, top); // 1.9111 = 109.5
}

static void _add_hydrogen_atoms_side_chains_CYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for CYS residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmSG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.822, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmSG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.822, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmSG, atmCB,
			atmCA, angl_typ_dieh_180, 1.0, 1.9111, 1.822, top); // 1.9111 = 109.5
}

static void _add_hydrogen_atoms_side_chains_ASP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ASN residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
}

static void _add_hydrogen_atoms_side_chains_ASN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ASN residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD21, atmND2, atmCG,
			atmCB, angl_typ_dieh_0, 1.0, 2.0944, 1.328, top); // 2.0944 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD22, atmND2, atmCG,
			atmCB, angl_typ_dieh_180, 1.0, 2.0944, 1.328, top); // 2.0944 = 120.0
}

static void _add_hydrogen_atoms_side_chains_ARG(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ARG residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCG,
			atmCA, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCG,
			atmCA, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG1, atmCG, atmCD,
			atmCB, angl_typ_dieh_0, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHG2, atmCG, atmCB,
			atmN, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD1, atmCD, atmNE,
			atmCG, angl_typ_dieh_0, 1.0, 1.9106, 1.30, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHD2, atmCD, atmNE,
			atmCG, angl_typ_dieh_180, 1.0, 1.9106, 1.30, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHE, atmNE, atmCD,
			atmCZ, angl_typ_dieh_0, 1.0, 1.9120, 1.30, top); // 1.9120 = 109.50
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH11, atmNH1, atmCZ,
			atmNE, angl_typ_dieh_0, 1.0, 2.0944, 1.30, top); // 2.0944 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH12, atmNH1, atmCZ,
			atmNE, angl_typ_dieh_180, 1.0, 2.0944, 1.30, top); // 2.0944 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH21, atmNH2, atmCZ,
			atmNE, angl_typ_dieh_0, 1.0, 2.0944, 1.30, top); // 2.0944 = 120.0
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHH22, atmNH2, atmCZ,
			atmNE, angl_typ_dieh_180, 1.0, 2.0944, 1.30, top); // 2.0944 = 120.0
}

static void _add_hydrogen_atoms_side_chains_ALA(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top){
	/*Adds all Hydrogen atoms for ALA residue in side chains*/
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB1, atmCB, atmCA,
			atmN, angl_typ_dieh_180, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB2, atmCB, atmCA,
			atmN, angl_typ_dieh_trans_120, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
	_add_atom_z_matriz(z_matrix, index_z_matrix, res_id, atmHB3, atmCB, atmCA,
			atmN, angl_typ_dieh_trans_240, 1.0, 1.9106, 1.53, top); // 1.9106 = 109.47
}
