#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "defines.h"
#include "enums.h"
#include "messages.h"
#include "futil.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"
#include "topology_charmm27_parameters.h"
#include "topology_atom_parameters.h"
#include "pdbatom.h"
#include "nerf.h"
#include "string_owner.h"
#include "nerfio.h"
#include "osutil.h"
#include "randomlib.h"

const topol_residues_t *get_residue_from_topol(type_aminos_t amino_id){
	/*Receives an amino_id and returns its topol_residues_t as constant pointer*/
	if (amino_id >= aNR){
		fatal_error("Residue index is wrong when try to execute get_residue_from_topol function!!!");
	}
	const topol_residues_t * aux=NULL;
	int TAM = asize(topol_residues_ff);
	for (int i=0; i<TAM;i++){
		if (topol_residues_ff[i].res_id == amino_id){
			aux = &topol_residues_ff[i];
		}
	}
	if (aux==NULL){
		fatal_error("amino_id was not find when try to execute get_residue_from_topol function!!!");
	}
	return aux;
}

const topol_residues_t *get_residue_from_topol_N_Terminal(type_aminos_t amino_id){
	/*Receives an amino_id and returns its topol_residues_t as constant pointer*/
	if (amino_id >= aNR){
		fatal_error("Residue index is wrong when try to execute get_residue_from_topol_N_Terminal function!!!");
	}
	const topol_residues_t * aux=NULL;
	int TAM = asize(topol_residues_ff_N_Terminal);
	for (int i=0; i<TAM;i++){
		if (topol_residues_ff_N_Terminal [i].res_id == amino_id){
			aux = &topol_residues_ff_N_Terminal[i];
		}
	}
	if (aux==NULL){
		fatal_error("amino_id was not find when try to execute get_residue_from_topol_N_Terminal function!!!");
	}
	return aux;
}

const topol_residues_t *get_residue_from_topol_C_Terminal(type_aminos_t amino_id){
	/*Receives an amino_id and returns its topol_residues_t as constant pointer*/
	if (amino_id >= aNR){
		fatal_error("Residue index is wrong when try to execute get_residue_from_topol_C_Terminal function!!!");
	}
	const topol_residues_t * aux=NULL;
	int TAM = asize(topol_residues_ff_C_Terminal);
	for (int i=0; i<TAM;i++){
		if (topol_residues_ff_C_Terminal [i].res_id == amino_id){
			aux = &topol_residues_ff_C_Terminal[i];
		}
	}
	if (aux==NULL){
		fatal_error("amino_id was not find when try to execute get_residue_from_topol_C_Terminal function!!!");
	}
	return aux;
}


const top_global_atom_t* _get_top_global_atom_t_from_topol(const int *num_res, const type_atoms_t atm, const top_global_t *top){
	/* Receives residue number, atom_id and topology
	 * Returns a constant pointer for top_global_atom_t struct which contains all information about atom_id.
	 * Normally, this function is called by internal function from topology.c file such as get_num_atom_from_topol
	 * */
	top_global_atom_t* ret=NULL;
	/* 1) index of residue is *num_res-1 because index starts with 0 and *num_res is started with 1
	 * 2) top->top_global_res_atm[*num_res-1].atom_first -1 is necessary because atom_first is started with 1 and i
	 *     must be 0 because it represents the index for top_global which is started with 0.
	 */
	for (int i = top->top_global_res_atm[*num_res-1].atom_first -1; i < top->top_global_res_atm[*num_res-1].atom_last;i++ ){
		if (top->top_global_atom[i].atom_id == atm){
			ret = &top->top_global_atom[i];
			break;
		}
	}
	if (ret==NULL){
		char msg[200];
		sprintf(msg,"Index of atoms is wrong when try to execute get_atom_from_topol function!!!");
		fatal_error(msg);
	}
	return ret;
}

const special_atom_parameters_t * _get_special_atom_parameters(type_atoms_t atm){
	int index = _get_index_special_atom_parameters(atm);
	return &special_atom_parameters[index];
}

static int _get_index_bond_angles_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2,type_atoms_t atom_id_3){
	/* Receives three atoms which represent a bond angle. These angles must be informed
	 * in topol_atoms_bond_angles_parameters.
	 * This function returns the index of topol_atoms_bond_angles_parameters only.
	 * Access the angle value you have to call topol_atoms_bond_angles_parameters[index].angle_value
	 * like  _get_bond_angle_from_atom_parameters function.
	 */
	int aux = -1;
	int TAM = asize(topol_atoms_bond_angles_parameters);
	for (int i=0;i<TAM;i++){
		if ( (topol_atoms_bond_angles_parameters[i].atom_id_1 == atom_id_1) &&
			(topol_atoms_bond_angles_parameters[i].atom_id_2 == atom_id_2) &&
			(topol_atoms_bond_angles_parameters[i].atom_id_3 == atom_id_3) ){
            aux = i;
            break;
		}else if ( (topol_atoms_bond_angles_parameters[i].atom_id_2 == atom_id_1) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_1 == atom_id_2) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_3 == atom_id_3) ){
	            aux = i;
	            break;
		}else if ( (topol_atoms_bond_angles_parameters[i].atom_id_3 == atom_id_1) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_2 == atom_id_2) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_1 == atom_id_3) ){
	            aux = i;
	            break;
		}else if ( (topol_atoms_bond_angles_parameters[i].atom_id_1 == atom_id_2) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_2 == atom_id_3) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_3 == atom_id_1) ){
	            aux = i;
	            break;
		}else if ( (topol_atoms_bond_angles_parameters[i].atom_id_1 == atom_id_1) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_2 == atom_id_3) &&
				(topol_atoms_bond_angles_parameters[i].atom_id_3 == atom_id_2) ){
	            aux = i;
	            break;
		}
	}
	if (aux==-1){
		fatal_error("Index of bond angle is wrong when try to execute _get_index_bond_angles_parameters function!!!.\nCheck your parameter sequence.");
	}
	return aux;
}

static int _get_index_bond_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2){
	/*Receives two atoms which must be connected in topology.
	 * Returns the index of this bond.
	 * This function returns the index only. The value depends what function
	 * called it.
	 * */
	int aux=-1;
	int TAM = asize(topol_atoms_bond_parameters);
    for (int i=0; i< TAM;i++){
    	if ( /*atmC and atmCa like amtCa and atmC*/
    		( (topol_atoms_bond_parameters[i].atom_id_1 == atom_id_1) &&
    		  (topol_atoms_bond_parameters[i].atom_id_2 == atom_id_2) )    ||

      		( (topol_atoms_bond_parameters[i].atom_id_1 == atom_id_2) &&
      		  (topol_atoms_bond_parameters[i].atom_id_2 == atom_id_1) )
    	   ){
    		aux = i;
    		break;
    	}
    }
	if (aux==-1){
		fatal_error("Index of bond is wrong when try to execute _get_index_bond_parameters function!!!");
	}
	return aux;
}

float _get_bond_length_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2){
	int index = _get_index_bond_parameters(atom_id_1, atom_id_2);
	return topol_atoms_bond_parameters[index].bond_length;
}

float _get_bond_angle_from_atom_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2, type_atoms_t atom_id_3){
	int index = _get_index_bond_angles_parameters(atom_id_1, atom_id_2, atom_id_3);
	return topol_atoms_bond_angles_parameters[index].angle_value;
}

void check_number_atoms(const int *numatom, const top_global_t *top){
	/*Check the number of atoms between topology and other*/
	if (*numatom != top->numatom-1){
		fatal_error("Error was found when program was trying to build the topology... The number of atoms are different");
	}
}


void set_numatom_from_topol(int *numatom, const int *nr_atm_ff){
    *numatom = *numatom + *nr_atm_ff;
}

void set_num_bond_angle_from_topol(int *num_bond_angles,
		const int * number){
	*num_bond_angles = *num_bond_angles + *number;
}

void set_num_side_chains_from_topol(int *num_protein_side_chains,
		const int *number){
	/*
	 * Computes the number of side chains for protein.
	 * Example: ARG(R) has 5 side chains. GLY(G) there is not.
	 * The primary sequence RGGR has 10 side chains number.
	 * */
	*num_protein_side_chains = *num_protein_side_chains + *number;
}

void set_num_dihedral_angles(int *num_dihedral_angles, const int *number){
	/* Computes the number of dihedral angles of proteins. Those angle are
	 * formed by four atoms.
	*/
	*num_dihedral_angles = *num_dihedral_angles +  *number;
}

void set_has_his(boolean_t *has_his, const type_aminos_t *aminoid){
	/* Check if amino is HIS */
	if (*aminoid == aHIS){
		*has_his = btrue;
	}
}

void set_amino_from_topol(amino_t *prot_prin, int index_seq_prin, const topol_residues_t top_res[]){
	prot_prin[index_seq_prin].id = top_res->res_id;
	strcpy(prot_prin[index_seq_prin].aminoacido, top_res->res_name_1);
	strcpy(prot_prin[index_seq_prin].idl3, top_res->res_name);
	prot_prin[index_seq_prin].number_late = top_res->nr_side_chains;
	prot_prin[index_seq_prin].HP = top_res->HP;
	prot_prin[index_seq_prin].late = NULL;//Malloc(float, prot_prin[i].number_late);
}


static int _get_index_special_atom_parameters(type_atoms_t atm){
	/* Receives an atom pointer
	 * Returns its index
	 * if aux value is -1 it means that the atom is not a special atom
	 *
	 */
	   int TAM = asize(special_atom_parameters);
	   int aux = -1;
	   for (int i = 0; i < TAM;i++){
		   if (special_atom_parameters[i].sp_atom_id == atm){
			   aux = i;
			   break;
		   }
	   }
       return aux;
}

boolean_t _is_special_atom(type_atoms_t atm){
	/*Receives an atom
	 * Returns if it is a special atom or not*/
	boolean_t aux = bfalse;
	int index = _get_index_special_atom_parameters(atm);
	if ( index > -1){
		aux = btrue;
	}
	return aux;
}

float _get_bond_angle_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const int *atom_3, const top_global_t *top){
	/* Receives the atoms that make a bond angle.
	 * Returns the value of bond angle.
	 * It is used for build z matrix. See at _build_z_matrix function.
	*/
	float aux = -1;
	for (int i = 0; i < top->number_bond_angles;i++){
		if (top->top_global_res_atms_bond_angle[i].res_number == *res){
			if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_2) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_3)) {
				aux = top->top_global_res_atms_bond_angle[i].angle_value;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_1) &&
				(top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_2) &&
				(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_3) ){
				aux = top->top_global_res_atms_bond_angle[i].angle_value;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_1) &&
							(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_2) &&
							(top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_3)){
				aux = top->top_global_res_atms_bond_angle[i].angle_value;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_2) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_3) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_1)){
		         aux = top->top_global_res_atms_bond_angle[i].angle_value;
		         break;
	        }else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_3) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_2)){
		         aux = top->top_global_res_atms_bond_angle[i].angle_value;
		         break;
	        }else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_3) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_2)){
		         aux = top->top_global_res_atms_bond_angle[i].angle_value;
		         break;
	        }
		}
	}
	if (aux == -1){
		char msg[200];
		sprintf(msg,"Error when try to obtain bond angle of atoms %i %i %i from topology. The residue number is %i\n",
				*atom_1, *atom_2, *atom_3, *res);
		fatal_error(msg);
	}
	return aux;
}


float _get_bond_len_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const top_global_t *top){
	float aux = -1;
	for (int i = 0; i < top->numatom;i++){
		if (top->top_global_res_atms_bond[i].res_number == *res){
			if (
			( (top->top_global_res_atms_bond[i].atom_number1 == *atom_1) &&
			   (top->top_global_res_atms_bond[i].atom_number2 == *atom_2) ) ||
			  ( (top->top_global_res_atms_bond[i].atom_number1 == *atom_2) &&
				(top->top_global_res_atms_bond[i].atom_number2 == *atom_1) )   ){
			    aux = top->top_global_res_atms_bond[i].bond_value;
				break;
			}
		}
	}
	if (aux == -1){
		char msg[200];
		sprintf(msg,"Error when try to obtain bond length from topology. Atom numbers %i %i. Residue number %i\n",
				*atom_1, *atom_2, *res);
		fatal_error(msg);
		aux = 0;
	}
	return aux;

}

static float _get_diedhral_angle_from_protein_based_on_z_matrix(
		const z_matrix_global_t *z_matrix,const protein *prot,
		const int *i_z, const int *res){
	/*Receives z_matrix and protein
	 *
	 *Returns the value of diehdral angle based on z_matrix. Remember that
	 *z_matrix stores the kind of diehdral angle.The value of diedhral is
	 *obtained from protein
	 */
	if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_phi){
		return prot->residuo[*res].phi;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_psi){
		return prot->residuo[*res].psi;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_psi_){
		return prot->residuo[*res-1].psi;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_omega){
		return prot->residuo[*res-1].omega;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_chi1){
		return prot->residuo[*res].late[0];
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_chi2){
		return prot->residuo[*res].late[1];
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_chi3){
		return prot->residuo[*res].late[2];
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_chi4){
		return prot->residuo[*res].late[3];
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_0){
		return 0;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_180){
		return PI;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_90){
		return PI/2;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_117_){
		return -2.0420; // -2.0420 radian means -117 degrees
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_1){
		/*PI means 180 degrees*/
		if (prot->residuo[*res].psi < 0){
			return prot->residuo[*res].psi + PI;
		}else{
			return prot->residuo[*res].psi - PI;
		}
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_2){
		/*-2.1468 radian means -123 degrees*/
		return -2.1468;
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_3){
		int a = 10;
		int r = _get_int_random_number(&a);
		float RD = _get_double_gauss(&r);
		if (RD > 0){
			return RD - PI;
		}else{
			return RD + PI;
		}
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_trans_120){
		return PI + 2.0944; // 2.0944 = 120.0
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_trans_240){
		return PI + 4.1888; // 4.1888 = 240.0
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_trans_123_){
		return -2.1468; // -123.0 = -2.1468
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_4){
		if (prot->residuo[*res].late[1]  > 0){ // prot->residuo[*res].late[1] means chi2
			return prot->residuo[*res].late[1] - PI;
		}else{
			return prot->residuo[*res].late[1] + PI;
		}
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_5){
		if (prot->residuo[*res].late[2]  > 0){ // prot->residuo[*res].late[2] means chi3
			return prot->residuo[*res].late[2] - PI;
		}else{
			return prot->residuo[*res].late[2] + PI;
		}
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_6){
		return prot->residuo[*res].late[0] - 2.0944; // prot->residuo[*res].late[0] means chi1. 2.0944 = 120.0
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_7){
		return prot->residuo[*res].late[1] - 2.0944; // prot->residuo[*res].late[1] means chi2. 2.0944 = 120.0
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_29_6){
		return 0.5166; // 29.6 = 0.5166
	}else if (z_matrix->z_matrix_info[*i_z].tpAngle == angl_typ_dieh_37_4_){
		return -0.6056; // -37.4 = -0.6056
	}
	else{
		fatal_error("Error at _get_diedhral_angle_from_protein_based_on_z_matrix \n");
	}
}

void _pdbatoms2protein(protein *prot, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
	int atom_index_aux;
	type_aminos_t amino_id;
	for (int r = 0; r < top_global->numres; r++){
		prot->residuo[r].phi = _compute_phi(&r,pdb_atoms, top_global);
		prot->residuo[r].psi = _compute_psi(&r,pdb_atoms, top_global);
		//Atom from residue. Index atom means atom number -1
		atom_index_aux = top_global->top_global_dieh_psi[r].atom_number1-1;
		amino_id = _get_amino_id_3(top_global->
				top_global_atom[atom_index_aux].res_name);
		if (_has_side_chain(&amino_id) == btrue){
			for (int chi = 0; chi < _number_side_chains(&amino_id); chi++){
				prot->residuo[r].late[chi] = _compute_side_chains_angles(&r,
						&chi,pdb_atoms,top_global);
			}

		}
	}
}

float _compute_diehdral_angle(const own_vector_t *a1,
		const own_vector_t *a2,const own_vector_t *a3,	const own_vector_t *a4){
	/* Computes the dihedral angle	 */
	own_vector_t *P1, *P2, *M, *r1, *r2, *r3;
	double mod_P1, mod_P2;
	double W;

	P1 = Malloc(own_vector_t,1);
	P2 = Malloc(own_vector_t,1);
	M = Malloc(own_vector_t,1);
	r1 = Malloc(own_vector_t,1);
	r2 = Malloc(own_vector_t,1);
	r3 = Malloc(own_vector_t,1);

	//Computing distances
	sub_vector(r1,a1,a2);
	sub_vector(r2,a2,a3);
	sub_vector(r3,a3,a4);

	cross_product(P1,r1,r2);
	cross_product(P2,r2,r3);
	mod_P1 = mod_vector(P1);
	mod_P2 = mod_vector(P2);

	W = acos(scalar_prod(P1,P2)/(mod_P1*mod_P2));
	//Check if is necessary change the signal of W
	cross_product(M,P1,P2);
	if (scalar_prod(M,r2) < 0){
		W = W *(-1);
	}
	//Deallocating variables
	free(P1);
	free(P2);
	free(M);
	free(r1);
	free(r2);
	free(r3);

	return W;

}


static float _compute_phi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
	/* Receives the residue, the pdbatoms structure and the phi topology
	 * Returns the value of phi angle
	 * When residue is N-Terminal, there is not a phi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	if (top_global->top_global_dieh_phi[*r].atom_number1 == 0){//N-Terminal
		return 0;
	}
	const own_vector_t *a1, *a2, *a3, *a4;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float phi;
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number1-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number1-1].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number2-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number2-1].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number3-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number3-1].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number4-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number4-1].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;
	phi = _compute_diehdral_angle(a1,a2,a3,a4);
	return phi;
}

static float _compute_psi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
	/* Receives the residue, the pdbatoms structure and the psi topology
	 * Returns the value of psi angle
	 * When residue is C-Terminal, there is not a psi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	if (top_global->top_global_dieh_psi[*r].atom_number4 == 0){//C-Terminal
		return 0;
	}
	const own_vector_t *a1, *a2, *a3, *a4;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float psi;
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number1-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number1-1].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number2-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number2-1].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number3-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number3-1].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number4-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number4-1].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;

	psi = _compute_diehdral_angle(a1,a2,a3,a4);
	return psi;
}

static float  _compute_side_chains_angles(const int *r, const int *chi,
		const pdb_atom_t *pdb_atoms, const top_global_t *top_global) {
	/* Receives the residue, the index of chi, the pdbatoms structure and
	 * the side_chains topology
	 * The chi index is based on residue. Your value can be obtained in
	 * topol_residues_ff structure at topology_charmm27_parameters.h
	 *
	 * Returns the value of chi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	const own_vector_t *a1, *a2, *a3, *a4;
	int index;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float chi_angle;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(1,r,
			chi, top_global);
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(2,r,
			chi, top_global);
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(3,r,
			chi, top_global);
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(4,r,
			chi, top_global);
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;
	chi_angle = _compute_diehdral_angle(a1,a2,a3,a4);
	return chi_angle;

}

static int _get_atom_index_from_top_global_dihedral_side_chain_t(
		int atm_opt, const int *r, 	const int *chi,
		const top_global_t *top_global){
	/*
	 * Return the index of atom which is used to compute side chains. This
	 * function is used in _compute_side_chains_angles function.
	 * Must be (*r+1) because r is started with 0 and it is employed a lot of
	 * functions. Please look at _pdbatoms2protein function.
	 * Must be (*chi+1) because chi is started with 0. It is an index. Please
	 * look at _pdbatoms2protein function.
	 * The value -1 because index is atom number - 1
	 */
	int aux = -1;
	for (int i =0; i < top_global->number_protein_side_chains; i++){
		if (top_global->top_global_dieh_side_chains[i].res_number == (*r+1)){
			if (top_global->top_global_dieh_side_chains[i].chi == (*chi+1)){
				//Returns based on atom option
				if (atm_opt == 1){
					aux =  top_global->top_global_dieh_side_chains[i].atom_number1 - 1;
				}else if (atm_opt == 2){
					aux = top_global->top_global_dieh_side_chains[i].atom_number2 - 1;
				}else if (atm_opt == 3){
					aux = top_global->top_global_dieh_side_chains[i].atom_number3 - 1;
				}else if (atm_opt == 4){
					aux = top_global->top_global_dieh_side_chains[i].atom_number4 - 1;
				}
			}
		}
	}
	if (aux == -1){
		fatal_error("Index not found in _get_atom_index_from_top_global_dihedral_side_chain_t \n");
	}
	return aux;
}



void _protein2pdbatoms(pdb_atom_t *pdb_atoms, const protein *prot,
		const top_global_t *top_global, const z_matrix_global_t *z_matrix){
	/* Computes the conversion from dihedral space for cartesian space.
	 * res_number in topology starts with 1. However, the index must start in 0.
	 * Therefore, when need to obtain index value from res_number, you have to use
	 * res_number-1. Please, see
	 * top_global->top_global_atom[z_matrix->z_matrix_info[i_z].index_top].res_number-1
	 * */
	int res_aux;//Store the index of residue for getting the diehdral angle value
	int i_z = 0;
	int topol_index, num_atom_topol;
	double diedhral;
	own_vector_t* vr;
	float default_angle_alfa = PI/2;
	float default_angle_beta = PI;
	float alfa,beta;

	vr = Malloc(own_vector_t,1);

	topol_index = z_matrix->z_matrix_info[i_z].index_top;
	pdb_atoms[i_z].coord.x = 0;
	pdb_atoms[i_z].coord.y = 0;
	pdb_atoms[i_z].coord.z = 0;
	set_values_from_topology_2_pdbatoms(pdb_atoms, &i_z, top_global,
			z_matrix);

	i_z++;
	topol_index = z_matrix->z_matrix_info[i_z].index_top;
	pdb_atoms[i_z].coord.x = 0;
	pdb_atoms[i_z].coord.y = 0;
	pdb_atoms[i_z].coord.z = z_matrix->z_matrix_info[i_z].bond_len;
	set_values_from_topology_2_pdbatoms(pdb_atoms, &i_z, top_global,
			z_matrix);

	i_z++;
	//Third atom
	topol_index = z_matrix->z_matrix_info[i_z].index_top;
	pdb_atoms[i_z].coord.x = 0;
	alfa = z_matrix->z_matrix_info[i_z].bond_angle - default_angle_alfa;
	pdb_atoms[i_z].coord.y = z_matrix->z_matrix_info[i_z].bond_len*cos(alfa);
	beta = default_angle_beta - z_matrix->z_matrix_info[i_z].bond_angle;
	pdb_atoms[i_z].coord.z = z_matrix->z_matrix_info[i_z].bond_len*cos(beta) +
			z_matrix->z_matrix_info[i_z-1].bond_len;
	set_values_from_topology_2_pdbatoms(pdb_atoms, &i_z, top_global,
			z_matrix);
	for (i_z = 3; i_z < z_matrix->num_elements; i_z++){
		//Remember that res_number is started with 1.
		topol_index = z_matrix->z_matrix_info[i_z].index_top;
		res_aux = top_global->top_global_atom[topol_index].res_number-1;
		diedhral = _get_diedhral_angle_from_protein_based_on_z_matrix(z_matrix,
				prot,&i_z,&res_aux);
		NeRF(vr,&z_matrix->z_matrix_info[i_z].bond_len_2,
				&z_matrix->z_matrix_info[i_z].bond_len,
				&z_matrix->z_matrix_info[i_z].bond_angle,
				&diedhral,
				&pdb_atoms[z_matrix->z_matrix_info[i_z].dihedral_connect-1].coord,//A
				&pdb_atoms[z_matrix->z_matrix_info[i_z].atom_angle-1].coord, //B
				&pdb_atoms[z_matrix->z_matrix_info[i_z].atom_connected-1].coord); //C
		pdb_atoms[i_z].coord.x = vr->x;
		pdb_atoms[i_z].coord.y = vr->y;
		pdb_atoms[i_z].coord.z = vr->z;
		set_values_from_topology_2_pdbatoms(pdb_atoms, &i_z, top_global,
				z_matrix);
	}
	free(vr);
}

type_aminos_t _get_amino_id_from_res_id(const int *res_id,
		const top_global_t *top){
	/*Receives res_id which means residue number in protein
	 * Returns amino_id
	 * index_amino means res_id - 1. Because res_id is started 1.
	 */
	int index_first_atom;
	int index_amino;

	index_amino = *res_id -1;
	if ( index_amino > top->numres){
		char msg[1024];
		sprintf(msg,"res_id %d is not a valid residue number in Topology %d \n ",*res_id, top->numres);
		fatal_error(msg);
	}
	index_first_atom = top->top_global_res_atm[index_amino].atom_first;
	return top->top_global_atom[index_first_atom].amino_id;
}
boolean_t _has_side_chain(const type_aminos_t *amino_id){
	/*Check residue has or not side chains*/
	const topol_residues_t *residue_ff= get_residue_from_topol(*amino_id);
	if (residue_ff->nr_side_chains > 0){
		return btrue;
	}else{
		return bfalse;
	}
}

int _number_side_chains(const type_aminos_t *amino_id){
	/*Returns the number of side chains per residues*/
	const topol_residues_t *residue_ff= get_residue_from_topol(*amino_id);
	return residue_ff->nr_side_chains;
}


static void set_values_from_topology_2_pdbatoms(pdb_atom_t *pdb_atoms,
		const int *i_z,	const top_global_t *top_global,
		const z_matrix_global_t *z_matrix ){
/*Receives pdb_atoms, index of pdbatoms and z_matix, global topology and z_matrx
 * Copy the value from global topology to pdbatoms structure. These values are
 * based on z_matrix.
 * It is assumed that z_matrix and pdbatoms have the same value of index.
 */
	set_pdb_atom_generic_information(&pdb_atoms[*i_z],
			top_global->top_global_atom[z_matrix->z_matrix_info[*i_z].index_top].atom_name,
			top_global->top_global_atom[z_matrix->z_matrix_info[*i_z].index_top].res_name,
			NULL,
			&top_global->top_global_atom[z_matrix->z_matrix_info[*i_z].index_top].res_number,
			&top_global->top_global_atom[z_matrix->z_matrix_info[*i_z].index_top].atom_number
			);
}
type_aminos_t _get_amino_id_3(char *c){
/*Receives an amino (char) and returns your id*/
		type_aminos_t amino_id;
		char msg[30];
		if ( (strcmp(c,"GLY") == 0) ){
			amino_id =  aGLY;
		}else if ( (strcmp(c,"ARG") == 0) ){
			amino_id =  aARG;
		}else if ( (strcmp(c,"ALA") == 0) ){
			amino_id =  aALA;
		}else if ( (strcmp(c,"SER") == 0) ){
			amino_id =  aSER;
		}else if ( (strcmp(c,"THR") == 0) ){
			amino_id =  aTHR;
		}else if ( (strcmp(c,"CYS") == 0) ){
			amino_id =  aCYS;
		}else if ( (strcmp(c,"VAL") == 0) ){
			amino_id =  aVAL;
		}else if ( (strcmp(c,"LEU") == 0) ){
			amino_id =  aLEU;
		}else if ( (strcmp(c,"ILE") == 0) ){
			amino_id =  aILE;
		}else if ( (strcmp(c,"MET") == 0) ){
			amino_id =  aMET;
		}else if ( (strcmp(c,"PRO") == 0) ){
			amino_id =  aPRO;
		}else if ( (strcmp(c,"PHE") == 0) ){
			amino_id =  aPHE;
		}else if ( (strcmp(c,"TYR") == 0) ){
			amino_id =  aTYR;
		}else if ( (strcmp(c,"TRP") == 0) ){
			amino_id =  aTRP;
		}else if ( (strcmp(c,"ASP") == 0) ){
			amino_id =  aASP;
		}else if ( (strcmp(c,"GLU") == 0) ){
			amino_id =  aGLU;
		}else if ( (strcmp(c,"ASN") == 0) ){
			amino_id =  aASN;
		}else if ( (strcmp(c,"GLN") == 0) ){
			amino_id =  aGLN;
		}else if ( (strcmp(c,"HIS") == 0) ){
			amino_id =  aHIS;
		}else if ( (strcmp(c,"LYS") == 0) ){
			amino_id =  aLYS;
		}else{
			if (strcmp(c,"") == 0){
				sprintf(msg,"Amino not found, because amino variable is empty. Check it. \n");
			}else{
				sprintf(msg,"%s Amino not found. Check it. \n",c);
			}
			fatal_error(msg);
		}
		return amino_id;
}


void _check_pdb_atoms_topology(const pdb_atom_t *pdb_atoms,
		const top_global_t *top){
	/*This function check pdb_atoms and global topology */
	char message[300];
	for (int a = 0; a < top->numatom;a++){
		//Checking residue numbers
		if (pdb_atoms[a].resnum != top->top_global_atom[a].res_number){
			sprintf(message,"The residue number is not equal to atom %i : "
					"In PDB %i Topol %i \n",a+1,pdb_atoms[a].resnum,
					top->top_global_atom[a].res_number);
			fatal_error(message);
		}
		//Checking residue names
		if ( is_equal(pdb_atoms[a].resname,
				top->top_global_atom[a].res_name) == bfalse ) {
			sprintf(message,"The residue name is not equal to atom %i : "
					"In PDB %s Topol %s \n",a+1,pdb_atoms[a].resname,
					top->top_global_atom[a].res_name);
			fatal_error(message);
		}
		//Checking atom numbers
		if (pdb_atoms[a].atmnumber != top->top_global_atom[a].atom_number){
			sprintf(message,"The atom number is not equal to atom %i : "
					"In PDB %i Topol %i \n",a+1,pdb_atoms[a].atmnumber,
					top->top_global_atom[a].atom_number);
			fatal_error(message);
		}
		//Checking atom id
		if (pdb_atoms[a].atomid != top->top_global_atom[a].atom_id){
			sprintf(message,"The atom id (IDentification) is not equal to atom %i in PDB and Topology. "
					"It means that these atoms can not be the same type. \n",a+1);
			fatal_error(message);
		}
	}
}

int _get_four_atom_number_and_type_for_dihedral_angle(
		type_dihedral_angles_t *tp_dihedral_angle,
		const int *res,
		const int *atom_1, const int *atom_2, const int *atom_3,
		const top_global_t *top){
	/* Receives the residue number and three atoms of that residue.
	 * Returns not only the four atom number that makes a dihedral angle with atom_1,
	 * atom_2 and atom_3, but also type of dihedral angle
	 */
	int aux = -1;
	for (int i = 0; i < top->number_dihedral_angles_type;i++){
		if (top->top_global_dihedral_angles_type[i].res_number == *res){
			if ( ((top->top_global_dihedral_angles_type[i].atom_number1 == *atom_1) &&
			   (top->top_global_dihedral_angles_type[i].atom_number2 == *atom_2)) &&
			   (top->top_global_dihedral_angles_type[i].atom_number3 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number4;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number1 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number2 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number4 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number3;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number1 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number3 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number4 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number2;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number2 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number3 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number4 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number1;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number2 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number4 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number1 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number3;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number2 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number3 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number1 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number4;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number3 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number4 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number2 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number1;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number3 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number1 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number2 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number4;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number4 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number1 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number2 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number3;
				break;
			}else if ( ((top->top_global_dihedral_angles_type[i].atom_number4 == *atom_1) &&
					   (top->top_global_dihedral_angles_type[i].atom_number3 == *atom_2)) &&
					   (top->top_global_dihedral_angles_type[i].atom_number1 == *atom_3) ){
				aux = top->top_global_dihedral_angles_type[i].atom_number2;
				break;
			}
		}
		// Setting type of dihedral angle
		*tp_dihedral_angle = top->top_global_dihedral_angles_type[i].type_dihedral;
		/* printf command is used when try to discover a bug
		printf("Testing type of dihedral angle %d %d \n",
				*tp_dihedral_angle,
				top->top_global_dihedral_angles_type[i].type_dihedral);
		printf("_get_four_atom_number_and_type_for_dihedral_angle enum %i %i %i %i %i\n",top->top_global_dihedral_angles_type[i].atom_number1,
				top->top_global_dihedral_angles_type[i].atom_number2,
				top->top_global_dihedral_angles_type[i].atom_number3,
				top->top_global_dihedral_angles_type[i].atom_number4,
				top->top_global_dihedral_angles_type[i].res_number);
	   */
	}
	if (aux == -1){
		char msg[200];
		sprintf(msg,"Error when try to obtain the four atom for diehdral angle. The atoms are %i %i %i from topology. The residue number is %i\n",
				*atom_1, *atom_2,*atom_3, *res);
		fatal_error(msg);
	}
	return aux;
}

int _get_third_atom_number_for_bond_angle(const int *res,
		const int *atom_1, const int *atom_2, const top_global_t *top){
	/* Receives the residue number and two atoms of that residue.
	 * Returns the third atom number that bond (atom_1 and atom_2). This third
	 * atom is necessary because these tree atoms make a bond angle.
	 */
	int aux = -1;
	for (int i = 0; i < top->number_bond_angles;i++){
		if (top->top_global_res_atms_bond_angle[i].res_number == *res){
			if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_2) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number2;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_2) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_1) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number2;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_2) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number1;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_2) &&
					(top->top_global_res_atms_bond_angle[i].atom_number3 == *atom_1) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number1;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_1) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_2) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number3;
				break;
			}else if ( (top->top_global_res_atms_bond_angle[i].atom_number1 == *atom_2) &&
					(top->top_global_res_atms_bond_angle[i].atom_number2 == *atom_1) ){
				aux = top->top_global_res_atms_bond_angle[i].atom_number3;
				break;
			}

		}
	}
	if (aux == -1){
		char msg[200];
		sprintf(msg,"Error when try to obtain the third atom for bond angle. The atoms are %i %i from topology. The residue number is %i\n",
				*atom_1, *atom_2, *res);
		fatal_error(msg);
	}
	return aux;
}


int _get_number_atom_bond(const int *res, const int *atom_1,
		const top_global_t *top){
	/*Receives a residue number, a atom number of residue and the topology.
	 * Returns the atom number which is connected with atom number informed.
	 * This function is used at _build_z_matrix function to get the atom number
	 */
	int aux = -1;
	for (int i = 0; i < top->numatom; i++){
		if (top->top_global_res_atms_bond[i].res_number == *res){
			if (top->top_global_res_atms_bond[i].atom_number2 == *atom_1){
				aux = top->top_global_res_atms_bond[i].atom_number1;
				break;
			}
		}
	}
	if (aux == -1){
		char msg[200];
		sprintf(msg,"Error when try to obtain the atom bond. The atom is %i . The residue number is %i\n",
				*atom_1, *res);
		fatal_error(msg);
	}
	return aux;
}

void _type_of_diedhral_angle2str(char *str,
		const type_dihedral_angles_t *type_dihedral){
	/* Converts type_dihedral to string
	 * type_dihedral_angles_t is defined at enums.h
	*/
	if ( *type_dihedral == angl_phi){
		strcpy(str,"PHI");
	}else if (*type_dihedral == angl_psi){
		strcpy(str,"PSI");
	}else if (*type_dihedral == angl_chi1){
		strcpy(str,"CHI1");
	}else if (*type_dihedral == angl_chi2){
		strcpy(str,"CHI2");
	}else if (*type_dihedral == angl_chi3){
		strcpy(str,"CHI3");
	}else if (*type_dihedral == angl_chi4){
		strcpy(str,"CHI4");
	}else if (*type_dihedral == angl_typ_dieh_1){
		strcpy(str,"W1");
	}else if (*type_dihedral == angl_typ_dieh_2){
		strcpy(str,"W2");
	}else if (*type_dihedral == angl_typ_dieh_3){
		strcpy(str,"W3");
	}else if (*type_dihedral == angl_typ_dieh_180){
		strcpy(str,"dieh_180");
	}else if (*type_dihedral == angl_typ_dieh_0){
		strcpy(str,"dieh_0");
	}else if (*type_dihedral == angl_typ_dieh_90){
		strcpy(str,"dieh_90");
	}else if (*type_dihedral == angl_typ_dieh_117_){
		strcpy(str,"dieh_117_");
	}else if (*type_dihedral == angl_typ_dieh_trans_120){
		strcpy(str,"dieh_trans_120");
	}else if (*type_dihedral == angl_typ_dieh_trans_240){
		strcpy(str,"dieh_trans_240");
	}else if (*type_dihedral == angl_typ_dieh_trans_123_){
		strcpy(str,"angl_typ_dieh_trans_123_");
	}else{
		fatal_error("Type of dihedral not found. Please look at type_dihedral_angles_t enum.\n");
	}
}

int _get_number_atoms_from_res(const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top){
	/*This function returns how many atoms the residue has*/
	const topol_residues_t *residue_ff= NULL;
	if (*res_id == 1){//N-Terminal
		residue_ff = get_residue_from_topol_N_Terminal(*amino_id);
	}else if ( (*res_id > 1) && (*res_id < top->numres) ){
		residue_ff = get_residue_from_topol(*amino_id);
	}else{//C-Terminal
		residue_ff = get_residue_from_topol_C_Terminal(*amino_id);
	}
	return residue_ff->nr_atoms;
}

int _get_number_atoms_from_res_C_Terminal(const type_aminos_t *amino_id,
		const top_global_t *top){
	/*This function returns how many atoms the residue has at C Terminal*/
	const topol_residues_t *residue_ff= NULL;
	residue_ff = get_residue_from_topol_C_Terminal(*amino_id);
	return residue_ff->nr_atoms;
}

const topol_residues_t* _get_topol_residues_t_from_res(const int *res_id,
		const type_aminos_t *amino_id, const top_global_t *top){
	/*This function returns the topol_residues_t struct of a specific residue*/
	const topol_residues_t *residue_ff= NULL;
	if (*res_id == 1){//N-Terminal
		residue_ff = get_residue_from_topol_N_Terminal(*amino_id);
	}else if ( (*res_id > 1) && (*res_id < top->numres) ){
		residue_ff = get_residue_from_topol(*amino_id);
	}else{//C-Terminal
		residue_ff = get_residue_from_topol_C_Terminal(*amino_id);
	}
	if (residue_ff == NULL){
		fatal_error("An error was found when try to execute the _get_topol_residues_t_from_res function \n");
	}
	return residue_ff;
}

const topol_residue_atoms_t* _get_topol_residue_atoms_t_from_res(int *num_atom,
		const int *res_id, 	const type_aminos_t *amino_id,
		const top_global_t *top){
	const topol_residue_atoms_t* atom_ff=NULL;
	const topol_residues_t *residue_ff= NULL;

	residue_ff = _get_topol_residues_t_from_res(res_id, amino_id, top);
	*num_atom = residue_ff->nr_atoms;
	atom_ff = residue_ff->residue_atoms;
	if (atom_ff == NULL){
		fatal_error("An error was found when try to execute the _get_topol_residue_atoms_t_from_res function \n");
	}
	return atom_ff;
}
