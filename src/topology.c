#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "topology.h"
#include "functions.h"
#include "consts.h"
#include "messages.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"
#include "z_matrix_types.h"

#define TAM_BLOCO_PRM 80


top_global_t *allocateTop_Global(const int *numatom, const int *numres,
		const int *num_bond_angles, const int *num_side_chains,
		const int *number_dihedrals_type){
	top_global_t *top_aux;

	top_aux = Malloc(top_global_t,1);
	top_aux->numatom = *numatom;
	top_aux->numres = *numres;
	top_aux->number_bond_angles = *num_bond_angles;
	top_aux->number_protein_side_chains = *num_side_chains;
	top_aux->number_dihedral_angles_type = *number_dihedrals_type;
	top_aux->top_global_atom = Malloc(top_global_atom_t,*numatom);
	top_aux->top_global_dieh_phi = Malloc(top_global_dihedral_t,*numres);
	top_aux->top_global_dieh_psi = Malloc(top_global_dihedral_t,*numres);
	top_aux->top_global_res_atm = Malloc(top_global_res_atm_t,*numres);
	top_aux->top_global_res_atms_bond = Malloc(top_global_res_atms_bond_t,*numatom);
	top_aux->top_global_res_atms_bond_angle = Malloc(top_global_res_atms_bond_angle_t,*num_bond_angles);
	//top_aux->top_global_dieh_phi = Malloc(top_global_dihedral_t,*numres);
	//top_aux->top_global_dieh_psi = Malloc(top_global_dihedral_t,*numres);
	top_aux->top_global_dieh_side_chains = Malloc(top_global_dihedral_side_chain_t,*num_side_chains);
	top_aux->top_global_dihedral_angles_type = Malloc(top_global_dihedral_angles_type_t, *number_dihedrals_type);
	top_aux->top_global_dieh_omega = Malloc(top_global_dihedral_t,*numres);	
	return top_aux;
}
void  desAllocateTop_Global(top_global_t *top_global){
	//Falta criar um correto desallocate
	free(top_global->top_global_atom);
	free(top_global->top_global_dieh_phi);
	free(top_global->top_global_dieh_psi);
	free(top_global->top_global_res_atm);
	free(top_global->top_global_res_atms_bond);
	free(top_global->top_global_res_atms_bond_angle);
	//free(top_global->top_global_dieh_phi);
	//free(top_global->top_global_dieh_psi);
	free(top_global->top_global_dieh_side_chains);
	free(top_global->top_global_dihedral_angles_type);
	free(top_global->top_global_dieh_omega);
	free(top_global);
}

protein_backbone_t *allocateProtein_backbone(const int *numres){
	protein_backbone_t *aux;
	aux = Malloc(protein_backbone_t,*numres);
	return aux;
}

void  desAllocateProtein_backbone(protein_backbone_t *protein_backbone){
	//Falta criar um correto desallocate
	free(protein_backbone);
}

static void build_sections_atom_and_residue_atoms(const amino_t *sequence_primary, top_global_t *top){
	int index_atoms_topol_global=-1;
	type_aminos_t amino_id;
	for (int r=0; r < top->numres ;r++){
		if (r==0){
			top->top_global_res_atm[r].atom_first = 1;
		}else{
			top->top_global_res_atm[r].atom_first = top->top_global_atom[index_atoms_topol_global].atom_number +1;
		}
		amino_id = sequence_primary[r].id;
		if (r==0){//N-Terminal
			set_atoms_from_topol_ff_N_Terminal(top, &index_atoms_topol_global,
					amino_id, &r);
		}else if ( (r > 0) && (r < top->numres-1) ){
			set_atoms_from_topol_ff(top, &index_atoms_topol_global, amino_id, &r);
		}else{ //C-Terminal
			set_atoms_from_topol_ff_C_Terminal(top, &index_atoms_topol_global,
					amino_id, &r);
		}
		top->top_global_res_atm[r].res_number = top->top_global_atom[index_atoms_topol_global].res_number;
		top->top_global_res_atm[r].atom_last  = top->top_global_atom[index_atoms_topol_global].atom_number;
	}
	check_number_atoms(&index_atoms_topol_global,top);
}

static void set_residue_atoms_bond_from_topol_ff(
		int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	const topol_residues_t *residue_ff= get_residue_from_topol(amino_id);
	int TAM = residue_ff->nr_atoms;
	type_atoms_t aux_atom_1, aux_atom_2;
    int res_id_aux_1, res_id_aux_2;

	for (int a=0; a < TAM; a++){
		*last_index_atm_bond = *last_index_atm_bond +1;
		top->top_global_res_atms_bond[*last_index_atm_bond].res_number = *res_id;
		set_residue_atom_values(&aux_atom_1,&res_id_aux_1,residue_ff->residue_atoms_bond[a].atom_id_1,res_id);
		if (res_id_aux_1 == 0){
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number1 = 0;
		}else{
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
					aux_atom_1,top);
		}
		set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_bond[a].atom_id_2,res_id);
		if (res_id_aux_2 == 0){
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = 0;
		}else{
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
					aux_atom_2,top);
		}
		top->top_global_res_atms_bond[*last_index_atm_bond].bond_value = get_bond_length_from_atom_parameters(
				                                                      residue_ff->residue_atoms_bond[a].atom_id_1,
				                                                      residue_ff->residue_atoms_bond[a].atom_id_2);
	}
}

static void set_residue_atoms_bond_from_topol_ff_N_Terminal(
		int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	const topol_residues_t *residue_ff= get_residue_from_topol_N_Terminal(amino_id);
	int TAM = residue_ff->nr_atoms;
	type_atoms_t aux_atom_2;
    int res_id_aux_2;
	for (int a=0; a < TAM; a++){
		*last_index_atm_bond = *last_index_atm_bond +1;
		top->top_global_res_atms_bond[*last_index_atm_bond].res_number = *res_id;
		top->top_global_res_atms_bond[*last_index_atm_bond].atom_number1 = get_num_atom_from_topol(res_id,residue_ff->residue_atoms_bond[a].atom_id_1,top);
		set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_bond[a].atom_id_2,res_id);
		if (res_id_aux_2 == 0){
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = 0;
		}else{
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
					aux_atom_2,top);
		}
		top->top_global_res_atms_bond[*last_index_atm_bond].bond_value = get_bond_length_from_atom_parameters(
				                                                      residue_ff->residue_atoms_bond[a].atom_id_1,
				                                                      residue_ff->residue_atoms_bond[a].atom_id_2);
	}
}

static void set_residue_atoms_bond_from_topol_ff_C_Terminal(
		int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	const topol_residues_t *residue_ff= get_residue_from_topol_C_Terminal(amino_id);
	int TAM = residue_ff->nr_atoms;
	type_atoms_t aux_atom_1, aux_atom_2;
    int res_id_aux_1, res_id_aux_2;
	for (int a=0; a < TAM; a++){
		*last_index_atm_bond = *last_index_atm_bond +1;
		top->top_global_res_atms_bond[*last_index_atm_bond].res_number = *res_id;
		set_residue_atom_values(&aux_atom_1,&res_id_aux_1,residue_ff->residue_atoms_bond[a].atom_id_1,res_id);
		if (res_id_aux_1 == 0){
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number1 = 0;
		}else{
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
					aux_atom_1,top);
		}
		set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_bond[a].atom_id_2,res_id);
		if (res_id_aux_2 == 0){
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = 0;
		}else{
			top->top_global_res_atms_bond[*last_index_atm_bond].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
					aux_atom_2,top);
		}
		top->top_global_res_atms_bond[*last_index_atm_bond].bond_value = get_bond_length_from_atom_parameters(
				                                                      residue_ff->residue_atoms_bond[a].atom_id_1,
				                                                      residue_ff->residue_atoms_bond[a].atom_id_2);
	}
}

static void set_residue_atoms_angle_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	/*
	 * The value 0 for atom_number2 and atom_number3 mean that there are not these atoms.
	 * In N_terminal case, the res_id -1 is impossible. Therefore, C_terminal case res_id + 1 is
	 * impossible too.
	 * Remember atom_number1 is the atom of reference residue. Hence, it exists always.
	 * */
	type_atoms_t aux_atom_2, aux_atom_3;
    int res_id_aux_2, res_id_aux_3;
	const topol_residues_t *residue_ff= NULL;
	if (*res_id == 1){//N-Terminal
		residue_ff = get_residue_from_topol_N_Terminal(amino_id);
	}else if ( (*res_id > 1) && (*res_id < top->numres) ){
		residue_ff = get_residue_from_topol(amino_id);
	}else{//C-Terminal
		residue_ff = get_residue_from_topol_C_Terminal(amino_id);
	}
	for (int a=0; a < residue_ff->number_bond_angle; a++){
		*last_index_atm_angle = *last_index_atm_angle +1;
		top->top_global_res_atms_bond_angle[*last_index_atm_angle].res_number = *res_id;
		top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(res_id,
						                  residue_ff->residue_atoms_bonds_angles[a].atom_id_1,top);
		set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_bonds_angles[a].atom_id_2,res_id);
		set_residue_atom_values(&aux_atom_3,&res_id_aux_3,residue_ff->residue_atoms_bonds_angles[a].atom_id_3,res_id);
		if (*res_id == 1){ //N_terminal
			if (res_id_aux_2 == 0){
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number2 = 0;
			}else{
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
									aux_atom_2,top);
			}
			if (res_id_aux_3 == 0){
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number3 = 0;
			}else{
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
						aux_atom_3,top);
			}
		}else if (*res_id <= top->numres -1){ //internal residues
			top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
					aux_atom_2,top);
			top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		} else if (*res_id == top->numres){ // C_terminal
			if (res_id_aux_2 > top->numres){
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number2 = 0;
			}else{
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
									aux_atom_2,top);
			}
			if (res_id_aux_3 > top->numres){
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number3 = 0;
			}else{
				top->top_global_res_atms_bond_angle[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
						aux_atom_3,top);
			}
		}
		top->top_global_res_atms_bond_angle[*last_index_atm_angle].angle_value = get_angle_from_atom_parameters(
				                 residue_ff->residue_atoms_bonds_angles[a].atom_id_1,
				                 residue_ff->residue_atoms_bonds_angles[a].atom_id_2,
				                 residue_ff->residue_atoms_bonds_angles[a].atom_id_3);
	}
}

static void set_residue_atoms_dihedral_angle_type_from_topol_ff(int *last_index_atm_angle,
		const int *res_id, type_aminos_t amino_id,	const top_global_t *top){
	/* This function stores the number of atoms for dihedral angle type
	 * The value 0 for atom_number1, atom_number2,  atom_number3 and atom_number4 mean that there are not these atoms.
	 * In N_terminal case, the res_id -1 is impossible. Therefore, C_terminal case res_id + 1 is
	 * impossible too.
     */
	type_atoms_t aux_atom_1, aux_atom_2, aux_atom_3, aux_atom_4;
    int res_id_aux_1, res_id_aux_2, res_id_aux_3, res_id_aux_4;
    const topol_residues_t *residue_ff= NULL;

	if (*res_id == 1){//N_terminal
		residue_ff = get_residue_from_topol_N_Terminal(amino_id);
	}else if (*res_id <= top->numres -1){ //internal residues
		residue_ff = get_residue_from_topol(amino_id);
	}else if (*res_id == top->numres){// C_terminal
		residue_ff = get_residue_from_topol_C_Terminal(amino_id);
	}
	for (int i = 0; i < residue_ff->number_dihedral_angle_type;i++){
		*last_index_atm_angle = *last_index_atm_angle +1;
		top->top_global_dihedral_angles_type[*last_index_atm_angle].res_number = *res_id;
		top->top_global_dihedral_angles_type[*last_index_atm_angle].type_dihedral = residue_ff->residue_atoms_dihedral_angles_type[i].type_dihedral;
		set_residue_atom_values(&aux_atom_1,&res_id_aux_1,
				residue_ff->residue_atoms_dihedral_angles_type[i].atom_id_1 ,res_id);
		set_residue_atom_values(&aux_atom_2,&res_id_aux_2,
				residue_ff->residue_atoms_dihedral_angles_type[i].atom_id_2 ,res_id);
		set_residue_atom_values(&aux_atom_3,&res_id_aux_3,
				residue_ff->residue_atoms_dihedral_angles_type[i].atom_id_3 ,res_id);
		set_residue_atom_values(&aux_atom_4,&res_id_aux_4,
				residue_ff->residue_atoms_dihedral_angles_type[i].atom_id_4 ,res_id);
		if (res_id_aux_1 == 0){
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number2 = 0;

		}else{
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
                   					aux_atom_2,top);
		}
		if (res_id_aux_3 == 0){
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number3 = 0;

		}else{
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
                   					aux_atom_3,top);
		}
		if (res_id_aux_4 == 0){
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dihedral_angles_type[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
                   					aux_atom_4,top);
		}

	}
}

static void set_residue_atoms_dihedral_phi_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	/*
	 * The value 0 for atom_number1, atom_number2,  atom_number3 and atom_number4 mean that there are not these atoms.
	 * In N_terminal case, the res_id -1 is impossible. Therefore, C_terminal case res_id + 1 is
	 * impossible too.
	 * */

	type_atoms_t aux_atom_1, aux_atom_2, aux_atom_3, aux_atom_4;
    int res_id_aux_1, res_id_aux_2, res_id_aux_3, res_id_aux_4;
    const topol_residues_t *residue_ff = NULL;

    if (*res_id == 1){ //N_terminal
    	residue_ff= get_residue_from_topol_N_Terminal(amino_id);
    }else if (*res_id <= top->numres -1){ //internal residues{
    	residue_ff= get_residue_from_topol(amino_id);
    }else if (*res_id == top->numres){ // C_terminal{
    	residue_ff= get_residue_from_topol_C_Terminal(amino_id);
    }

	*last_index_atm_angle = *last_index_atm_angle +1;
	top->top_global_dieh_phi[*last_index_atm_angle].res_number = *res_id;

	set_residue_atom_values(&aux_atom_1,&res_id_aux_1,residue_ff->residue_atoms_phi[0].atom_id_1,res_id);
	set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_phi[0].atom_id_2,res_id);
	set_residue_atom_values(&aux_atom_3,&res_id_aux_3,residue_ff->residue_atoms_phi[0].atom_id_3,res_id);
	set_residue_atom_values(&aux_atom_4,&res_id_aux_4,residue_ff->residue_atoms_phi[0].atom_id_4,res_id);

	if (*res_id == 1){ //N_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}

	}else if (*res_id <= top->numres -1){ //internal residues
		top->top_global_dieh_phi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                aux_atom_1,top);
		top->top_global_dieh_phi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
				aux_atom_2,top);
		top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
				aux_atom_3,top);
		top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
				aux_atom_4,top);
	} else if (*res_id == top->numres){ // C_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 == 0){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}
	}
}

static void set_residue_atoms_dihedral_side_chains_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	type_atoms_t aux_atom_1, aux_atom_2, aux_atom_3, aux_atom_4;

	const topol_residues_t *residue_ff= get_residue_from_topol(amino_id);

    if (residue_ff->residue_atoms_side_chains != NULL){
    	for (int chi=1; chi <= residue_ff->nr_side_chains; chi++){
    		*last_index_atm_angle = *last_index_atm_angle +1;
    		top->top_global_dieh_side_chains[*last_index_atm_angle].res_number = *res_id;
        	top->top_global_dieh_side_chains[*last_index_atm_angle].chi = chi;
        	top->top_global_dieh_side_chains[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(res_id,
        			residue_ff->residue_atoms_side_chains[chi-1].atom_id_1,top);
        	top->top_global_dieh_side_chains[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(res_id,
        			residue_ff->residue_atoms_side_chains[chi-1].atom_id_2,top);
        	top->top_global_dieh_side_chains[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(res_id,
        			residue_ff->residue_atoms_side_chains[chi-1].atom_id_3,top);
        	top->top_global_dieh_side_chains[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(res_id,
        			residue_ff->residue_atoms_side_chains[chi-1].atom_id_4,top);
    	}
    }
}


static void set_residue_atoms_dihedral_psi_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	/*
	 * The value 0 for atom_number1, atom_number2,  atom_number3 and atom_number4 mean that there are not these atoms.
	 * In N_terminal case, the res_id -1 is impossible. Therefore, C_terminal case res_id + 1 is
	 * impossible too.
	 * */

	type_atoms_t aux_atom_1, aux_atom_2, aux_atom_3, aux_atom_4;
    int res_id_aux_1, res_id_aux_2, res_id_aux_3, res_id_aux_4;
	*last_index_atm_angle = *last_index_atm_angle +1;
	top->top_global_dieh_psi[*last_index_atm_angle].res_number = *res_id;

	const topol_residues_t *residue_ff= NULL;
	if (*res_id == 1){//N-Terminal
		residue_ff = get_residue_from_topol_N_Terminal(amino_id);
	}else if ( (*res_id > 1) && (*res_id < top->numres) ){
		residue_ff = get_residue_from_topol(amino_id);
	}else{//C-Terminal
		residue_ff = get_residue_from_topol_C_Terminal(amino_id);
	}


	set_residue_atom_values(&aux_atom_1,&res_id_aux_1,residue_ff->residue_atoms_psi[0].atom_id_1,res_id);
	set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_psi[0].atom_id_2,res_id);
	set_residue_atom_values(&aux_atom_3,&res_id_aux_3,residue_ff->residue_atoms_psi[0].atom_id_3,res_id);
	set_residue_atom_values(&aux_atom_4,&res_id_aux_4,residue_ff->residue_atoms_psi[0].atom_id_4,res_id);

	if (*res_id == 1){ //N_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 == 0){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 == 0){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}

	}else if (*res_id <= top->numres -1){ //internal residues
		top->top_global_dieh_psi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                aux_atom_1,top);
		top->top_global_dieh_psi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
				aux_atom_2,top);
		top->top_global_dieh_psi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
				aux_atom_3,top);
		top->top_global_dieh_psi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
				aux_atom_4,top);
	} else if (*res_id == top->numres){ // C_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 > top->numres){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 > top->numres){
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 > top->numres){
			top->top_global_dieh_psi[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_phi[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}
	}
}

static void build_section_bonds(const amino_t *primary_sequence,top_global_t *top){
	/*Receives a primary sequence and topology. Build the bond section of protein
	 *
	 * 1) Have to start with 1, because residue number is started with 1.
	 * 2) r-1 means the index of sequence_primary starts in 0.
	 * 3) index_atm_bond starts with -1 because I increase this value in
	 * set_residue_atoms_bond_from_topol_ff
	 */
	type_aminos_t amino_id;
	int index_atm_bond = -1;
	for (int r=1; r<=top->numres;r++){
        amino_id = primary_sequence[r-1].id;
        if (r == 1){//N-Terminal
        	set_residue_atoms_bond_from_topol_ff_N_Terminal(&index_atm_bond,&r,
        			amino_id, top);
        }else if ((r > 1) && (r<=top->numres-1) ){
        	set_residue_atoms_bond_from_topol_ff(&index_atm_bond,&r,amino_id,
        			top);
        }else{//C-Terminal
        	set_residue_atoms_bond_from_topol_ff_C_Terminal(&index_atm_bond,&r,
        			amino_id, top);
        }
	}
}

static void build_section_angle_bond(const amino_t *primary_sequence,
		top_global_t *top){
	/*Receives a primary sequence and topology. Build the angles section of
	 * protein
	 *
	 * 1) Have to start with 1, because residue number is started with 1.
	 * 2) r-1 means the index of sequence_primary starts in 0.
	 * 3) index_atm_angle is started with -1 because I is increased at
	 * set_residue_atoms_bond_from_topol_ff
	 */
	int index_atm_angle = -1;
	type_aminos_t amino_id;
	for (int r=1; r<=top->numres;r++){
        amino_id = primary_sequence[r-1].id;
		set_residue_atoms_angle_from_topol_ff(&index_atm_angle,&r,amino_id, top);
	}
}

static void build_section_dihedral(const amino_t *primary_sequence,top_global_t *top){
	/*Receives a primary sequence and topology. Build the dihedral angles
	 * section of protein.
	 *
	 * Here, dihedral angles means: torsional angles (phi and psi) and side chains
	 * angles
	 *
	 */
	int index_atm_angle = -1;
	int index_atm_angle_psi = -1;
	int index_atm_angle_omega = -1;
	int index_atm_angle_side_chains = -1;
	type_aminos_t amino_id;
	for (int r=1; r<=top->numres;r++){
        amino_id = primary_sequence[r-1].id;
        if (r == 37){
        	printf("%i\n", r);
        }
        if (amino_id != aX){
	        //Phi angle
	        set_residue_atoms_dihedral_phi_from_topol_ff(&index_atm_angle,&r,amino_id, top);
	        //Psi angle
	        set_residue_atoms_dihedral_psi_from_topol_ff(&index_atm_angle_psi,&r,amino_id, top);
	        //Omega angle
	        set_residue_atoms_dihedral_omega_from_topol_ff(&index_atm_angle_omega,&r,amino_id, top);        
	        // side chains angle
	        set_residue_atoms_dihedral_side_chains_from_topol_ff(&index_atm_angle_side_chains,&r,amino_id, top);
        }else{
        	if (r == 1){
        		//It means that ACE has 6 atoms. Therefore, It is N_terminal. So 6 first atoms
        		index_atm_angle = 5;
        	}else if (r == top->numres){
        		//It means that NME has 6 atoms. Therefore, It is C_terminal. So 6 last atoms
        		index_atm_angle = index_atm_angle + 5; 
        	}
        	
        }
	}
}

static void build_section_dihedral_angle_types(const amino_t *primary_sequence,
		top_global_t *top){
	/* This function build dihedral angle type section */
	int index_atm_angle = -1;
	type_aminos_t amino_id;
	for (int r=1; r<=top->numres;r++){
        amino_id = primary_sequence[r-1].id;
        set_residue_atoms_dihedral_angle_type_from_topol_ff(&index_atm_angle,&r,amino_id, top);
	}
}

void set_protein_charge_in_topol(top_global_t *top){
	/*Computes the charge of protein and store it in protein_charge record of
	 * topology.
	 * This charge can be used in genion program.
	 */
	float protein_charge = 0;
	for (int a = 0; a < top->numatom; a ++){
		protein_charge = protein_charge + top->top_global_atom[a].charge;
	}
	top->protein_charge = protein_charge;
}


void build_top_global(const amino_t *sequence_primary, top_global_t *top){
	/*Build the global topology*/
    build_sections_atom_and_residue_atoms(sequence_primary,top);    
    build_section_dihedral(sequence_primary,top);
    /* I commented because my focus is built the Z matrix which has all these
     * parameters. My time is shorter to review this part of program.
     * The Z matrix has all these parameters where I put manually.
    build_section_bonds(sequence_primary,top);
    build_section_angle_bond(sequence_primary,top);
     */
    //build_section_dihedral_angle_types(sequence_primary,top);
    set_protein_charge_in_topol(top);
}

void build_top_global_without_dihedral(const amino_t *sequence_primary,
		top_global_t *top){
	/*Build the global topology without dihedral section*/
    build_sections_atom_and_residue_atoms(sequence_primary,top);
    //build_section_dihedral(sequence_primary,top);
    /* I commented because my focus is built the Z matrix which has all these
     * parameters. My time is shorter to review this part of program.
     * The Z matrix has all these parameters where I put manually.
    build_section_bonds(sequence_primary,top);
    build_section_angle_bond(sequence_primary,top);
     */
    //build_section_dihedral_angle_types(sequence_primary,top);
    set_protein_charge_in_topol(top);
}
int get_num_atom_from_topol(const int *num_res, const type_atoms_t atm, const top_global_t *top){
	const top_global_atom_t* aux;
	aux = _get_top_global_atom_t_from_topol(num_res,atm,top);
	return aux->atom_number;
}

void set_atoms_from_topol_ff(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index){
	const topol_residues_t *residue_ff= get_residue_from_topol(amino_id);
	for (int a = 0; a<residue_ff->nr_atoms; a++){
		*last_index_atoms_top_glo = *last_index_atoms_top_glo +1;
		top->top_global_atom[*last_index_atoms_top_glo].atom_id = residue_ff->residue_atoms[a].atom_id;
		top->top_global_atom[*last_index_atoms_top_glo].charge = residue_ff->residue_atoms[a].charge;

		strcpy(top->top_global_atom[*last_index_atoms_top_glo].atom_name,residue_ff->residue_atoms[a].atom_name);
		strcpy(top->top_global_atom[*last_index_atoms_top_glo].res_name,residue_ff->res_name);
		top->top_global_atom[*last_index_atoms_top_glo].atom_number = *last_index_atoms_top_glo +1;//because index started in 0. atom started in 1.
		top->top_global_atom[*last_index_atoms_top_glo].res_number = *res_index +1;
		top->top_global_atom[*last_index_atoms_top_glo].amino_id = residue_ff->res_id;
	}
}

void set_atoms_from_topol_ff_N_Terminal(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index){
	const topol_residues_t *residue_ff= get_residue_from_topol_N_Terminal(amino_id);
	for (int a = 0; a<residue_ff->nr_atoms; a++){
		*last_index_atoms_top_glo = *last_index_atoms_top_glo +1;
		top->top_global_atom[*last_index_atoms_top_glo].atom_id = residue_ff->residue_atoms[a].atom_id;
		top->top_global_atom[*last_index_atoms_top_glo].charge = residue_ff->residue_atoms[a].charge;

		strcpy(top->top_global_atom[*last_index_atoms_top_glo].atom_name,residue_ff->residue_atoms[a].atom_name);
		strcpy(top->top_global_atom[*last_index_atoms_top_glo].res_name,residue_ff->res_name);
		top->top_global_atom[*last_index_atoms_top_glo].atom_number = *last_index_atoms_top_glo +1;//because index started in 0. atom started in 1.
		top->top_global_atom[*last_index_atoms_top_glo].res_number = *res_index +1;
		top->top_global_atom[*last_index_atoms_top_glo].amino_id = residue_ff->res_id;
	}
}

void set_atoms_from_topol_ff_C_Terminal(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index){
	const topol_residues_t *residue_ff= get_residue_from_topol_C_Terminal(amino_id);
	for (int a = 0; a<residue_ff->nr_atoms; a++){
		*last_index_atoms_top_glo = *last_index_atoms_top_glo +1;
		top->top_global_atom[*last_index_atoms_top_glo].atom_id = residue_ff->residue_atoms[a].atom_id;
		top->top_global_atom[*last_index_atoms_top_glo].charge = residue_ff->residue_atoms[a].charge;

		strcpy(top->top_global_atom[*last_index_atoms_top_glo].atom_name,residue_ff->residue_atoms[a].atom_name);
		strcpy(top->top_global_atom[*last_index_atoms_top_glo].res_name,residue_ff->res_name);
		top->top_global_atom[*last_index_atoms_top_glo].atom_number = *last_index_atoms_top_glo +1;//because index started in 0. atom started in 1.
		top->top_global_atom[*last_index_atoms_top_glo].res_number = *res_index +1;
		top->top_global_atom[*last_index_atoms_top_glo].amino_id = residue_ff->res_id;
	}
}

void save_topology(const char *path, const char *file_name, const top_global_t *top){
	_save_topology_file(path,file_name,top);
}

void create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein){
	_create_fasta_pdb(prot_name, chain_name, prot_seq, file_name_protein);
}

amino_t *load_protein2(const char *file_name_protein, int *n_residues,
		int *numatom, int *bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles, const input_parameters_t *in_para){
	amino_t *aux_amino;
	boolean_t has_his = bfalse;

	aux_amino = _load_fasta_pdb(file_name_protein, n_residues, numatom, bond_angles,
			num_protein_side_chains,num_dihedral_angles, &has_his);
	/* Commented for while. I have to research more details. I use -ignh option
	 * to ignorate the Hydrogen atoms for HIS
	if (has_his == btrue){
		replace_his(aux_amino, n_residues, numatom,bond_angles,
				num_protein_side_chains, num_dihedral_angles, in_para->path_database,
				in_para->path_local_execute, in_para->path_program_HIS_protonation,
				in_para->path_gromacs_programs);
	}
	*/
	return aux_amino;
}

/*Convert a protein structure (dihedral space) to array of pdb_atoms
 * structure (Cartesian space)
 */
void protein2pdbatoms(pdb_atom_t *pdb_atoms, const protein *prot,
		const top_global_t *top_global, const z_matrix_global_t *z_matrix) {
	_protein2pdbatoms(pdb_atoms,prot,top_global, z_matrix);
}

/*Convert an array of pdb_atoms structure (Cartesian space) to protein
 * structure (dihedral space)
 */
void pdbatoms2protein(protein *prot, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
   _pdbatoms2protein(prot, pdb_atoms,top_global);
}


static float get_bond_length_from_atom_parameters(type_atoms_t atom_id_1,
		type_atoms_t atom_id_2){
	return _get_bond_length_parameters(atom_id_1, atom_id_2);
}


static float get_angle_from_atom_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2, type_atoms_t atom_id_3){
	return _get_bond_angle_from_atom_parameters(atom_id_1, atom_id_2, atom_id_3);
}

static boolean_t is_special_atom(type_atoms_t atm){
	return _is_special_atom(atm);
}

static const special_atom_parameters_t * get_special_atom_parameters(type_atoms_t atm){
	return _get_special_atom_parameters(atm);
}

void set_residue_atom_values(type_atoms_t *aux_atom, int *res_id_aux,type_atoms_t atm,
		const int *res_id){
	/* This function has goal set values for residue and atom
	 * Because of special atoms of bond angles and bond length
	 * it is necessary check if atom (atm) is a special atom or not.
	 * */
	const special_atom_parameters_t *sp_atoms;
    if (is_special_atom(atm) == btrue){
    	sp_atoms = get_special_atom_parameters(atm);
    	*aux_atom = sp_atoms->atom_id;
    	*res_id_aux = *res_id + sp_atoms->residue_position;
    }else{
    	*aux_atom = atm;
    	*res_id_aux = *res_id;
    }

}


void build_protein_backbone(protein_backbone_t *protback,
		const top_global_t *top){
	/* Receives array of protein_backbone_t and the topology
	 * Store into array the index of atoms which compose the
	 * protein backbone (Carbon, Alfa Carbon, Nitrogen, HA and O or OT2)
	 * OT2 is necessary when residue is C terminal there is OT2 instead of
	 * Oxygen (O)
	 */
	type_atoms_t atm_aux;
	int res_id_aux;
	for (int i = 1; i <= top->numres;i++){
		protback[i-1].res_number = i;
		protback[i-1].atom_C = get_num_atom_from_topol(&i,atmC,top);
		protback[i-1].atom_Ca = get_num_atom_from_topol(&i,atmCA,top);
		protback[i-1].atom_N = get_num_atom_from_topol(&i,atmN,top);
		if (i == 1){//means that residue is N terminal.
			protback[i-1].atom_C_ = 0;
		}else{
			set_residue_atom_values(&atm_aux,&res_id_aux,atmC_,&i);
			protback[i-1].atom_C_ = get_num_atom_from_topol(&res_id_aux,atm_aux,top);
		}
		if (i == top->numres){//means that residue is C terminal.
			protback[i-1].atom_N_plus = 0;
		}else{
			set_residue_atom_values(&atm_aux,&res_id_aux,atmN_plus,&i);
			protback[i-1].atom_N_plus = get_num_atom_from_topol(&res_id_aux,atm_aux,top);
		}
	}
}

void show_protein_backbone(const protein_backbone_t *prot_back,
		const top_global_t *top){
	_show_protein_backbone(prot_back, top);
}

float get_bond_angle_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const int *atom_3, const top_global_t *top){
	return _get_bond_angle_from_topol(res,atom_1,atom_2,atom_3, top);
}

float get_bond_len_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const top_global_t *top){
	return _get_bond_len_from_topol(res,atom_1,atom_2, top);
}

type_aminos_t get_amino_id_from_res_id(const int *res_id,
		const top_global_t *top){
	return _get_amino_id_from_res_id(res_id,top);
}

boolean_t has_side_chain(const type_aminos_t *amino_id){
	return _has_side_chain(amino_id);
}

int number_side_chains(const type_aminos_t *amino_id){
	return _number_side_chains(amino_id);
}

type_aminos_t get_amino_id_3(char *c){
	return _get_amino_id_3(c);
}

void check_pdb_atoms_topology(const pdb_atom_t *pdb_atoms,
		const top_global_t *top){
	_check_pdb_atoms_topology(pdb_atoms,top);
}

int get_four_atom_number_for_bond_angle(
		type_dihedral_angles_t *tp_dihedral_angle,
		const int *res,
		const int *atom_1, const int *atom_2, const int *atom_3,
		const top_global_t *top){
    return _get_four_atom_number_and_type_for_dihedral_angle(tp_dihedral_angle,
    		res, atom_1, atom_2, atom_3, top);

}
int get_third_atom_number_for_bond_angle(const int *res,
		const int *atom_1, const int *atom_2, const top_global_t *top){
	return _get_third_atom_number_for_bond_angle(res, atom_1, atom_2, top);
}

int get_number_atom_bond(const int *res, const int *atom_1,
		const top_global_t *top){
	return _get_number_atom_bond(res, atom_1, top);
}

float get_charge_topol(const top_global_t *top){
	/* This function returns the value of charge of system
	 */
	return top->protein_charge;
}

int get_charge_topol_int(const top_global_t *top){
	/* This function returns the value of charge of system in integer
	 */
	return (int) top->protein_charge;
}

int get_number_atoms_from_res(const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top){
	return _get_number_atoms_from_res(res_id, amino_id, top);
}

int get_number_atoms_from_res_C_Terminal(const type_aminos_t *amino_id,
		const top_global_t *top){
	return _get_number_atoms_from_res_C_Terminal(amino_id, top);
}

const topol_residue_atoms_t* get_topol_residue_atoms_t_from_res(int *num_atom,
		const int *res_id, 	const type_aminos_t *amino_id,
		const top_global_t *top){
	return _get_topol_residue_atoms_t_from_res(num_atom, res_id, amino_id,top);
}

static void set_residue_atoms_dihedral_omega_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top){
	/*
	 * The value 0 for atom_number1, atom_number2,  atom_number3 and atom_number4 mean that there are not these atoms.
	 * In N_terminal case, the res_id -1 is impossible. Therefore, C_terminal case res_id + 1 is
	 * impossible too.
	 * */

	type_atoms_t aux_atom_1, aux_atom_2, aux_atom_3, aux_atom_4;
    int res_id_aux_1, res_id_aux_2, res_id_aux_3, res_id_aux_4;
    const topol_residues_t *residue_ff = NULL;

    if (*res_id == 1){ //N_terminal
    	residue_ff= get_residue_from_topol_N_Terminal(amino_id);
    }else if (*res_id <= top->numres -1){ //internal residues{
    	residue_ff= get_residue_from_topol(amino_id);
    }else if (*res_id == top->numres){ // C_terminal{
    	residue_ff= get_residue_from_topol_C_Terminal(amino_id);
    }

	*last_index_atm_angle = *last_index_atm_angle +1;
	top->top_global_dieh_omega[*last_index_atm_angle].res_number = *res_id;

	set_residue_atom_values(&aux_atom_1,&res_id_aux_1,residue_ff->residue_atoms_omega[0].atom_id_1,res_id);
	set_residue_atom_values(&aux_atom_2,&res_id_aux_2,residue_ff->residue_atoms_omega[0].atom_id_2,res_id);
	set_residue_atom_values(&aux_atom_3,&res_id_aux_3,residue_ff->residue_atoms_omega[0].atom_id_3,res_id);
	set_residue_atom_values(&aux_atom_4,&res_id_aux_4,residue_ff->residue_atoms_omega[0].atom_id_4,res_id);

	if (*res_id == 1){ //N_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}

	}else if (*res_id <= top->numres -1){ //internal residues
		top->top_global_dieh_omega[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                aux_atom_1,top);
		top->top_global_dieh_omega[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
				aux_atom_2,top);
		top->top_global_dieh_omega[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
				aux_atom_3,top);
		top->top_global_dieh_omega[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
				aux_atom_4,top);
	} else if (*res_id == top->numres){ // C_terminal
		if (res_id_aux_1 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number1 = 0;

		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number1 = get_num_atom_from_topol(&res_id_aux_1,
                   					aux_atom_1,top);
		}
		if (res_id_aux_2 == 0){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number2 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number2 = get_num_atom_from_topol(&res_id_aux_2,
								aux_atom_2,top);
		}
		if (res_id_aux_3 > top->numres){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number3 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number3 = get_num_atom_from_topol(&res_id_aux_3,
					aux_atom_3,top);
		}
		if (res_id_aux_4 > top->numres){
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number4 = 0;
		}else{
			top->top_global_dieh_omega[*last_index_atm_angle].atom_number4 = get_num_atom_from_topol(&res_id_aux_4,
					aux_atom_4,top);
		}
	}
}


void population2pdb_atom(pdb_atom_t **pop_pdb, protein **pop,
	const int *popSize, const top_global_t *top_global){
	/* This function converts a population of protein to population of pdb_atom_t */
	for (int i = 0; i < *popSize; i++){
       protein2pdbatoms(pop_pdb[i], pop[i], top_global, pop[i]->z_matrix);
	}
}
