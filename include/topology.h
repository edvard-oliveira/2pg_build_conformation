#ifndef OLD_TOPOLOGY_H
#define OLD_TOPOLOGY_H

#include "protein.h"
#include "pdbatom.h"
#include "topology_types.h"
#include "z_matrix_types.h"
#include "parameters_type.h"

void protein2pdbatoms(pdb_atom_t *pdb_atoms, const protein *prot,
		const top_global_t *top_global, const z_matrix_global_t *z_matrix);
/*Convert a protein structure (dihedral space) to array of pdb_atoms structure (Cartesian space) */
void pdbatoms2protein(protein *prot, const pdb_atom_t *pdb_atoms, const top_global_t *top_global);
/*Convert an array of pdb_atoms structure (Cartesian space) to protein structure (dihedral space)*/

void build_top_global(const amino_t *sequence_primary, top_global_t *top);
void build_top_global_without_dihedral(const amino_t *sequence_primary,
		top_global_t *top);
top_global_t *allocateTop_Global(const int *numatom, const int *numres,
		const int *num_bond_angles, const int *num_side_chains,
		const int *number_dihedrals_type);
void  desAllocateTop_Global(top_global_t *top_global);
protein_backbone_t *allocateProtein_backbone(const int *numres);
void  desAllocateProtein_backbone(protein_backbone_t *protein_backbone);
amino_t *load_protein2(const char *file_name_protein, int *n_residues,
		int *numatom, int *bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles, const input_parameters_t *in_para);
void set_atoms_from_topol_ff(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index);
void set_atoms_from_topol_ff_N_Terminal(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index);
void set_atoms_from_topol_ff_C_Terminal(top_global_t *top,
		int *last_index_atoms_top_glo,type_aminos_t amino_id,
		const int *res_index);
void save_topology(const char *path, const char *file_name, const top_global_t *top);
void create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein);
void check_number_atoms(const int *numatom, const top_global_t *top);
int get_num_atom_from_topol(const int *num_res, const type_atoms_t atm, const top_global_t *top);
void build_protein_backbone(protein_backbone_t *protback, const top_global_t *top);
void show_protein_backbone(const protein_backbone_t *prot_back,
		const top_global_t *top);
float get_bond_angle_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const int *atom_3, const top_global_t *top);
float get_bond_len_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const top_global_t *top);
type_aminos_t get_amino_id_from_res_id(const int *res_id,
		const top_global_t *top);
boolean_t has_side_chain(const type_aminos_t *amino_id);
int number_side_chains(const type_aminos_t *amino_id);
type_aminos_t get_amino_id_3(char *c);
void check_pdb_atoms_topology(const pdb_atom_t *pdb_atoms,
		const top_global_t *top);
int get_four_atom_number_for_bond_angle(
		type_dihedral_angles_t *tp_dihedral_angle,
		const int *res,
		const int *atom_1, const int *atom_2, const int *atom_3,
		const top_global_t *top);
int get_third_atom_number_for_bond_angle(const int *res,
		const int *atom_1, const int *atom_2, const top_global_t *top);
int get_number_atom_bond(const int *res, const int *atom_1,
		const top_global_t *top);
void set_residue_atom_values(type_atoms_t *aux_atom, int *res_id_aux,type_atoms_t atm, const int *res_id);
float get_charge_topol(const top_global_t *top);
int get_charge_topol_int(const top_global_t *top);
int get_number_atoms_from_res(const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top);
int get_number_atoms_from_res_C_Terminal(const type_aminos_t *amino_id,
		const top_global_t *top);
const topol_residue_atoms_t* get_topol_residue_atoms_t_from_res(int *num_atom,
		const int *res_id, 	const type_aminos_t *amino_id,
		const top_global_t *top);
void population2pdb_atom(pdb_atom_t **pop_pdb, protein **pop,
	const int *popSize, const top_global_t *top_global);

/*static functions*/
static float get_bond_length_from_atom_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2);
static float get_angle_from_atom_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2, type_atoms_t atom_id_3);
static void set_residue_atoms_bond_from_topol_ff(int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_bond_from_topol_ff_N_Terminal(
		int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_bond_from_topol_ff_C_Terminal(
		int *last_index_atm_bond,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
/*static void set_residue_atoms_bond_from_topol_ff_N_Terminal(
		int *last_index_atm_bond, const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
	*/
static void set_residue_atoms_angle_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void build_sections_atom_and_residue_atoms(const amino_t *sequence_primary, top_global_t *top);
static void build_section_bonds(top_global_t *top);
static void build_section_angle_bond(const amino_t *primary_sequence,top_global_t *top);
static void build_section_dihedral(const amino_t *primary_sequence,top_global_t *top);
static boolean_t is_special_atom(type_atoms_t atm);
static const special_atom_parameters_t * get_special_atom_parameters(type_atoms_t atm);
static void set_residue_atoms_dihedral_angle_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_dihedral_angle_type_from_topol_ff(int *last_index_atm_angle,
		const int *res_id, type_aminos_t amino_id,	const top_global_t *top);
static void set_residue_atoms_dihedral_phi_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_dihedral_psi_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_dihedral_side_chains_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);
static void set_residue_atoms_dihedral_omega_from_topol_ff(int *last_index_atm_angle,const int *res_id, type_aminos_t amino_id,
		const top_global_t *top);

#endif
