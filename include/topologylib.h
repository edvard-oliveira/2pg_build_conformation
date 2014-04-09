#include "topology_types.h"
#include "pdbatom.h"
#include "z_matrix_types.h"

const top_global_atom_t* _get_top_global_atom_t_from_topol(const int *num_res, const type_atoms_t atm, const top_global_t *top);
void check_number_atoms(const int *numatom, const top_global_t *top);
const topol_residues_t *get_residue_from_topol(type_aminos_t amino_id);
const topol_residues_t *get_residue_from_topol_N_Terminal(type_aminos_t amino_id);
const topol_residues_t *get_residue_from_topol_C_Terminal(type_aminos_t amino_id);
void set_amino_from_topol(amino_t *prot_prin, int index_seq_prin, const topol_residues_t top_res[]);
void set_numatom_from_topol(int *numatom, const int *nr_atm_ff);
void set_num_bond_angle_from_topol(int *num_bond_angles,
		const int * number);
void set_num_side_chains_from_topol(int *num_protein_side_chains,
		const int *number);
void set_num_dihedral_angles(int *num_dihedral_angles, const int *number);
void set_has_his(boolean_t *has_his, const type_aminos_t *aminoid);
static int _get_index_bond_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2);
static int _get_index_bond_angles_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2,type_atoms_t atom_id_3);
float _get_bond_length_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2);
float _get_bond_angle_from_atom_parameters(type_atoms_t atom_id_1, type_atoms_t atom_id_2, type_atoms_t atom_id_3);
static int _get_index_special_atom_parameters(type_atoms_t atm);
boolean_t _is_special_atom(type_atoms_t atm);
const special_atom_parameters_t * _get_special_atom_parameters(type_atoms_t atm);
float _get_bond_angle_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const int *atom_3, const top_global_t *top);
float _get_bond_len_from_topol(const int *res,const int *atom_1,
		const int *atom_2, const top_global_t *top);
void _protein2pdbatoms(pdb_atom_t *pdb_atoms, const protein *prot,
		const top_global_t *top_global, const z_matrix_global_t *z_matrix);
type_aminos_t _get_amino_id_from_res_id(const int *res_id,
		const top_global_t *top);
boolean_t _has_side_chain(const type_aminos_t *amino_id);
int _number_side_chains(const type_aminos_t *amino_id);
void _pdbatoms2protein(protein *prot, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global);
float _compute_diehdral_angle(const own_vector_t *a1, const own_vector_t *a2,
		const own_vector_t *a3, const own_vector_t *a4);
type_aminos_t _get_amino_id_3(char *c);
void _check_pdb_atoms_topology(const pdb_atom_t *pdb_atoms,
		const top_global_t *top);
int _get_four_atom_number_and_type_for_dihedral_angle(
		type_dihedral_angles_t *tp_dihedral_angle,
		const int *res,
		const int *atom_1, const int *atom_2, const int *atom_3,
		const top_global_t *top);
int _get_third_atom_number_for_bond_angle(const int *res,
		const int *atom_1, const int *atom_2, const top_global_t *top);
int _get_number_atom_bond(const int *res, const int *atom_1,
		const top_global_t *top);
void _type_of_diedhral_angle2str(char *str,
		const type_dihedral_angles_t *type_dihedral);
int _get_number_atoms_from_res(const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top);
int _get_number_atoms_from_res_C_Terminal(const type_aminos_t *amino_id,
		const top_global_t *top);
const topol_residues_t* _get_topol_residues_t_from_res(const int *res_id,
		const type_aminos_t *amino_id, const top_global_t *top);
const topol_residue_atoms_t* _get_topol_residue_atoms_t_from_res(int *num_atom,
		const int *res_id, 	const type_aminos_t *amino_id,
		const top_global_t *top);
static float _compute_phi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global);
static float _compute_psi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global);
static float _get_diedhral_angle_from_protein_based_on_z_matrix(
		const z_matrix_global_t *z_matrix,const protein *prot,
		const int *i_z, const int *res);
static void set_values_from_topology_2_pdbatoms(pdb_atom_t *pdb_atoms,
		const int *i_z,	const top_global_t *top_global,
		const z_matrix_global_t *z_matrix );
static int _get_atom_index_from_top_global_dihedral_side_chain_t(
		int atm_opt, const int *r, 	const int *chi,
		const top_global_t *top_global);
static float  _compute_side_chains_angles(const int *r, const int *chi,
		const pdb_atom_t *pdb_atoms, const top_global_t *top_global);

