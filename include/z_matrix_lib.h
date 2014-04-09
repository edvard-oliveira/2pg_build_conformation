#include "z_matrix_types.h"
#include "topology_types.h"
#include "enums.h"

void _build_z_matrix(z_matrix_global_t *z_matrix, const top_global_t *top);
void _build_z_matrix_database(z_matrix_global_t *z_matrix,
		const top_global_t *top, const char *path);
void _get_type_of_dihedral(char stype[], const z_matrix_global_t *z_matrix,
		const int *i);

static void _set_type_dihedral_angle(z_matrix_global_t *z_matrix,
		const int *index_z_matrix, const type_dihedral_angles_t *type_dihedral);

static int get_index_atom_reference(const int *i, const top_global_t *top);
static void _set_information_atom_z_matrix(z_matrix_global_t *z_matrix,
		int *index_atom_reference,	int *index_z_matrix,
		const type_atoms_t *atom_reference, const int *res_id,
		const top_global_t *top);
static void _add_atom_z_matriz(z_matrix_global_t *z_matrix, int *index_z_matrix,
		const int *res_id,
		type_atoms_t atom_reference, type_atoms_t atom_conected,
		type_atoms_t atom_angle,type_atoms_t atom_dihedral,
		type_dihedral_angles_t type_dihedral,
		float bond_len, float bond_angle, float bond_len_2,
		const top_global_t *top);
static void _add_atoms_backbone(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_backbone_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_backbone_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_backbone_N_Terminal_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_backbone_C_terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_backbone_C_terminal_DATABASE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
const top_global_t *top);
static void _add_atoms_side_chains_HIS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_VAL(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_TYR(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_TRP(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_THR(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_PRO(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_PHE(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_MET(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_SER(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_LYS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_LEU(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_ILE(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_GLU(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_GLN(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_CYS(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_ASP(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_ASN(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_ALA(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_atoms_side_chains_ARG(z_matrix_global_t *z_matrix,
int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top);
static void _add_hydrogen_atoms_backbone_N_Terminal_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_N_Terminal_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_VAL(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_TYR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_TRP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_THR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_PHE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_MET(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_SER(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_LYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_LEU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_ILE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_GLU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_GLN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_CYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_ASP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_ASN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_GLY(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_GLY_N_Terminal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_ARG(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_backbone_ALA(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const type_aminos_t *amino_id,
		const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_VAL(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_HSE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_HSD(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_HIS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_TYR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_TRP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_THR(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_PRO(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_PHE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_MET(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_SER(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_LYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_LEU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_ILE(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_GLU(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_GLN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_CYS_Internal(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_CYS(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_ASP(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_ASN(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_ARG(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
static void _add_hydrogen_atoms_side_chains_ALA(z_matrix_global_t *z_matrix,
		int *index_z_matrix, const int *res_id, const top_global_t *top);
