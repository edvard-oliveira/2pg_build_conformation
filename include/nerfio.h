
#include"vector_types.h"

void save_nerf_file(const char *path, const char *file_name,
		const int *atom_index, const int *atm_num_topol_A, const char *atm_name_A,
		__const char *atm_name_B,
		__const char *atm_name_C, __const char *atm_name_D,
		const own_vector_t *D, const double *bond_len_BC,
		const double *bond_len_CD, const double *bond_angle_BCD,
		const double *torsion_BC, const own_vector_t *A,
		const own_vector_t *B, const own_vector_t *C);

static void write_header_nerf(FILE *nerf_file);
static void write_line_nerf(FILE *nerf_file, const int *atom_index,
		const int *atm_num_topol_A,
		const char *atm_name_A, __const char *atm_name_B, __const char *atm_name_C,
		__const char *atm_name_D, 	const own_vector_t *D,
		const double *bond_len_BC,	const double *bond_len_CD,
		const double *bond_angle_BCD, 	const double *torsion_BC,
		const own_vector_t *A, 	const own_vector_t *B, const own_vector_t *C);
