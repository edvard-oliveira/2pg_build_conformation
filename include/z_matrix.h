#include "z_matrix_types.h"
#include "topology_types.h"

z_matrix_global_t * allocateZ_matrix(const int *num_atoms);
void initialize_z_matrix_info(z_matrix_global_t * aux, const int *num_atoms);
void desAllocateZ_matrix(z_matrix_global_t * z_matrix);
void build_z_matrix(z_matrix_global_t *z_matrix, const top_global_t *top);
void build_z_matrix_DATABASE(z_matrix_global_t *z_matrix,
		const top_global_t *top, const char *path);
void show_z_matrix(const z_matrix_global_t *z_matrix);
void save_z_matrix_file(const char *path, const char *file_name,
			const z_matrix_global_t *z_matrix);
void copy_z_matrix_residue(z_matrix_global_t *dest, const z_matrix_global_t *source, 
	const int *first_atom_residue, const int *last_atom_residue);
void copy_z_matrix(z_matrix_global_t *dest, const z_matrix_global_t *source);
static void set_values_z_matriz_information(z_matrix_information_t *dest, 
	const z_matrix_information_t *source);

