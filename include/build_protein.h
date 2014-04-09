#include "protein.h"
#include "build_protein_database_types.h"
#include "parameters_type.h"
#include "z_matrix_types.h"
#include "topology_types.h"

void initialize_build_protein(library_dihedral_info_t* lib_dihe_info,
		const amino_t *primary_sequence, const int *nrresidues,
		const char *path_database,int *nr_kind_res);
void build_random_protein_rotamer_library(protein ** pop,const int *index_pop,
		const amino_t *primary_sequence,const int *nr_kind_res,
		const library_dihedral_info_t* lib_dihe_info,
		const int *num_res);
void build_random_protein_gsl(protein ** pop,const int *index_pop,
		const amino_t *primary_sequence,const int *nr_kind_res,
		const library_dihedral_info_t* lib_dihe_info,
		const int *num_res);
void build_random_amino_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
void build_random_amino_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res );

static void load_database_diehdral_angles(library_dihedral_info_t* lib_dihe_info,
		const amino_t *primary_sequence, const int *nrresidues,
		const char *path_database,int *nr_kind_res);


