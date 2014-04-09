#ifndef OLD_BUILD_PROTEIN_LIB
#define OLD_BUILD_PROTEIN_LIB

#include <stdio.h>

#include "enums.h"
#include "build_protein_database_types.h"
#include "protein.h"

void set_database_started(boolean_t b);
database_initialized_t * get_database_started();

void _load_amino_database(library_dihedral_info_t *lib_dihe, type_aminos_t amino,
		const int *index_res, const char *database);
library_dihedral_info_t* allocate_library_dihedral_info(int num_res);
void desAllocate_library_dihedral_info(library_dihedral_info_t * lib_dihe_info, const int *nr_res);
library_dihedral_info_tors_t * allocate_library_dihedral_info_tors(int num_angles);
library_dihedral_info_side_chains_t * allocate_library_dihedral_info_side_chains(int num_angles);
void _build_random_amino_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
void _build_random_amino_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
int _get_index_library_dihedral_info(const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
boolean_t _search_amino(const amino_t *aminos, const type_aminos_t* amino,
		const int *len_aminos);
amino_t * _get_unique_res(int *nr_kind_res,const amino_t *primary_sequence,
		const int *nrresidues);
float _get_side_chain_angle_from_lib(const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino);

static void desAllocate_library_dihedral_info_tors (library_dihedral_info_t * lib_dihe_info, const int *nr_res);
static void desAllocate_library_dihedral_info_side_chains();
static void _check_max_number(const long int *max_param,const long int *max_file);
static void check_amino_type(type_aminos_t *amino);
static int get_index_amino_database_parameters(type_aminos_t amino);
static void _get_line_values(FILE *side_chain_file, const type_aminos_t *amino, const char *format,
		float *chi1, float *chi2, float *chi3, float *chi4, float *chi5,
		float *freq, float *despad, int *fscanfError);
static void _side_chain_database_line2lib_dieh_info(library_dihedral_info_t *lib_dihe, const int *index_res,
		const type_aminos_t* amino,const float *chi1,const float *chi2,const float *chi3,const float *chi4,
		const float *chi5, const float *freq,const float *despad, const int *index_angles);
static void _build_random_dihedral_angles_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
static void _check_side_chain_angle_library(const int *index_amino,
		const library_dihedral_info_t* lib_dihe_info);
static void _build_random_dihedral_angles_gsl(amino_t *amino_aux,
		const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, const int *nr_kind_res);
static double _get_random_angle_gsl();
static void _set_random_side_chains_rotamer_library(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino);
static void _set_random_side_chains_gsl(amino_t *amino_aux, const type_aminos_t *amino,
		const library_dihedral_info_t* lib_dihe_info, int *index_lib,
		const int *index_amino);

#endif
