#ifndef OLD_BUILD_PROTEIN_DATABASE_TYPES
#define OLD_BUILD_PROTEIN_DATABASE_TYPES

#include "enums.h"
#include "protein.h"


/* Informs if database has already initialized or not.
 * This is used in load_database_diehdral_angles function, for example
 */
typedef struct sdatabase_initialized{
	boolean_t initialized;
}database_initialized_t;

/*Represents parameters to load
 * all possible dihedral angles
 * Load from dihedral angles library
 */
typedef struct samino_database_parameters{
	type_aminos_t aminoid;
	long int max_torsional_angles;
	const char *file_torsion;
	long int max_side_chains_angles;
	int num_side_angles;
	const char *file_side_chains;
	const char *format_file;
}amino_database_parameters_t;

/*Represents the information of
 * torsional dihedral angles
 * from database.
 */
typedef struct slibrary_dihedral_info_tors{
	int id;
	float phi;
	float psi;
}library_dihedral_info_tors_t;

/*Represents the information of
 * side chains dihedral angles
 * from database.
 */
typedef struct slibrary_dihedral_info_side_chains{
	int id;
	side_chains_t* side_chains_angles;
	float freq;
	float despad;
}library_dihedral_info_side_chains_t;


/*Represents the informations of
 * diehdral angles (torsional and side chains)
 * from database
 */
typedef struct slibrary_dihedral_info{
	type_aminos_t aminoid;
	int num_torsional;
	library_dihedral_info_tors_t *torsional;
	int num_side_chains;
	int num_angles_angles;
	library_dihedral_info_side_chains_t *side_chains;
}library_dihedral_info_t;


#endif
