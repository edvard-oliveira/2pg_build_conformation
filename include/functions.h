#ifndef OLD_FUNCTIONS_H
#define OLD_FUNCTIONS_H

#include "protein.h"
#include "enums.h"
#include "topology_types.h"

type_fitness_energies_t str2type_fitness_energies(char *name_fitness_energies);
void type_fitness_energies2str(char *name_fitness_energies,
		const type_fitness_energies_t *type_fitness );
void readPopulationFile(protein **pop, const int *innumberpopulation,
		const char *path, const char *file_name_pop,
		const amino_t *primary_seq, const int *nresiduos, 
		const top_global_t *top_global);
void readPopulationFile_Diehdral(char *chBloco, FILE *arq_name_pop, protein **pop, 
	const int *innumberpopulation, const char *path, const char *file_name_pop,
		const amino_t *primary_seq, const int *nresiduos, const top_global_t *top_global);
void readPopulationFile_PDB(char *chBloco, FILE *arq_name_pop, protein **pop, 
	const int *innumberpopulation, const amino_t *primary_seq, 	
	const top_global_t *top_global);
int get_model_number(const char *s);
int search_amino(char chave, int max);
//amino_t *load_protein(char *file_name_protein, int *nresiduos, int *numatom);
protein** allocatePopulation(int inPopSize, int nresiduos, int nr_fitness, int numatom);
void deAllocatePopulation(protein** population, int inPopSize);
//void set_seqtyp(int seqtyp[], const amino_t *prot_prin, int nresiduo);
type_aminos_t get_amino_id(char c);
void build_pdb_file_name(char *pdb_file_name, const char *aux_name,
		const  char *__restrict prefix);
void set_population_file_name_with_pop_file_name(char *pop_file_name,
		const int *generation);
void set_numatom_from_topol(int *numatom, const int *nr_atm_ff);
void set_do_dssp_percentage(float *percentage_helix, float *percentage_beta,
		float *percentage_total, const char *path_local_execute,
		const char *file_name);

/*Set seqtyp array that is used for buildpdb*/


//void fromDiedralToCartezian(protein * prot);
/* Cast from Diedhral for cartezian space.
 * Atoms are numered from 1, but vector is
 * started from 0. Thus, we use n+1
 * rather than n.
 * We use n rather than n-1.
*/

#endif

