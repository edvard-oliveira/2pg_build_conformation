#ifndef OLD_PROTEIN_H
#define OLD_PROTEIN_H

#include "enums.h"
#include "topology_types.h"
#include "z_matrix.h"

typedef struct sside_chains{
	float chi1;
	float chi2;
	float chi3;
	float chi4;
	float chi5;
}side_chains_t;


typedef struct samino {
   type_aminos_t id;
   int number_late;
   float phi;
   float psi;
   float omega;
   int pos_late;
   float *late;
   char aminoacido[2]; // code for 1 Letter
   char idl3[4]; // code for 3 Letters
   int HP;
   char chainid[2];
 }amino_t;

typedef struct sprotein {
   long double *Fitness; //MONO Fitness[0]
   amino_t *residuo;
   int nr_fitness;
   int nr_residues;
   z_matrix_global_t *z_matrix;
 }protein;

protein* allocateProtein(int nresiduos, int nr_fitness, int numatom);
protein * allocateProteinVector(const int *size, const int *nresiduos,
      const int *nr_fitness, const int *numatom);
void deAllocateProtein(protein* prot);
void deAllocateProteinArray(protein* prot,const int *size);
side_chains_t* allocate_side_chains(int num);
void desAllocate_side_chains();
void copy_population_without_allocating(protein **pop_dest, protein **pop_source,
		const int *sizepop);
void copy_population(protein **pop_dest, protein **pop_source,
		const int *sizepop);
void copy_population_specific_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop, const int *fit);
void copy_population_without_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop);
void copy_population_unchange_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop);
void swap_population(protein **pop1, protein **pop2);
void copy_protein(protein *p_dest, const protein *p_source);
void copy_protein_values_without_allocating(protein *p_dest,
		const protein *p_source);
void copy_protein_specific_fitness(protein *p_dest, const protein *p_source,
		const int *fit);
void copy_protein_without_fitness(protein *p_dest, const protein *p_source);
amino_t* allocateAmino(int nr_res);
void copy_protein_unchange_fitness(protein *p_dest, const protein *p_source);
long double get_specific_fitnes(protein *p, const int *fit);
long double get_oposite_specific_fitnes(protein *p, const int *fit);
void deAllocateAmino(amino_t* aux_amino, int nr_res);
void copy_amino(amino_t *amino_dest, const amino_t *amino_source);
void copy_amino_values_without_allocating(amino_t *amino_dest,
		const amino_t *amino_source);
void copy_amino_allocating(amino_t *amino_dest, const amino_t *amino_source);
void copy_amino_allocating_only(amino_t *amino_dest, const amino_t *amino_source);
void set_amino_protein(protein *prot, const amino_t *amino_aux,
		const int *index_res);
void set_amino_protein_values_without_angles(protein *prot,
		const amino_t *amino_aux,const int *index_res);
void show_protein(const protein *prot);
void show_population(protein **pop, const int *len_pop);
void show_fitness_protein(const protein *prot);
void show_fitness_population(protein **pop, const int *sizepop);
void protein_array2population(protein **prot, protein *prot_v,
		const int *sizeProt);
void force_torsion_angles_values(protein *prot);
void protein_population2array(protein *prot_v, protein **prot,
		const int *sizeProt);

static void initialize_amino(protein * protein_aux, const int *nresiduos);
static void initialize_amino_values(amino_t *amino);
static void initialize_amino_vector(amino_t *amino, const int *nresiduos);
static void set_general_information_protein(protein * protein_aux,
		const int *nresiduos, const int *nr_fitness);

#endif
