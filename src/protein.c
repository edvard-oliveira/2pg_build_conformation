#include <stdio.h>
#include <string.h>

#include "protein.h"
#include "defines.h"
#include "messages.h"

static void initialize_amino_values(amino_t *amino)
{
	amino->id = aNR;
	amino->late = NULL;
	amino->number_late = 0;
	amino->phi = 0;
	amino->psi = 0;
	amino->omega = PI;
	amino->pos_late = 0;
	strcpy(amino->aminoacido, "");
	strcpy(amino->idl3, "");
	amino->HP = -1;
	strcpy(amino->chainid, "A");
}

static void initialize_amino_vector(amino_t *amino, const int *nresiduos)
{
	for (int j=0; j < *nresiduos; j++) {
	        initialize_amino_values(&(amino[j]));
	}
}



static void initialize_amino(protein * protein_aux, const int *nresiduos){
	/*Initialize amino for protein */
	initialize_amino_vector(protein_aux->residuo, nresiduos);
}

static void set_general_information_protein(protein * protein_aux,
		const int *nresiduos, const int *nr_fitness){
	/*Set the general information for protein*/
	protein_aux->nr_fitness = *nr_fitness;
	protein_aux->nr_residues = *nresiduos;
	for (int f=0; f < *nr_fitness;f++){
		protein_aux->Fitness[f] = 0;
	}
}

protein * allocateProteinVector(const int *size, const int *nresiduos,
		const int *nr_fitness, const int *numatom){
	/*Allocates Protein as vector instead of matrix.
	 * It is used in order_pop_mono function, for example.
	 */
	protein *protein_aux;
	protein_aux = Malloc(protein,*size);
	for (int i = 0; i < *size; i++){
		protein_aux[i].residuo = allocateAmino(*nresiduos);
		protein_aux[i].Fitness = Malloc(long double, *nr_fitness);
		protein_aux->z_matrix = allocateZ_matrix(numatom);
		set_general_information_protein(&protein_aux[i],nresiduos,nr_fitness);		
		//initialize_amino(&protein_aux[i],nresiduos);
	}
    return protein_aux;
}

protein* allocateProtein(int nresiduos, int nr_fitness, int numatom){
	protein *protein_aux;
	protein_aux = Malloc(protein,1);
	protein_aux->Fitness = Malloc(long double, nr_fitness);
	protein_aux->residuo = allocateAmino(nresiduos);
	protein_aux->z_matrix = allocateZ_matrix(&numatom);
	set_general_information_protein(protein_aux,&nresiduos,&nr_fitness);
	//initialize_amino(protein_aux,&nresiduos);
    return protein_aux;
}

void copy_population_without_allocating(protein **pop_dest, protein **pop_source,
		const int *sizepop){
	/* This function copies the population pop_source to pop_dest
	 * For each individual is called copy_protein_values_without_allocating
	 * function because protein represents an individual of population and they
	 * have already allocated.
	 */
	const protein *prot_aux;
	for (int i = 0; i <*sizepop;i++){
		prot_aux = pop_source[i];
		copy_protein_values_without_allocating(pop_dest[i],prot_aux);
	}
}

void copy_population(protein **pop_dest, protein **pop_source,
		const int *sizepop){
	/* This function copies the population pop_source to pop_dest
	 * For each individual is called copy_protein function because protein
	 * represents an individual of population.
	 */
	const protein *prot_aux;
	for (int i = 0; i <*sizepop;i++){
		prot_aux = pop_source[i];
		copy_protein(pop_dest[i],prot_aux);
	}
}

void copy_population_specific_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop, const int *fit){
	/* This function copies the population pop_source to pop_dest. However,
	 * it copies one value of fitness only. This fitness is specific at fit
	 * For each individual is called copy_protein_specific_fitness function
	 * because protein represents an individual of population.
	 */
	const protein *prot_aux;
	for (int i = 0; i <*sizepop;i++){
		prot_aux = pop_source[i];
		copy_protein_specific_fitness(pop_dest[i],prot_aux, fit);
	}
}

void copy_population_without_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop){
	/* This function copies the population pop_source to pop_dest. However,
	 * it does not copy value of fitness. The fitness value will be obtained
	 * by other way. See build_fitness_population function in ea_meat.c, for
	 * example.
	 * For each individual is called opy_protein_without_fitness function
	 * because protein represents an individual of population.
	 */
	const protein *prot_aux;
	for (int i = 0; i <*sizepop;i++){
		prot_aux = pop_source[i];
		copy_protein_without_fitness(pop_dest[i],prot_aux);
	}
}

void copy_population_unchange_fitness(protein **pop_dest, protein **pop_source,
		const int *sizepop){
	/* This function copies the population pop_source to pop_dest
	 * For each individual is called copy_protein_unchange_fitness function
	 * because protein represents an individual of population.
	 */
	const protein *prot_aux;
	for (int i = 0; i <*sizepop;i++){
		prot_aux = pop_source[i];
		copy_protein_unchange_fitness(pop_dest[i],prot_aux);
	}
}

void swap_population(protein **pop1, protein **pop2){
	/* This function swaps pop2 with pop1
	 * Warning: It is possible because both populations has the same size.
	 * This is used at reproduce_population function.
	 */
	protein ** pop_aux = NULL;

	pop_aux = pop2;
	pop2 = pop1;
	pop1 = pop_aux;
}
void copy_protein(protein *p_dest, const protein *p_source){
	/* This function copies a protein from p_source to p_dest.
	 * In copy_amino is allocated side-chains. See
	 * copy_protein_values_without_allocating function that does not allocate
	 * anything.
	 */
	p_dest->nr_fitness = p_source->nr_fitness;
	p_dest->nr_residues = p_source->nr_residues;
	copy_z_matrix(p_dest->z_matrix, p_source->z_matrix);
	for (int f=0; f < p_dest->nr_fitness;f++){
		p_dest->Fitness[f] = p_source->Fitness[f];
	}
	for (int r = 0; r < p_dest->nr_residues; r++){
		copy_amino(&p_dest->residuo[r],&p_source->residuo[r]);
	}	
}

void copy_protein_values_without_allocating(protein *p_dest,
		const protein *p_source){
	/* This function copies a protein from p_source to p_dest.
	 * In this function is not allocated anything. See copy_protein that
	 * allocates side chains to p_dest because it uses copy_amino function.
	 */
	p_dest->nr_fitness = p_source->nr_fitness;
	p_dest->nr_residues = p_source->nr_residues;
	copy_z_matrix(p_dest->z_matrix, p_source->z_matrix);
	for (int f=0; f < p_dest->nr_fitness;f++){
		p_dest->Fitness[f] = p_source->Fitness[f];
	}
	for (int r = 0; r < p_dest->nr_residues; r++){
		copy_amino_values_without_allocating(&p_dest->residuo[r],
				&p_source->residuo[r]);
	}
}

void copy_protein_specific_fitness(protein *p_dest, const protein *p_source,
		const int *fit){
	/* This function copies a protein from p_source to p_dest. However,
	 * the fitness value copied what is fit value. Therefore the number of
	 * fitness must be 1
	 * In this function is not allocated anything. Therefore, both have to be
	 * allocated before.
	 *
	 */
	p_dest->nr_fitness = 1;
	p_dest->nr_residues = p_source->nr_residues;
	p_dest->Fitness[0] = p_source->Fitness[*fit];
	copy_z_matrix(p_dest->z_matrix, p_source->z_matrix);
	for (int r = 0; r < p_dest->nr_residues; r++){
		copy_amino(&p_dest->residuo[r],&p_source->residuo[r]);
	}
}

void copy_protein_without_fitness(protein *p_dest, const protein *p_source){
	/* This function copies a protein from p_source to p_dest. However,
	 * it does not copy value of fitness. The fitness value will be obtained
	 * by other way. See build_fitness_population function in ea_meat.c, for
	 * example.
	 * In this function is not allocated anything. Therefore, both have to be
	 * allocated before.
	 *
	 */
	p_dest->nr_fitness = 1;
	p_dest->nr_residues = p_source->nr_residues;
	p_dest->Fitness[0] = 0;
	copy_z_matrix(p_dest->z_matrix, p_source->z_matrix);
	for (int r = 0; r < p_dest->nr_residues; r++){
		copy_amino(&p_dest->residuo[r],&p_source->residuo[r]);
	}
}

void copy_protein_unchange_fitness(protein *p_dest, const protein *p_source){
	/* This function copies a protein from p_source to p_dest. However,
	 * it does change value of fitness. The fitness value will be obtained
	 * by other way.
	 * This function has the goal to copy two protein (individuals) but without
	 * changing fitness values. So, it is used when you want to copy values of
	 * amino acids (dihdral angles, for example)
	 * This function is called at select_individuals in ea_meat.c file.
	 * In this function is not allocated anything. Therefore, it is necessary to
	 * allocate before call it.
	 */
	p_dest->nr_residues = p_source->nr_residues;
	copy_z_matrix(p_dest->z_matrix, p_source->z_matrix);
	for (int r = 0; r < p_dest->nr_residues; r++){
		copy_amino(&p_dest->residuo[r],&p_source->residuo[r]);
	}
}

long double get_specific_fitnes(protein *p, const int *fit){
	/* Returns a specific fitness value
	 */
	return p->Fitness[*fit];
}

long double get_oposite_specific_fitnes(protein *p, const int *fit){
	/* Returns an opposite specific fitness value.
	 * There are fitness that must be maximized such as H_Bond and Hydrophilic.
	 * So, these fitness are stored in opposite value. Please see functions that
	 * calculate them such as compute_hbond_mult function.
	 */
	return p->Fitness[*fit]*(-1);
}

amino_t* allocateAmino(int nr_res){
	    amino_t* aux_amino;
	    aux_amino = Malloc(amino_t, nr_res);
	    initialize_amino_vector(aux_amino, &nr_res);
    	return aux_amino;
}

void deAllocateAmino(amino_t* aux_amino, int nr_res){
	//free(aux_amino);
	//printf("PRECISO IMPLEMENTAR deAllocateAmino \n");

	for (int i=0;i<nr_res;i++){
        if(aux_amino[i].late!=NULL){
            free(aux_amino[i].late);
        }
	}
	free(aux_amino);
}

void copy_amino(amino_t *amino_dest, const amino_t *amino_source){
	/* This function copies an amino from amino_source to amino_dest
	 *  The difference this function and copy_amino function is the last does
	 *  not allocate.
	 */
	amino_dest->HP = amino_source->HP;
	strcpy(amino_dest->aminoacido, amino_source->aminoacido);
	strcpy(amino_dest->chainid,amino_source->chainid);
	strcpy(amino_dest->idl3, amino_source->idl3);
	amino_dest->id = amino_source->id;
	amino_dest->number_late = amino_source->number_late;
	if (amino_dest->late != NULL){
		free(amino_dest->late);
	}
	amino_dest->late = Malloc(float,amino_dest->number_late);
	for (int i=0;i<amino_dest->number_late;i++){
		amino_dest->late[i] = amino_source->late[i];
	}
	amino_dest->pos_late = amino_source->pos_late;
	amino_dest->phi = amino_source->phi;
	amino_dest->psi = amino_source->psi;
	amino_dest->omega = amino_source->omega;
}

void copy_amino_values_without_allocating(amino_t *amino_dest,
		const amino_t *amino_source){
	/* This function copies an amino from amino_source to amino_dest
	 *  The difference this function and copy_amino function is the first does
	 *  not allocate.
	 */
	int i;
	amino_dest->HP = amino_source->HP;
	strcpy(amino_dest->aminoacido, amino_source->aminoacido);
	strcpy(amino_dest->chainid,amino_source->chainid);
	strcpy(amino_dest->idl3, amino_source->idl3);
	amino_dest->id = amino_source->id;
	amino_dest->number_late = amino_source->number_late;
	for (i=0;i<amino_dest->number_late;i++){
		amino_dest->late[i] = amino_source->late[i];
	}
	amino_dest->pos_late = amino_source->pos_late;
	amino_dest->phi = amino_source->phi;
	amino_dest->psi = amino_source->psi;
	amino_dest->omega = amino_source->omega;
}

void copy_amino_allocating(amino_t *amino_dest, const amino_t *amino_source){
	/* This function copies an amino from amino_source to amino_dest allocating
	 * if necessary.
	 *  The difference this function and copy_amino function is the last does
	 *  not allocate.
	 *  This function is used in build_initial_population_lib.c file
	 */
	amino_dest->HP = amino_source->HP;
	strcpy(amino_dest->aminoacido, amino_source->aminoacido);
	strcpy(amino_dest->chainid,amino_source->chainid);
	strcpy(amino_dest->idl3, amino_source->idl3);
	amino_dest->id = amino_source->id;
	amino_dest->number_late = amino_source->number_late;
	if (amino_source->late != NULL){
		if (amino_dest->late != NULL){
			fatal_error("amino_dest in copy_amino function can not be allocated before its late");
		}
		amino_dest->late = Malloc(float,amino_dest->number_late);
		for (int i=0;i<amino_dest->number_late;i++){
			amino_dest->late[i] = amino_source->late[i];
		}
	}
	amino_dest->pos_late = amino_source->pos_late;
	amino_dest->phi = amino_source->phi;
	amino_dest->psi = amino_source->psi;
	amino_dest->omega = amino_source->omega;
}

void copy_amino_allocating_only(amino_t *amino_dest, const amino_t *amino_source){
	/* This function copies an amino from amino_source to amino_dest allocating
	 * it.
	 * Important: THIS FUNCTION ONLY ALLOCATING. VALUES WILL BE OBTAINED IN
	 * OTHER PLACE.
	 *  The difference this function and copy_amino function is the last does
	 *  not allocate.
	 *  This function is used in build_initial_population_lib.c file
	 */
	amino_dest->HP = amino_source->HP;
	strcpy(amino_dest->aminoacido, amino_source->aminoacido);
	strcpy(amino_dest->chainid,amino_source->chainid);
	strcpy(amino_dest->idl3, amino_source->idl3);
	amino_dest->id = amino_source->id;
	amino_dest->number_late = amino_source->number_late;
	if (amino_source->late != NULL){
		if (amino_dest->late != NULL){
			fatal_error("amino_dest in copy_amino function can not be allocated before its late");
		}
		amino_dest->late = Malloc(float,amino_dest->number_late);
	}
	amino_dest->pos_late = amino_source->pos_late;
}

void deAllocateProtein(protein* prot){
	deAllocateAmino(prot->residuo, prot->nr_residues);
	if (prot->z_matrix != NULL){
		desAllocateZ_matrix(prot->z_matrix);
	}
    free(prot->Fitness);
    free(prot);
}

void deAllocateProteinArray(protein* prot,const int *size){
	free(prot);
	display_msg("Must be implemented - deAllocateProteinArray \n");
}

side_chains_t* allocate_side_chains(int num){
	side_chains_t* aux;
	aux = Malloc(side_chains_t, num);
	for (int i = 0; i < num;i++){
		aux[i].chi1 = 0;
		aux[i].chi2 = 0;
		aux[i].chi3 = 0;
		aux[i].chi4 = 0;
		aux[i].chi5 = 0;
	}
	return aux;
}


void desAllocate_side_chains(side_chains_t *aux){
        free(aux);
	printf("desAllocate_side_chains Must be implemented\n");
}

void set_amino_protein(protein *prot, const amino_t *amino_aux,
		const int *index_res){
	/* Receives a protein and an amino.
	 * This function copy amino's values to protein.
	 * See set_amino_protein_values_without_angles function there is difference.
	 */
	//Saving torsional angle
	prot->residuo[*index_res].id = amino_aux->id;
	strcpy(prot->residuo[*index_res].aminoacido, amino_aux->aminoacido);
	strcpy(prot->residuo[*index_res].idl3,  amino_aux->idl3);
	prot->residuo[*index_res].phi = amino_aux->phi;
	prot->residuo[*index_res].psi = amino_aux->psi;
	prot->residuo[*index_res].omega = amino_aux->omega;
	//Saving side chains angles
	prot->residuo[*index_res].number_late = amino_aux->number_late;
    if ( (prot->residuo[*index_res].late == NULL) &&
    		(amino_aux->number_late > 0)){
    	prot->residuo[*index_res].late = Malloc(float,amino_aux->number_late);
    }
    if (prot->residuo[*index_res].number_late > 0){
    	for (int s=0; s < amino_aux->number_late ;s++){
    		prot->residuo[*index_res].late[s] = amino_aux->late[s];
    	}
    }
}

void set_amino_protein_values_without_angles(protein *prot,
		const amino_t *amino_aux,const int *index_res){
	/* Receives a protein and amino.
	 * This function copies values from amino_aux to prot.
	 * The difference set_amino_protein and set_amino_protein_values_without_angles
	 * is the first copies all angles. The second does not copy angles values.
	 */
	prot->residuo[*index_res].id = amino_aux->id;
	strcpy(prot->residuo[*index_res].aminoacido, amino_aux->aminoacido);
	strcpy(prot->residuo[*index_res].idl3,  amino_aux->idl3);
	prot->residuo[*index_res].number_late = amino_aux->number_late;
	prot->residuo[*index_res].HP = amino_aux->HP;
	prot->residuo[*index_res].late = Malloc(float,amino_aux->number_late);
}

void show_population(protein **pop, const int *len_pop){
	char msg[50];
	const protein *prot_aux;
	for(int i = 0; i <*len_pop; i++){
		sprintf(msg,"Display values of individual %i \n",i+1 );
		display_msg(msg);
		prot_aux = pop[i];
		show_protein(prot_aux);
		sprintf(msg,"End of values of individual %i \n",i+1 );
		display_msg(msg);
	}
}

void show_protein(const protein *prot){
	/*This function receives a protein and shows it.
	 */
	char msg_aux[200];
	char msg[200];
	display_msg("Fitness value of individual \n");
    for (int f =0; f<prot->nr_fitness;f++){
    	sprintf(msg," %i %Lf \n",f, prot->Fitness[f] );
    	display_msg(msg);
    }
    for (int res = 0; res < prot->nr_residues;res++){
    	sprintf(msg," %s \n",prot->residuo[res].aminoacido );
    	display_msg(msg);
    	sprintf(msg_aux,"PHI: %f PSI: %f Omega %f ",prot->residuo[res].phi,
    			prot->residuo[res].psi,  prot->residuo[res].omega);
    	strcpy(msg,msg_aux);
		for (int rot = 0; rot < prot->residuo[res].number_late; rot++){
			sprintf(msg_aux,"Chi%i: %f ",rot, prot->residuo[res].late[rot]);
			strcat(msg, msg_aux);
		}

    	strcat(msg,"\n");
    	display_msg(msg);
    }
}

void show_fitness_protein(const protein *prot){
	/*Shows the fitness of protein (individual)*/
	char msg[100];
	for (int f = 0; f < prot->nr_fitness;f++){
		sprintf(msg,"Fitness %i is %Lf \n",f,prot->Fitness[f]);
		display_msg(msg);
	}
}

void show_fitness_population(protein **pop, const int *sizepop){
	char msg[100];
	const protein *prot_aux;
	for (int i = 0; i < *sizepop;i++){
		sprintf(msg,"Individual %i \n",i);
		prot_aux = pop[i];
		show_fitness_protein(prot_aux);
	}
}

void protein_population2array(protein *prot_v, protein **prot,
		const int *sizeProt){
	/* Converts protein matrix to protein vector
	 * prot_v represents a vector of protein
	 */
	protein *prot_aux;
	for (int i = 0; i < *sizeProt; i++){
		prot_aux = prot[i];
		copy_protein(&prot_v[i],prot_aux);
	}
}

void protein_array2population(protein **prot, protein *prot_v,
		const int *sizeProt){
	/* Converts protein vector to protein matrix
	 * prot_v represents a vector of protein
	 */
	protein *prot_aux;
	for (int i = 0; i < *sizeProt; i++){
		prot_aux = prot[i];
		copy_protein(prot_aux,&prot_v[i]);
	}
}

void force_torsion_angles_values(protein *prot){
	/* This function forces the values of torsion angles
	 * Torsion angles of backbone has a special values
	 * phi of N-Terminal and psi of C-Terminal must be zero.
	 * This function is used after genetic operations.
	 */
	prot->residuo[0].phi = 0;
	prot->residuo[prot->nr_residues-1].psi = 0;
}
