#include<stdio.h>
#include<stdlib.h>

#include"futil.h"
#include"fitnessio.h"

static void write_header_generation(FILE *fit_file, const int *generation){
	fprintf (fit_file,";Generation %d \n", *generation);
	fprintf (fit_file,";Individual Fitness \n");

}

static void write_fitness_values(FILE *fit_file, const int *fit, protein ** pop,
		const int *pop_size){
	int i;
	for (i=0; i < *pop_size; i++){
		fprintf (fit_file,"%d %Lf \n", i, get_specific_fitnes(pop[i],fit));
	}
}

static void write_oposite_fitness_values(FILE *fit_file, const int *fit, protein ** pop,
		const int *pop_size){
	int i;
	for (i=0; i < *pop_size; i++){
		fprintf (fit_file,"%d %Lf \n", i, get_oposite_specific_fitnes(pop[i],fit));
	}
}

void _save_fitness_file(const char *path, const char *file_name,
		const int *fit, protein ** pop, const int *pop_size,
		const int *generation, const input_parameters_t *in_para){
        char *fname = path_join_file(path,file_name);
	FILE *fit_file = open_file(fname,fWRITE);
	write_header_generation(fit_file, generation);
	/* Fitness that must be maximized. When they are obtained these values
	 * are multiplied by -1 because ParadisEO works with these opposite values.
	 * However, when they will be stored, they must be written in original
	 * value */
	if ( (in_para->fitness_energies[*fit] == fit_hbond) ||
		(in_para->fitness_energies[*fit] == fit_hydrophilic)||
		(in_para->fitness_energies[*fit] == fit_hbond_main)	||
		(in_para->fitness_energies[*fit] == fit_stride_total) ||
		(in_para->fitness_energies[*fit] == fit_stride_helix)	||
		(in_para->fitness_energies[*fit] == fit_stride_beta) ) {
		write_oposite_fitness_values(fit_file,fit,pop,pop_size);
	}else{
		write_fitness_values(fit_file,fit,pop,pop_size);
	}
	fclose(fit_file);
	free(fname);
}
