#ifndef OLD_FITNESSIO_H
#define OLD_FITNESSIO_H

#include"protein.h"
#include"parameters_type.h"

void _save_fitness_file(const char *path, const char *file_name,
		const int *fit, protein ** pop, const int *pop_size,
		const int *generation, const input_parameters_t *in_para);
#endif
