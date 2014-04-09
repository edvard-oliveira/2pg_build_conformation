#include "parameters_type.h"
#include "topology.h"
#include "z_matrix.h"

void build_initial_population(const input_parameters_t *in_para);

static void save_pop_file(const char *path, const char *file_name,
		protein ** pop, int pop_size);
