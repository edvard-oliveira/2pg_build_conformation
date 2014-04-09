#include <stdlib.h>
#include <string.h>


#include "load_parameters.h"
#include "build_initial_population.h"
#include "parameters_type.h"
#include "messages.h"

int main(int argc, char *argv[]){
	input_parameters_t in_param;

	display_msg("Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	build_initial_population(&in_param);

	deAllocateload_parameters(&in_param);

	display_msg("Done Build Initial Population !!! \n");
	return 0;
}
