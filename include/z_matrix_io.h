#ifndef OLD_Z_MATRIX_IO_H
#define OLD_Z_MATRIX_IO_H

#include "z_matrix_types.h"

void _show_z_matrix(const z_matrix_global_t *z_matrix);
void _save_z_matrix_file(const char *path, const char *file_name,
		const z_matrix_global_t *z_matrix);
void read_file_distance_parameters_z_matrix(z_matrix_distance_parameters_t *dist_para,
		const char *path, const int *res_num);
void read_file_angle_parameters_z_matrix(z_matrix_angle_parameters_t *ange_para,
		const char *path, const int *res_num);
//Static functions
static void get_type_of_diedhral(const z_matrix_global_t *z_matrix,
		const int *i);
static void write_z_matrix_file(FILE *z_matrix_file,
		const z_matrix_global_t *z_matrix);
static void build_line_output_z_matrix(char *line, const int *i,
		const z_matrix_global_t *z_matrix);
static void build_header_z_matrix(char *line);
#endif
