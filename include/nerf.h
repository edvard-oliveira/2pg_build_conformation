#ifndef OLD_NERF_H
#define OLD_NERF_H

#include "vector_math.h"

void NeRF(own_vector_t *D,const double *bond_len_BC, const double *bond_len_CD,
		const double *bond_angle_BCD, const double *torsion_BC,
		const own_vector_t *A, const own_vector_t *B, const own_vector_t *C);

#endif
