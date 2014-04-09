/* This file is the implementation of NeRF algorithm. This algorithm
 * converts from Torsion Space to Cartesian Space
 *
 * This implementation was based on:
 * authores: JEROD PARSONS, J. BRADLEY HOLMES, J. MAURICE ROJAS, JERRY TSAI,
 * CHARLIE E. M. STRAUSS
 * Title: Practical Conversion from Torsion Space to Cartesian Space for
 * In Silico Protein Synthesis
 * Published: online in Wiley InterScience (www.interscience.wiley.com)
 * DOI: 10.1002/jcc.20237
 * I has implemented NeRF algorithm.
 * */

#include<math.h>

#include"nerf.h"

void NeRF(own_vector_t *D,const double *bond_len_BC, const double *bond_len_CD,
		const double *bond_angle_BCD, const double *torsion_BC,
		const own_vector_t *A, const own_vector_t *B, const own_vector_t *C){
 /*
 * 1) The letters M, n, BC were based on paper above.
 * 2) All angles must be radius.
 */
	own_vector_long_t *v_res;
	own_vector_long_t vD2;
	own_vector_long_t vBC;
	own_vector_long_t vAB;
	own_vector_long_t vN;
	own_vector_long_t vM;
	long double mod_vet_n;

	//Computes bc
	sub_vector_long(&vBC,B,C);
	vBC.x = vBC.x/(*bond_len_BC);
	vBC.y = vBC.y/(*bond_len_BC);
	vBC.z = vBC.z/(*bond_len_BC);

	//Computes n
    sub_vector_long(&vAB,A,B);
    cross_product_long(&vN,&vAB,&vBC);
    mod_vet_n = mod_vector_long(&vN);
    vN.x = vN.x/mod_vet_n;
    vN.y = vN.y/mod_vet_n;
    vN.z = vN.z/mod_vet_n;

    //Compute M
    cross_product_long(&vM,&vN,&vBC);

	//Computes D2
    vD2.x = *bond_len_CD*cos(*bond_angle_BCD);
    vD2.y = *bond_len_CD*cos(*torsion_BC)*sin(*bond_angle_BCD);
    vD2.z = *bond_len_CD*sin(*torsion_BC)*sin(*bond_angle_BCD);

    //Compute final position (x,y,z)
    D->x = vBC.x*vD2.x + vM.x*vD2.y + vN.x*vD2.z + C->x;
    D->y = vBC.y*vD2.x + vM.y*vD2.y + vN.y*vD2.z + C->y;
    D->z = vBC.z*vD2.x + vM.z*vD2.y + vN.z*vD2.z + C->z;
}
