#ifndef OLD_NERF_TOPOL_TYPES_H
#define OLD_NERF_TOPOL_TYPES_H

#include "enums.h"
#include "topology_types.h"

typedef struct sz_matrix_information{
	type_atoms_t atomid;
	int index_top;//Index of reference atom at global topology
	char atomname[5];//name of reference atom
	int atom_reference;//number of atom reference atom
	int atom_connected;//number of atom connected with atom reference
	double bond_len;
	double bond_len_2; //it is used in nerf algorithm
	int atom_angle; //number of atom which build bond angle with atom_connected and atom_reference
	double bond_angle;
	int dihedral_connect;// number of atom which build a dihedral angle
	float dihedral_angle;
	type_dihedral_angles_t tpAngle;
}z_matrix_information_t;

typedef struct sz_matrix_global{
	int num_elements; // Number of elements
	z_matrix_information_t *z_matrix_info;
	//const protein_backbone_t *prot_back;
}z_matrix_global_t;

typedef struct sz_matrix_distance_parameters{
	int res_num;
	double N_CA;
	double CA_C;
	double C_Nplus;
}z_matrix_distance_parameters_t;


typedef struct sz_matrix_angle_parameters{
	int res_num;
	double N_CA_C;
	double CA_C_Nplus;
	double C_Nplus_CAplus;
}z_matrix_angle_parameters_t;



#endif
