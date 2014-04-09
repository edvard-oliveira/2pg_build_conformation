/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_ILE_N_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.21},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
//		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.1},
		{atmCB, "CB","CT2",-0.18},
		{atmCG1,"CG1","CT2",-0.18},
		{atmCG2, "CG2",	"CT3", -0.27},
		{atmCD, "CD",	"CT3" ,	-0.27},
		{atmHB,  "HB", "HA",	0.09},
		{atmHG11, "HG11", "HA",	0.09},
		{atmHG12, "HG12", "HA",	0.09},
		{atmHG21, "HG21", "HA",	0.09},
		{atmHG22, "HG22", "HA",	0.09},
		{atmHG23, "HG23", "HA",	0.09},
		{atmHD1, "HD1", "HA",	0.09},
		{atmHD2, "HD2", "HA",	0.09},
		{atmHD3, "HD3", "HA",	0.09}
};

topol_residue_atoms_t residue_atoms_ILE [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG1,"CG1","CT2",-0.18},
		{atmCG2, "CG2",	"CT3", -0.27},
		{atmCD, "CD",	"CT3" ,	-0.27},
		{atmHB,  "HB", "HA",	0.09},
		{atmHG11, "HG11", "HA",	0.09},
		{atmHG12, "HG12", "HA",	0.09},
		{atmHG21, "HG21", "HA",	0.09},
		{atmHG22, "HG22", "HA",	0.09},
		{atmHG23, "HG23", "HA",	0.09},
		{atmHD1, "HD1", "HA",	0.09},
		{atmHD2, "HD2", "HA",	0.09},
		{atmHD3, "HD3", "HA",	0.09}
};


topol_residue_atoms_t residue_atoms_ILE_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.34},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG1,"CG1","CT2",-0.18},
		{atmCG2, "CG2",	"CT3", -0.27},
		{atmCD, "CD",	"CT3" ,	-0.27},
		{atmHB,  "HB", "HA",	0.09},
		{atmHG11, "HG11", "HA",	0.09},
		{atmHG12, "HG12", "HA",	0.09},
		{atmHG21, "HG21", "HA",	0.09},
		{atmHG22, "HG22", "HA",	0.09},
		{atmHG23, "HG23", "HA",	0.09},
		{atmHD1, "HD1", "HA",	0.09},
		{atmHD2, "HD2", "HA",	0.09},
		{atmHD3, "HD3", "HA",	0.09}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_ILE [] = {
		{atmN, atmCA, atmCB, atmCG1}, //chi1
		{atmCA, atmCB, atmCG1, atmCD}, //chi2
};
