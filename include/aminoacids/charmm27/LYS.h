/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_LYS_N_Terminal [] = {
		{atmN, "N","NH1",-0.3},
		{atmCA, "CA","CT1",0.21},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
		//{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG, "CG","CT2",-0.18},
		{atmCD, "CD","CT2",0.20},
		{atmCE, "CE", "CT2",0.21},
		{atmNZ,	"NZ", "NH3", -0.30},
		{atmHB1,"HB1","HA",	0.09},
		{atmHB2,"HB2","HA",	0.09},
		{atmHG1,"HG1","HA",	0.09},
		{atmHG2,"HG2","HA",	0.09},
		{atmHD1,"HD1","HA",	0.09},
		{atmHD2,"HD2","HA",	0.09},
		{atmHE1,"HE1","HA",	0.05},
		{atmHE2,"HE2","HA",	0.05},
		{atmHZ1,"HZ1","HA",	0.33},
		{atmHZ2,"HZ2","HA",	0.33},
		{atmHZ3,"HZ3","HA",	0.33}
};

topol_residue_atoms_t residue_atoms_LYS [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG, "CG","CT2",-0.18},
		{atmCD, "CD","CT2",0.20},
		{atmCE, "CE", "CT2",0.21},
		{atmNZ,	"NZ", "NH3", -0.30},
		{atmHB1,"HB1","HA",	0.09},
		{atmHB2,"HB2","HA",	0.09},
		{atmHG1,"HG1","HA",	0.09},
		{atmHG2,"HG2","HA",	0.09},
		{atmHD1,"HD1","HA",	0.09},
		{atmHD2,"HD2","HA",	0.09},
		{atmHE1,"HE1","HA",	0.05},
		{atmHE2,"HE2","HA",	0.05},
		{atmHZ1,"HZ1","HA",	0.33},
		{atmHZ2,"HZ2","HA",	0.33},
		{atmHZ3,"HZ3","HA",	0.33}
};

topol_residue_atoms_t residue_atoms_LYS_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG, "CG","CT2",-0.18},
		{atmCD, "CD","CT2",0.20},
		{atmCE, "CE", "CT2",0.21},
		{atmNZ,	"NZ", "NH3", -0.30},
		{atmHB1,"HB1","HA",	0.09},
		{atmHB2,"HB2","HA",	0.09},
		{atmHG1,"HG1","HA",	0.09},
		{atmHG2,"HG2","HA",	0.09},
		{atmHD1,"HD1","HA",	0.09},
		{atmHD2,"HD2","HA",	0.09},
		{atmHE1,"HE1","HA",	0.05},
		{atmHE2,"HE2","HA",	0.05},
		{atmHZ1,"HZ1","HA",	0.33},
		{atmHZ2,"HZ2","HA",	0.33},
		{atmHZ3,"HZ3","HA",	0.33}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_LYS [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmCD}, //chi2
		{atmCB, atmCG, atmCD, atmCE}, //chi3
		{atmCG, atmCD, atmCE, atmNZ}, //chi4
//		{atmCD, atmNE, atmCZ, atmNH1} //chi5
};
