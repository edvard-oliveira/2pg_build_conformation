/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_SER_N_Terminal [] = {
		{atmN, "N","NH1",-0.3},
		{atmCA, "CA","CT1",0.21},
		{atmC, "C","C",0.51},
		{atmO, "O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
		//{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.1},
		{atmCB, "CB","CT3",0.05},
		{atmOG,	"OG", "OH1",-0.66},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHG1, "HG1",	"H", 0.43}
};


topol_residue_atoms_t residue_atoms_SER [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC, "C","C",0.51},
		{atmO, "O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT3",0.05},
		{atmOG,	"OG", "OH1",-0.66},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHG1, "HG1",	"H", 0.43}
};

topol_residue_atoms_t residue_atoms_SER_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC, "C","C",0.51},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT3",0.05},
		{atmOG,	"OG", "OH1",-0.66},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHG1, "HG1",	"H", 0.43}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_SER [] = {
		{atmN,atmCA,atmCB,atmOG}, //chi1
};

