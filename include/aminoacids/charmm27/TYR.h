/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_TYR_N_Terminal [] = {
		{atmN, "N","NH1",-0.3},
		{atmCA, "CA","CT1",0.21},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
		//{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.1},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CA",0.0},
		{atmCD1, "CD1",	"CA", -0.115},
		{atmCD2, "CD2",	"CA", -0.115},
		{atmCE1, "CE1",	"CA", -0.115},
		{atmCE2, "CE2",	"CA", -0.115},
		{atmCZ,	"CZ", "CA",	0.11},
		{atmOH, "OH", "OH1", -0.54},
		{atmHB1, "HB1",	"HA",	0.09},
		{atmHB2, "HB2",	"HA",	0.09},
		{atmHD1, "HD1",	"HP",	0.115},
		{atmHD2, "HD2",	"HP", 0.115},
		{atmHE1, "HE1",	"HP",	0.115},
		{atmHE2, "HE2",	"HP", 0.115},
		{atmHH,	"HH", "H", 0.43}
};

topol_residue_atoms_t residue_atoms_TYR [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CA",0.0},
		{atmCD1, "CD1",	"CA", -0.115},
		{atmCD2, "CD2",	"CA", -0.115},
		{atmCE1, "CE1",	"CA", -0.115},
		{atmCE2, "CE2",	"CA", -0.115},
		{atmCZ,	"CZ", "CA",	0.11},
		{atmOH, "OH", "OH1", -0.54},
		{atmHB1, "HB1",	"HA",	0.09},
		{atmHB2, "HB2",	"HA",	0.09},
		{atmHD1, "HD1",	"HP",	0.115},
		{atmHD2, "HD2",	"HP", 0.115},
		{atmHE1, "HE1",	"HP",	0.115},
		{atmHE2, "HE2",	"HP", 0.115},
		{atmHH,	"HH", "H", 0.43}
};

topol_residue_atoms_t residue_atoms_TYR_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.34},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CA",0.0},
		{atmCD1, "CD1",	"CA", -0.115},
		{atmCD2, "CD2",	"CA", -0.115},
		{atmCE1, "CE1",	"CA", -0.115},
		{atmCE2, "CE2",	"CA", -0.115},
		{atmCZ,	"CZ", "CA",	0.11},
		{atmOH, "OH", "OH1", -0.54},
		{atmHB1, "HB1",	"HA",	0.09},
		{atmHB2, "HB2",	"HA",	0.09},
		{atmHD1, "HD1",	"HP",	0.115},
		{atmHD2, "HD2",	"HP", 0.115},
		{atmHE1, "HE1",	"HP",	0.115},
		{atmHE2, "HE2",	"HP", 0.115},
		{atmHH,	"HH", "H", 0.43}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_TYR [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmCD1} //chi2
};
