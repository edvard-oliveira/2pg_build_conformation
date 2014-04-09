/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

/*HSE is one of HIS protonation conformation*/

topol_residue_atoms_t residue_atoms_HSE_N_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
		//{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.08},
		{atmCG,"CG","CPH1",0.22},
		{atmND1, "ND1",	"NR2", -0.70},
		{atmCE1,"CE1",	"CPH2",	0.25},
		{atmNE2, "NE2",	"NR2", -0.36},
		{atmCD2, "CD2",	"CPH1",	-0.05},
		{atmHB1, "HB1",	"HA", 0.09},
		{atmHB2, "HB2",	"HA", 0.09},
		//{atmHD1, "HD1",	"H", 0.32},
		{atmHE1, "HE1",	"HR1", 0.13},
		{atmHE2, "HE2",	"H", 0.32},
		//{atmH, "H",	"H", 0.32},
		{atmHD2, "HD2",	"HR3",	0.09}
};

topol_residue_atoms_t residue_atoms_HSE [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.08},
		{atmCG,"CG","CPH1",0.22},
		{atmND1, "ND1",	"NR2", -0.70},
		{atmCE1,"CE1",	"CPH2",	0.25},
		{atmNE2, "NE2",	"NR2", -0.36},
		{atmCD2, "CD2",	"CPH1",	-0.05},
		{atmHB1, "HB1",	"HA", 0.09},
		{atmHB2, "HB2",	"HA", 0.09},
		//{atmHD1, "HD1",	"H", 0.32},
		{atmHE1, "HE1",	"HR1", 0.13},
		{atmHE2, "HE2",	"H", 0.32},
		//{atmH, "H",	"H", 0.32},
		{atmHD2, "HD2",	"HR3",	0.09}
};


topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_HSE [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmND1} //chi2
};
