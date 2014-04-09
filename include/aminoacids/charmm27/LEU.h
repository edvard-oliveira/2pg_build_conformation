/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_LEU_N_Terminal [] = {
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
		{atmCG,"CG","CT1",-0.09},
		{atmCD1, "CD1",	"CT3" ,	-0.27},
		{atmCD2, "CD2",	"CT3" ,	-0.27},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHG,	"HG", "HA",	0.09},
		{atmHD11, "HD11", "HA",	0.09},
		{atmHD12, "HD12", "HA",	0.09},
		{atmHD13, "HD13", "HA",	0.09},
		{atmHD21, "HD21", "HA",	0.09},
		{atmHD22, "HD22", "HA",	0.09},
		{atmHD23, "HD23", "HA",	0.09}
};

topol_residue_atoms_t residue_atoms_LEU [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CT1",-0.09},
		{atmCD1, "CD1",	"CT3" ,	-0.27},
		{atmCD2, "CD2",	"CT3" ,	-0.27},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHG,	"HG", "HA",	0.09},
		{atmHD11, "HD11", "HA",	0.09},
		{atmHD12, "HD12", "HA",	0.09},
		{atmHD13, "HD13", "HA",	0.09},
		{atmHD21, "HD21", "HA",	0.09},
		{atmHD22, "HD22", "HA",	0.09},
		{atmHD23, "HD23", "HA",	0.09}
};

topol_residue_atoms_t residue_atoms_LEU_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CT1",-0.09},
		{atmCD1, "CD1",	"CT3" ,	-0.27},
		{atmCD2, "CD2",	"CT3" ,	-0.27},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHG,	"HG", "HA",	0.09},
		{atmHD11, "HD11", "HA",	0.09},
		{atmHD12, "HD12", "HA",	0.09},
		{atmHD13, "HD13", "HA",	0.09},
		{atmHD21, "HD21", "HA",	0.09},
		{atmHD22, "HD22", "HA",	0.09},
		{atmHD23, "HD23", "HA",	0.09}
};
topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_LEU [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmCD1}//chi2
};
