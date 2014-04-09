/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */


topol_residue_atoms_t residue_atoms_PRO_N_Terminal [] = {
		{atmN, "N","N",-0.07},
		{atmCA, "CA","CP1",0.16},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.24},
		{atmH2, "H2","H",0.24},
		//{atmH3, "H3","H",0.33},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CT1",-0.18},
		{atmCD,	"CD", "CP3", 0.16},
		{atmHB1, "HB1",	"HA", 0.09},
		{atmHB2, "HB2",	"HA", 0.09},
		{atmHD1, "HD1",	"HA", 0.09},
		{atmHD2, "HD2",	"HA", 0.09},
		{atmHG1, "HG1",	"HA", 0.09},
		{atmHG2, "HG2",	"HA", 0.09}
};

topol_residue_atoms_t residue_atoms_PRO [] = {
		{atmN, "N","N",-0.29},
		{atmCA, "CA","CP1",0.02},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CT1",-0.18},
		{atmCD,	"CD", "CP3", 0.00},
		{atmHB1, "HB1",	"HA", 0.09},
		{atmHB2, "HB2",	"HA", 0.09},
		{atmHD1, "HD1",	"HA", 0.09},
		{atmHD2, "HD2",	"HA", 0.09},
		{atmHG1, "HG1",	"HA", 0.09},
		{atmHG2, "HG2",	"HA", 0.09}
};

topol_residue_atoms_t residue_atoms_PRO_C_Terminal [] = {
		{atmN, "N","N",-0.29},
		{atmCA, "CA","CP1",0.02},
		{atmC,"C","C",0.34},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG,"CG","CT1",-0.18},
		{atmCD,	"CD", "CP3", 0.00},
		{atmHB1, "HB1",	"HA", 0.09},
		{atmHB2, "HB2",	"HA", 0.09},
		{atmHD1, "HD1",	"HA", 0.09},
		{atmHD2, "HD2",	"HA", 0.09},
		{atmHG1, "HG1",	"HA", 0.09},
		{atmHG2, "HG2",	"HA", 0.09}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_PRO [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmCD} //chi2
};
