/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_ASN_N_Terminal [] = {
		{atmN,"N", "NH1",-0.47},
		{atmCA,"CA", "CT1",	0.07},
		{atmC,"C", "C",	0.51},
		{atmO,"O", "O",	-0.51},
		{atmH1, "H1","H",0.31},
		{atmH2, "H2","H",0.31},
		{atmH3, "H3","H",0.31},
		//{atmHN,"HN","H", 0.31},
		{atmHA,"HA", "HB", 0.09},
		{atmCB,"CB", "CT2", -0.18},
		{atmCG,"CG", "CC", 0.55},
		{atmOD1,"OD1", "O", -0.55},
		{atmND2,"ND2", "NH2", -0.62},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHD21,"HD21", "H", 0.32},
		{atmHD22,"HD22", "H", 0.30}
};

topol_residue_atoms_t residue_atoms_ASN [] = {
		{atmN,"N", "NH1",-0.47},
		{atmCA,"CA", "CT1",	0.07},
		{atmC,"C", "C",	0.51},
		{atmO,"O", "O",	-0.51},
		{atmHN,"HN","H", 0.31},
		{atmHA,"HA", "HB", 0.09},
		{atmCB,"CB", "CT2", -0.18},
		{atmCG,"CG", "CC", 0.55},
		{atmOD1,"OD1", "O", -0.55},
		{atmND2,"ND2", "NH2", -0.62},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHD21,"HD21", "H", 0.32},
		{atmHD22,"HD22", "H", 0.30}
};

topol_residue_atoms_t residue_atoms_ASN_C_Terminal [] = {
		{atmN,"N", "NH1",-0.47},
		{atmCA,"CA", "CT1",	0.07},
		{atmC,"C", "C",	0.51},
		{atmOT1,"O", "O",	-0.51},
		{atmOT2,"O", "O",	-0.51},
		//{atmO,"O", "O",	-0.51},
		{atmHN,"HN","H", 0.31},
		{atmHA,"HA", "HB", 0.09},
		{atmCB,"CB", "CT2", -0.18},
		{atmCG,"CG", "CC", 0.55},
		{atmOD1,"OD1", "O", -0.55},
		{atmND2,"ND2", "NH2", -0.62},
		{atmHB1,"HB1", "HA", 0.09},
		{atmHB2,"HB2", "HA", 0.09},
		{atmHD21,"HD21", "H", 0.32},
		{atmHD22,"HD22", "H", 0.30}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_ASN [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmOD1} //chi2
};

/*Represents the dihedral angle and its type of residue.*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ASN [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmN_plus, atmN, angl_psi},
	{atmN, atmCA, atmC, atmCB, angl_typ_dieh_trans_123_},
	{atmN, atmCA, atmCB, atmCG, angl_chi1 },
	{atmCA, atmCB, atmCG, atmOD1, angl_chi2 },
};
