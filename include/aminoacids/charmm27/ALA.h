/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_ALA_N_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.21},
		{atmC, "C","C",0.51},
		{atmO, "O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
//		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.1},
		{atmCB, "CB","CT3",-0.27},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHB3, "HB3","HA",0.09}
};

topol_residue_atoms_t residue_atoms_ALA [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC, "C","C",0.51},
		{atmO, "O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT3",-0.27},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHB3, "HB3","HA",0.09}
};


topol_residue_atoms_t residue_atoms_ALA_C_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC, "C","C",0.51},
		{atmOT1,"O","O",-0.67},
		{atmOT2,"O","O",-0.67},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT3",-0.27},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2, "HB2","HA",0.09},
		{atmHB3, "HB3","HA",0.09}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_ALA [] = {
		{atmC_,atmN},//in aminoacids.rtp was wrote N+
		//{atmN,atmHN},
		{atmN,atmCA},
		//{atmCA,atmHA},
		{atmCA,atmCB},
		//{atmCB,atmHB1},
		//{atmCB,atmHB2},
		//{atmCB,atmHB3},
		{atmCA,atmC},
		{atmC,atmO}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_ALA_N_Terminal [] = {
		//{atmN,atmH1},
		//{atmN,atmH2},
		//{atmN,atmH3},
		{atmN,atmCA},
		{atmCA,atmCB},
		//{atmCA,atmHA},
		//{atmCB,atmHB1},
		//{atmCB,atmHB2},
		//{atmCB,atmHB3},
		{atmCA,atmC},
		{atmC,atmO},
		{atmC,atmN_plus}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_ALA_C_Terminal [] = {
		{atmC_,atmN},//in aminoacids.rtp was wrote N+
		//{atmN,atmHN},
		{atmN,atmCA},
		//{atmCA,atmHA},
		{atmCA,atmCB},
		//{atmCB,atmHB1},
		//{atmCB,atmHB2},
		//{atmCB,atmHB3},
		{atmCA,atmC},
		{atmC,atmOT1},
		{atmC,atmOT2}
};

topol_residue_atoms_bonds_angles_t residue_atoms_bonds_angles_ALA [] = {
		{atmN,atmC_,atmCA},
		//{atmHA,atmN,atmCA},
		//{atmHB1,atmCA,atmCB},
		//{atmHB2,atmCA,atmCB},
		//{atmHB3,atmCA,atmCB},
		{atmCA,atmN,atmC_},
		{atmC,atmCA,atmN},
		{atmCB,atmCA,atmN},
		{atmC,atmCA,atmCB},
		//{atmHN,atmN,atmC},
		{atmCA,atmO,atmC},
		{atmO,atmC,atmN_plus}
};

topol_residue_atoms_bonds_angles_t residue_atoms_bonds_angles_ALA_N_Terminal [] = {
		{atmN,atmCB,atmCA},
		//{atmHA,atmN,atmCA},
		//{atmH1,atmCA,atmN},
		//{atmH2,atmCA,atmN},
		//{atmH3,atmCA,atmN},
		//{atmHB1,atmCA,atmCB},
		//{atmHB2,atmCA,atmCB},
		//{atmHB3,atmCA,atmCB},
		{atmC,atmCA,atmN},
		{atmCB,atmC,atmCA},
		{atmCA,atmO,atmC},
		{atmO,atmC,atmN_plus}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_ALA [] = {
		{atmCB,atmC,atmCA,atmN}, //chi1
};

/*Represents the dihedral angle and its type of residue.*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ALA [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmN_plus, atmN, angl_psi}
};

/*Represents the dihedral angle and its type of residue when it is N-Terminal*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ALA_N_Terminal [] ={
	{atmN, atmC, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmN_plus, atmN, angl_psi}
};

topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ALA_C_Terminal [] ={
	{atmN, atmC, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmOT1, atmN, angl_psi} //atmN_plus
};
