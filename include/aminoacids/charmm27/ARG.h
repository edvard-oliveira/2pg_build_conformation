/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_ARG_N_Terminal [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.21},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmH1, "H1","H",0.33},
		{atmH2, "H2","H",0.33},
		{atmH3, "H3","H",0.33},
//		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.10},
		{atmCB, "CB","CT2",-0.18},
		{atmCG, "CG","CT2",-0.18},
		{atmCD, "CD","CT2",0.20},
		{atmNE, "NE","NC2",-0.70},
		{atmCZ, "CZ","C",0.64},
		{atmNH1, "NH1","NC2",-0.80},
		{atmNH2, "NH2","NC2",-0.80},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2,"HB2","HA",0.09},
		{atmHG1, "HG1","HA",0.09},
		{atmHG2, "HG2","HA",0.09},
		{atmHD1, "HD1","HA",0.09},
		{atmHD2, "HD2","HA",0.09},
		{atmHE, "HE","HC",0.44},
		{atmHH11, "HH11","HC",0.46},
		{atmHH12, "HH12","HC",0.46},
		{atmHH21, "HH21","HC",0.46},
		{atmHH22, "HH22","HC",0.46}
		//{atmOT1,"OT1","OT1",-0.51}
};

topol_residue_atoms_t residue_atoms_ARG [] = {
		{atmN, "N","NH1",-0.47},
		{atmCA, "CA","CT1",0.07},
		{atmC,"C","C",0.51},
		{atmO,"O","O",-0.51},
		{atmHN, "HN","H",0.31},
		{atmHA, "HA","HB",0.09},
		{atmCB, "CB","CT2",-0.18},
		{atmCG, "CG","CT2",-0.18},
		{atmCD, "CD","CT2",0.20},
		{atmNE, "NE","NC2",-0.70},
		{atmCZ, "CZ","C",0.64},
		{atmNH1, "NH1","NC2",-0.80},
		{atmNH2, "NH2","NC2",-0.80},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2,"HB2","HA",0.09},
		{atmHG1, "HG1","HA",0.09},
		{atmHG2, "HG2","HA",0.09},
		{atmHD1, "HD1","HA",0.09},
		{atmHD2, "HD2","HA",0.09},
		{atmHE, "HE","HC",0.44},
		{atmHH11, "HH11","HC",0.46},
		{atmHH12, "HH12","HC",0.46},
		{atmHH21, "HH21","HC",0.46},
		{atmHH22, "HH22","HC",0.46}
		//{atmOT1,"OT1","OT1",-0.51}
};

topol_residue_atoms_t residue_atoms_ARG_C_Terminal [] = {
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
		{atmNE, "NE","NC2",-0.70},
		{atmCZ, "CZ","C",0.64},
		{atmNH1, "NH1","NC2",-0.80},
		{atmNH2, "NH2","NC2",-0.80},
		{atmHB1, "HB1","HA",0.09},
		{atmHB2,"HB2","HA",0.09},
		{atmHG1, "HG1","HA",0.09},
		{atmHG2, "HG2","HA",0.09},
		{atmHD1, "HD1","HA",0.09},
		{atmHD2, "HD2","HA",0.09},
		{atmHE, "HE","HC",0.44},
		{atmHH11, "HH11","HC",0.46},
		{atmHH12, "HH12","HC",0.46},
		{atmHH21, "HH21","HC",0.46},
		{atmHH22, "HH22","HC",0.46}

};

topol_residue_atoms_bonds_t residue_atoms_bonds_ARG [] = {
		{atmC_,atmN}, //in aminoacids.rtp was wrote N+
		{atmCA,atmCB,},
		{atmCB,atmCG},
		{atmCG,atmCD},
		{atmCD,atmNE},
		{atmNE,atmCZ},
		//{atmCZ,atmNH2},
		{atmN,atmHN},
		{atmN,atmCA},
		{atmCA,atmC},
		{atmCA,atmHA},
		{atmCB,atmHB1},
		{atmCB,atmHB2},
		{atmCG,atmHG1},
		{atmCG,atmHG2},
		{atmCD,atmHD1},
		{atmCD,atmHD2},
		{atmNE,atmHE},
		{atmNH1,atmHH11},
		{atmNH1,atmHH12},
		{atmNH2,atmHH21},
		{atmNH2,atmHH22},
		{atmOT1,atmC},
		{atmC,atmO}
		//{atmCZ,atmNH1}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_ARG_N_Terminal [] = {
		{atmC,atmN}, //in aminoacids.rtp was wrote N+
		{atmCA,atmCB,},
		{atmCB,atmCG},
		{atmCG,atmCD},
		{atmCD,atmNE},
		{atmNE,atmCZ},
		//{atmCZ,atmNH2},
		//{atmN,atmHN},
		{atmN_plus,atmCA},
		{atmCA,atmC},
		//{atmCA,atmHA},
		//{atmCB,atmHB1},
		//{atmCB,atmHB2},
		//{atmCG,atmHG1},
		//{atmCG,atmHG2},
		//{atmCD,atmHD1},
		//{atmCD,atmHD2},
		//{atmNE,atmHE},
		//{atmNH1,atmHH11},
		//{atmNH1,atmHH12},
		//{atmNH2,atmHH21},
		//{atmNH2,atmHH22},
		{atmC,atmO}
		//{atmCZ,atmNH1}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_ARG_C_Terminal [] = {
		{atmC_,atmN}, //in aminoacids.rtp was wrote N+
		{atmCA,atmCB},
		{atmCB,atmCG},
		{atmCG,atmCD},
		{atmCD,atmNE},
		{atmNE,atmCZ},
		{atmCZ,atmNH2},
		{atmN,atmHN},
		{atmN,atmCA},
		{atmCA,atmC},
		{atmCA,atmHA},
		{atmCB,atmHB1},
		{atmCB,atmHB2},
		{atmCG,atmHG1},
		{atmCG,atmHG2},
		{atmCD,atmHD1},
		{atmCD,atmHD2},
		{atmNE,atmHE},
		{atmNH1,atmHH11},
		{atmNH1,atmHH12},
		{atmNH2,atmHH21},
		{atmNH2,atmHH22},
		{atmC,atmOT1},
		{atmC,atmOT2}
		//{atmCZ,atmNH1}
};
topol_residue_atoms_bonds_angles_t residue_atoms_bonds_angles_ARG [] = {
	{atmN,atmC_,atmCA}, //atmCA_
	{atmHN,atmN,atmCA},
	{atmHA,atmN,atmCA},
	{atmCA,atmN,atmC_},
	{atmCB,atmCA,atmN},
	{atmCG,atmCB,atmCA},
	{atmCD,atmCG,atmCB},
	{atmNE,atmCD,atmCG},
	{atmCZ,atmNE,atmCD},
	{atmNH2,atmCZ,atmNE},
	{atmNH1,atmCZ,atmNE},
	{atmNH1,atmCZ,atmNH2},
	{atmC,atmCA,atmN},
	{atmC,atmCA,atmCB},
	{atmO,atmC,atmCA},
	{atmHB1,atmCB,atmCA},
	{atmHB2,atmCB,atmCA},
	{atmHG1,atmCG,atmCB},
	{atmHG2,atmCG,atmCB},
	{atmHD1,atmCD,atmCG},
	{atmHD2,atmCD,atmCG},
	{atmHE,atmNE,atmCD},
	{atmHH11,atmNH1,atmCZ},
	{atmHH12,atmNH1,atmCZ},
	{atmHH21,atmNH2,atmCZ},
	{atmHH22,atmNH2,atmCZ},
	{atmO,atmC,atmN_plus}
};

topol_residue_atoms_bonds_angles_t residue_atoms_bonds_angles_ARG_C_Terminal [] = {
	{atmN,atmC_,atmCA},
	{atmCA,atmN,atmC_},
	{atmCB,atmCA,atmN},
	{atmCG,atmCB,atmCA},
	{atmCD,atmCG,atmCB},
	{atmNE,atmCD,atmCG},
	{atmCZ,atmNE,atmCD},
	{atmNH2,atmCZ,atmNE},
	{atmNH1,atmCZ,atmNE},
	{atmNH1,atmCZ,atmNH2},
	{atmC,atmCA,atmN},
	{atmC,atmCA,atmCB},
	{atmOT1,atmC,atmCA},
	{atmHN,atmN,atmC},
	{atmHA,atmN,atmCA},
	{atmHB1,atmCB,atmCA},
	{atmHB2,atmCB,atmCA},
	{atmHG1,atmCG,atmCB},
	{atmHG2,atmCG,atmCB},
	{atmHD1,atmCD,atmCG},
	{atmHD2,atmCD,atmCG},
	{atmHE,atmNE,atmCD},
	{atmHH11,atmNH1,atmCZ},
	{atmHH12,atmNH1,atmCZ},
	{atmHH21,atmNH2,atmCZ},
	{atmHH22,atmNH2,atmCZ},
	{atmOT2,atmC,atmOT1}
//	{atmOT1,atmC,atmN_plus}
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_rot_ARG [] = {
		{atmN, atmCA, atmCB, atmCG}, //chi1
		{atmCA, atmCB, atmCG, atmCD}, //chi2
		{atmCB, atmCG, atmCD, atmNE}, //chi3
		{atmCG, atmCD, atmNE, atmCZ}, //chi4
//		{atmCD, atmNE, atmCZ, atmNH1} //chi5
};


/*Represents the dihedral angle and its type of residue.*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ARG [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmN_plus, atmN, angl_psi},
	{atmNE,atmNH2,atmNH1,atmCZ, angl_typ_dieh_1} //VERFICAR SE ESTA CORRETO
};

/*Represents the dihedral angle and its type of residue when it is C-Terminal.*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_ARG_C_Terminal [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmOT1, atmN, angl_psi}, //atmN_plus
	{atmNE, atmNH2, atmNH1, atmCZ, angl_typ_dieh_1} //VERFICAR SE ESTA CORRETO
};
