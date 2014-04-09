/*Based on aminoacids.rtp from Gromacs. Force-field was Charmm27
 * We need to cite: Bjelkmar, Larsson, Cuendet, Hess, and Lindahl, JCTC 2010
 */

topol_residue_atoms_t residue_atoms_GLY_N_Terminal [] = {
                                                {atmN,  "N","NH1",-0.47},
                                                {atmCA, "CA","CT2",0.13},
                                                {atmC, "C","C",0.51},
                                                {atmO, "O","O",-0.51},
                                                //{atmHN, "HN","H", 0.31},
                                        		{atmH1, "H1","H",0.33},
                                        		{atmH2, "H2","H",0.33},
                                        		{atmH3, "H3","H",0.33},
                                                {atmHA1,"HA1","HB",0.09},
                                                {atmHA2,"HA2","HB",0.09}
                                           };

topol_residue_atoms_t residue_atoms_GLY [] = {
                                                {atmN,  "N","NH1",-0.47},
                                                {atmCA, "CA","CT2",-0.02},
                                                {atmC, "C","C",0.51},
                                                {atmO, "O","O",-0.51},
                                                {atmHN, "HN","H", 0.31},
                                                {atmHA1,"HA1","HB",0.09},
                                                {atmHA2,"HA2","HB",0.09}
                                           };

topol_residue_atoms_t residue_atoms_GLY_C_Terminal [] = {
                                                {atmN,  "N","NH1",-0.47},
                                                {atmCA, "CA","CT2",-0.02},
                                                {atmC, "C","C",0.51},
                                        		{atmOT1,"O","O",-0.67},
                                        		{atmOT2,"O","O",-0.67},
                                                {atmHN, "HN","H", 0.31},
                                                {atmHA1,"HA1","HB",0.09},
                                                {atmHA2,"HA2","HB",0.09}
                                           };

topol_residue_atoms_bonds_t residue_atoms_bonds_GLY_N_Terminal [] = {
		{atmC_,atmN}, //in aminoacids.rtp was wrote N+
		//{atmN, atmHN},
		//{atmN, atmH1},
		//{atmN, atmH2},
		//{atmN, atmH3},
		{atmN, atmCA},
		{atmCA,atmHA1},
		{atmCA,atmHA2},
		{atmCA, atmC},
		{atmC, atmO}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_GLY [] = {
		{atmC_,atmN}, //in aminoacids.rtp was wrote N+
		{atmN, atmHN},
		{atmN, atmCA},
		{atmCA,atmHA1},
		{atmCA,atmHA2},
		{atmCA, atmC},
		{atmC, atmO}
};

topol_residue_atoms_bonds_t residue_atoms_bonds_GLY_C_Terminal [] = {
		{atmC_,atmN}, //in aminoacids.rtp was wrote N+
		{atmN, atmHN},
		{atmN, atmCA},
		{atmCA,atmHA1},
		{atmCA,atmHA2},
		{atmCA, atmC},
		{atmC, atmOT1},
		{atmC, atmOT2}
};

topol_residue_atoms_bonds_angles_t residue_atoms_bonds_angles_GLY [] ={
		{atmN,atmC_,atmCA}, //atmCA_
		{atmHN,atmN,atmC_},
		{atmHA1,atmN,atmCA},
		{atmHA2,atmN,atmCA},
		{atmCA,atmN,atmC_},
		{atmC,atmCA,atmN},
		{atmO,atmC,atmCA},
		//{atmC,atmCA,atmCB}, This line was based on rotamers libray and it is different from Gromacs.
		{atmO,atmC,atmN_plus}
};


/*Represents the dihedral angle and its type of residue.*/
topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_GLY [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmN_plus, atmN, angl_psi}
};

topol_residue_atoms_dihedral_angles_type_t residue_atoms_dihedral_angles_type_GLY_C_Terminal [] ={
	{atmN, atmC_, atmCA, atmC, angl_phi},
	{atmC, atmCA, atmOT1, atmN, angl_psi} //atmN_plus
};
