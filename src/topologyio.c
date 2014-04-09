#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#include "functions.h"
#include "defines.h"
#include "protein.h"
#include "futil.h"
#include "topology_types.h"
#include "enums.h"
#include "topologylib.h"
#include "messages.h"
#include "string_owner.h"

#define MAX_LINE_FASTA 81
#define MAX_LEN_PROTEIN 1000

static void write_next_section(FILE *ftop){
	char next_section [] = "\n \n \n";
	fprintf(ftop,"%s",next_section);

}

static void write_general_topol_informations(FILE *ftop,const top_global_t *top){
	int charge = (int) top->protein_charge;
	fprintf(ftop,"[General information section] \n");
	fprintf(ftop,";num_atom  num_res num_bond_angles num_protein_side_chains protein_charge num_dihedral_angle_type \n");
	fprintf(ftop,"%i  %i %i %i %i %d \n",top->numatom, top->numres,
			top->number_bond_angles, top->number_protein_side_chains,
			charge,
			top->number_dihedral_angles_type);

}

static void write_atom_topol_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Atom section] \n");
	fprintf(ftop,";num_atom  atom_name  res_num  res_name charge\n");
	for (int a =0; a < (top->numatom);a++){
		fprintf(ftop,"%i     %5s     %i  %4s    %f\n",
				top->top_global_atom[a].atom_number,
				top->top_global_atom[a].atom_name,
				top->top_global_atom[a].res_number,
				top->top_global_atom[a].res_name,
				top->top_global_atom[a].charge);
	}
}

static void write_residues_atoms_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residue and its sequence atoms] \n");
	fprintf(ftop,";res_num  atom_first atom_end \n");
	for (int r=0; r<top->numres;r++){
		fprintf(ftop,"%i  %i    %i\n",top->top_global_res_atm[r].res_number,
				   top->top_global_res_atm[r].atom_first,top->top_global_res_atm[r].atom_last);
	}
}

static void write_residues_atoms_bond_length_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and their atom bond length] \n");
	fprintf(ftop,";res_num  atom_1 atom_2  bond_length\n");
	for (int i=0; i<top->numatom;i++){
		fprintf(ftop,"%i  %i    %i %f\n",top->top_global_res_atms_bond[i].res_number, top->top_global_res_atms_bond[i].atom_number1,
				top->top_global_res_atms_bond[i].atom_number2, top->top_global_res_atms_bond[i].bond_value);
	}

}

static void write_residues_atoms_bond_angles_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and their atom bond angle] \n");
	fprintf(ftop,";res_num  atom_1 atom_2  atom_3 angle_value\n");
	for (int i=0; i<top->number_bond_angles;i++){
		fprintf(ftop,"%i  %i  %i  %i   %f\n",top->top_global_res_atms_bond_angle[i].res_number,
				top->top_global_res_atms_bond_angle[i].atom_number1, top->top_global_res_atms_bond_angle[i].atom_number2,
				top->top_global_res_atms_bond_angle[i].atom_number3, top->top_global_res_atms_bond_angle[i].angle_value);
	}
}

static void write_residues_atoms_dihedral_phi_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and atoms for phi angle] \n");
	fprintf(ftop,";res_num  atom_1 atom_2  atom_3 atom_4\n");
	for (int i=0; i<top->numres;i++){
		fprintf(ftop,"%i  %i  %i  %i   %i\n",top->top_global_dieh_phi[i].res_number,
				top->top_global_dieh_phi[i].atom_number1, top->top_global_dieh_phi[i].atom_number2,
				top->top_global_dieh_phi[i].atom_number3, top->top_global_dieh_phi[i].atom_number4);
	}

}

static void write_residues_atoms_dihedral_omega_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and atoms for omega angle] \n");
	fprintf(ftop,";res_num  atom_1 atom_2  atom_3 atom_4\n");
	for (int i=0; i<top->numres;i++){
		fprintf(ftop,"%i  %i  %i  %i   %i\n",top->top_global_dieh_omega[i].res_number,
				top->top_global_dieh_omega[i].atom_number1, top->top_global_dieh_omega[i].atom_number2,
				top->top_global_dieh_omega[i].atom_number3, top->top_global_dieh_omega[i].atom_number4);
	}

}

static void write_residues_atoms_dihedral_psi_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and atoms for psi angle] \n");
	fprintf(ftop,";res_num  atom_1 atom_2  atom_3 atom_4\n");
	for (int i=0; i<top->numres;i++){
		fprintf(ftop,"%i  %i  %i  %i   %i\n",top->top_global_dieh_psi[i].res_number,
				top->top_global_dieh_psi[i].atom_number1, top->top_global_dieh_psi[i].atom_number2,
				top->top_global_dieh_psi[i].atom_number3, top->top_global_dieh_psi[i].atom_number4);
	}
}


static void write_residues_atoms_dihedral_side_chains_section(FILE *ftop,const top_global_t *top){
	fprintf(ftop,"[Residues and atoms for side chains angle] \n");
	fprintf(ftop,";res_num  chi atom_1 atom_2  atom_3 atom_4\n");
	for (int i=0; i<top->number_protein_side_chains;i++){
		fprintf(ftop,"%i  %i %i  %i  %i   %i\n",
				top->top_global_dieh_side_chains[i].res_number,
				top->top_global_dieh_side_chains[i].chi,
				top->top_global_dieh_side_chains[i].atom_number1,
				top->top_global_dieh_side_chains[i].atom_number2,
				top->top_global_dieh_side_chains[i].atom_number3,
				top->top_global_dieh_side_chains[i].atom_number4
				);
	}
}
static void write_residues_atoms_dihedral_angles_type(FILE *ftop,
		const top_global_t *top){
	char *str_type;
	str_type = Malloc(char,MAX_STR_DIHEDRAL_ANGLE_TYPE);
	fprintf(ftop,"[atoms and type of dihedral angles] \n");
	fprintf(ftop,";res_num atom_1 atom_2  atom_3 atom_4 type_dihedral\n");
	for (int i=0; i<top->number_dihedral_angles_type;i++){
		_type_of_diedhral_angle2str(str_type,
				&top->top_global_dihedral_angles_type[i].type_dihedral);
		fprintf(ftop,"%d %d  %d  %d  %d %s \n",
				top->top_global_dihedral_angles_type[i].res_number,
				top->top_global_dihedral_angles_type[i].atom_number1,
				top->top_global_dihedral_angles_type[i].atom_number2,
				top->top_global_dihedral_angles_type[i].atom_number3,
				top->top_global_dihedral_angles_type[i].atom_number4,
				str_type);


	}
	free(str_type);
}


void _save_topology_file(const char *path, const char *file_name,
		const top_global_t *top){
	FILE *ftop;
	char *fname = path_join_file(path,file_name);
	ftop = open_file(fname,fWRITE);

	write_general_topol_informations(ftop,top);
	write_next_section(ftop);
	write_atom_topol_section(ftop,top);
	write_next_section(ftop);
	write_residues_atoms_section(ftop,top);
	write_next_section(ftop);
	write_residues_atoms_dihedral_phi_section(ftop,top);
	write_next_section(ftop);
	write_residues_atoms_dihedral_psi_section(ftop,top);
	write_next_section(ftop);
	write_residues_atoms_dihedral_omega_section(ftop,top);
	write_next_section(ftop);		
	write_residues_atoms_dihedral_side_chains_section(ftop,top);
	write_next_section(ftop);
	/*
	write_residues_atoms_bond_length_section(ftop,top);
	write_next_section(ftop);
    write_residues_atoms_bond_angles_section(ftop,top);
	write_next_section(ftop);
	write_residues_atoms_dihedral_angles_type(ftop,top);
	*/
	free(fname);
	fclose(ftop);
}


static boolean_t _check_pdb_fasta_file(char *line){
	/*This function check if line contains pdbid which means that
	 * the file is pdb fasta file
	 * Example line: >1VII:A|PDBID|CHAIN|SEQUENCE
	 */
	char *pch;
	pch = strtok (line,"|");
	while (pch != NULL){
		if (strcmp(pch,"PDBID") == 0){
			return btrue;
		}
		pch = strtok (NULL, "|");
	}
	return bfalse;

}

void _create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein){
	/* Create a pdb fasta file.
	 * Example of pdb fasta file:
	 * >1VII:A|PDBID|CHAIN|SEQUENCE
     * MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF
	*/
	FILE *arq;
	int seq_index = -1;
	char *prot_seq_aux;
	int len = 0;
	char header[MAX_LINE_FASTA];

	arq = open_file(file_name_protein, fWRITE);
	//Writing header
	strcpy(header,">");
	strcat(header,prot_name);
	strcat(header,":");
	strcat(header,chain_name);
	strcat(header,"|PDBID|CHAIN|SEQUENCE");
	fprintf (arq, "%s\n", header);
	len = strlen(prot_seq);
	prot_seq_aux = Malloc(char, MAX_LINE_FASTA+1);

	strncpy(prot_seq_aux,prot_seq,MAX_LINE_FASTA);
	fprintf (arq, "%s\n",prot_seq_aux);
	/*
	while (seq_index < len){
		seq_index++;
		//breaking in lines
		while ( (seq_index < MAX_LINE_FASTA) &&
				seq_index < len){
			seq_index++;
			strcat(prot_seq_aux,&prot_seq[seq_index]);
		}
		fprintf (arq, "%s\n",prot_seq_aux);
	}*/
    free(prot_seq_aux);
	fclose(arq);

}

amino_t *_load_fasta_pdb(const char *file_name_protein, int *n_residues,
		int *numatom, int *num_bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles_type, boolean_t *has_his){
	/* load amino_t based on pdb fasta file.
	 * Example of pdb fasta file:
	 * >1VII:A|PDBID|CHAIN|SEQUENCE
     * MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF
	*/
	FILE *arq;
	amino_t * seq_prim;
	int fscanfError;
	type_aminos_t amino_id;
    int primary_seq_index;
    const topol_residues_t *residue_ff;
    char line[MAX_LINE_FASTA+1]; //represents Fasta lines
    char seq_line[MAX_LEN_PROTEIN];//represents the primary sequence of protein
    char c;
    int i;
    boolean_t first_line = btrue;
    boolean_t read_fasta_file = btrue;

	arq = open_file(file_name_protein, fREAD);
	fgets(line,MAX_LINE_FASTA,arq);
	if (_check_pdb_fasta_file(line) == bfalse){
		char msg[300];
	    sprintf(msg,"Fasta file %s is not a pdb fasta file \n",file_name_protein);
		fatal_error(msg);
	}
	//After checking Fasta file, initialize variables
	*numatom = 0; //Initialize number of atoms of protein
	*num_bond_angles = 0;//Initialize number of bond angles of protein
	*num_protein_side_chains = 0; //Initialize number of protein chains
	*num_dihedral_angles_type = 0; //Initialize number of dihedral angles of protein
	primary_seq_index = 0;
	/*Getting the primary sequence. It is represented by seq_line variable.
	 * So it contains residues from one chain of protein only. Because of this
	 * is checked if line starts with >*/
	fgets(line,MAX_LINE_FASTA,arq); //Here line contains residues of protein. Maybe it is first part like 1E8A
	trim(line);
	remove_character(line, '\n');
	strcpy(seq_line,line);
	while ( (fgets(line,MAX_LINE_FASTA,arq) != NULL) &&
			(read_fasta_file == btrue) ){
		remove_character(line,'\n');
		trim(line);
		if (strncmp(line,">",1) == 0){ //Force to work with one sequence only
			read_fasta_file = bfalse;
		}else{
			strcat(seq_line,line);
		}
	}
	fclose(arq);
	remove_character(seq_line, '\n');
	*n_residues = strlen(seq_line);
	seq_prim = allocateAmino(*n_residues);
	//Malloc(amino_t, *n_residues);
	for (i = 0; i < *n_residues; i++){
		c = seq_line[i];
		amino_id = get_amino_id(c);

		if (primary_seq_index == 0){//N_Terminal
			residue_ff = get_residue_from_topol_N_Terminal(amino_id);
		}else if ( (primary_seq_index > 0) && (primary_seq_index < *n_residues-1) ) {
			residue_ff = get_residue_from_topol(amino_id);
		}else{//C-Terminal
			residue_ff = get_residue_from_topol_C_Terminal(amino_id);
		}
		set_numatom_from_topol(numatom,&residue_ff->nr_atoms); //set the number of atoms from each amino acids of primary sequence.
		set_amino_from_topol(seq_prim,primary_seq_index,residue_ff);
		set_num_bond_angle_from_topol(num_bond_angles,&residue_ff->number_bond_angle);
		set_num_side_chains_from_topol(num_protein_side_chains,&residue_ff->nr_side_chains);
		set_num_dihedral_angles(num_dihedral_angles_type, &residue_ff->number_dihedral_angle_type);
		set_has_his(has_his,&amino_id);
        primary_seq_index = primary_seq_index + 1;
	}
	return seq_prim;

}

amino_t *_load_protein2(const char *file_name_protein, int *n_residues,
		int *numatom, int *num_bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles_type){
	FILE *arq;
	amino_t * seq_prim;
	int fscanfError;
	type_aminos_t amino_id;
    int primary_seq_index;
    char c;
    const topol_residues_t *residue_ff;

	arq = open_file(file_name_protein, fREAD);
	fscanfError = fscanf(arq, "%d\n", n_residues); //get number of residues from sequence file
	seq_prim = Malloc(amino_t, *n_residues);
	*numatom = 0; //Initialize number of atoms of protein
	*num_bond_angles = 0;//Initialize number of bond angles of protein
	*num_protein_side_chains = 0; //Initialize number of protein chains
	*num_dihedral_angles_type = 0; //Initialize number of dihedral angles of protein
	primary_seq_index = 0;

	fscanfError = fscanf(arq, "%c", &c); //starts to read protein sequence
	while (!feof(arq) && ( c != ' ' && c != '\n')) {
		amino_id = get_amino_id(c);
		if (primary_seq_index == 0){//N_Terminal
			residue_ff = get_residue_from_topol_N_Terminal(amino_id);
		}else if ( (primary_seq_index > 0) && (primary_seq_index < *n_residues-1) ) {
			residue_ff = get_residue_from_topol(amino_id);
		}else{//C-Terminal
			residue_ff = get_residue_from_topol_C_Terminal(amino_id);
		}
		set_numatom_from_topol(numatom,&residue_ff->nr_atoms); //set the number of atoms from each amino acids of primary sequence.
		set_amino_from_topol(seq_prim,primary_seq_index,residue_ff);
		set_num_bond_angle_from_topol(num_bond_angles,&residue_ff->number_bond_angle);
		set_num_side_chains_from_topol(num_protein_side_chains,&residue_ff->nr_side_chains);
		set_num_dihedral_angles(num_dihedral_angles_type, &residue_ff->number_dihedral_angle_type);
        primary_seq_index = primary_seq_index + 1;
        fscanfError = fscanf(arq,"%c", &c);
	}
	fclose(arq);
	return seq_prim;
}

void _show_protein_backbone(const protein_backbone_t *prot_back,
		const top_global_t *top){
	char msg [50];
    display_msg("Residue Carbon Alfa_Carbon Nitrogen Carbon_previous Nitrogen_next \n");
    for (int i = 1; i <= top->numres;i++){
    	sprintf(msg,"%d %d %d %d %d %d",prot_back[i-1].res_number,
    			prot_back[i-1].atom_C,
    			prot_back[i-1].atom_Ca,
    			prot_back[i-1].atom_N,
    			prot_back[i-1].atom_C_,
    			prot_back[i-1].atom_N_plus );
    	strcat(msg,"\n");
    	display_msg(msg);
    }

}



