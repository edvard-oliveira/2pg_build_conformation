#include <iostream>
#include <fstream>
#include <string.h>
#include <malloc.h>

#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "functions.h"
#include "messages.h"
#include "futil.h"
#include "string_owner.h"
#include "math_owner.h"
#include "pdbatom.h"
#include "pdbio.h"
#include "topology.h"
#include "topologylib.h"
#include "messages.h"

type_fitness_energies_t str2type_fitness_energies(char *name_fitness_energies){
	if (strcmp(name_fitness_energies,"Potential") == 0){
		return fit_ener_potential;
	}else if (strcmp(name_fitness_energies,"Van_der_Waals") == 0){
		return fit_ener_edw;
	}else if (strcmp(name_fitness_energies,"Electrostatic") == 0){
		return fit_ener_ele;
	}else if (strcmp(name_fitness_energies,"Hydrophobic") == 0){
		return fit_hydrophobic;
	}else if (strcmp(name_fitness_energies,"Hydrophilic") == 0){
		return fit_hydrophilic;
	}else if (strcmp(name_fitness_energies,"Total_Area") == 0){
		return fit_total_area;
	}else if (strcmp(name_fitness_energies,"Gyrate") == 0){
		return fit_gyrate;
	}else if (strcmp(name_fitness_energies,"H_Bond") == 0){
		return fit_hbond;
	}else if (strcmp(name_fitness_energies,"H_Bond_Main") == 0){
		return fit_hbond_main;
	}else if (strcmp(name_fitness_energies,"GBSA_Solvatation") == 0){
		return fit_GBSA_Solvatation;
	}else if (strcmp(name_fitness_energies,"Stride_total") == 0){
		return fit_stride_total;
	}else if (strcmp(name_fitness_energies,"Stride_helix") == 0){
		return fit_stride_helix;
	}else if (strcmp(name_fitness_energies,"Stride_beta") == 0){
		return fit_stride_beta;
	}else{
		char msg[1024];
		sprintf(msg,"The option %s did not find at type_fitness_energies_t. Please, check it!", name_fitness_energies);
		fatal_error(msg);
	}

}

void type_fitness_energies2str(char *name_fitness_energies,
		const type_fitness_energies_t *type_fitness ){
	if (*type_fitness == fit_ener_potential){
		strcpy(name_fitness_energies,"Potential");
	}else if (*type_fitness == fit_ener_edw){
		strcpy(name_fitness_energies,"Van_der_Waals");
	}else if (*type_fitness == fit_ener_ele){
		strcpy(name_fitness_energies,"Electrostatic");
	}else if (*type_fitness == fit_hydrophobic){
		strcpy(name_fitness_energies,"Hydrophobic");
	}else if (*type_fitness == fit_hydrophilic){
		strcpy(name_fitness_energies,"Hydrophilic");
	}else if (*type_fitness == fit_total_area){
		strcpy(name_fitness_energies,"Total_Area");
	}else if (*type_fitness == fit_gyrate){
		strcpy(name_fitness_energies,"Gyrate");
	}else if (*type_fitness == fit_hbond){
		strcpy(name_fitness_energies,"H_Bond");
	}else if (*type_fitness == fit_hbond_main){
		strcpy(name_fitness_energies,"H_Bond_Main");
	}else if (*type_fitness == fit_GBSA_Solvatation){
		strcpy(name_fitness_energies,"GBSA_Solvatation");
	}else if (*type_fitness == fit_stride_total){
		strcpy(name_fitness_energies,"Stride_total");
	}else if (*type_fitness == fit_stride_helix){
		strcpy(name_fitness_energies,"Stride_helix");
	}else if (*type_fitness == fit_stride_beta){
		strcpy(name_fitness_energies,"Stride_beta");
	}else{
		char msg[1024];
		sprintf(msg,"The option did not find at type_fitness_energies2str function. Please, check it!");
		fatal_error(msg);
	}
}

void readPopulationFile(protein **pop, const int *innumberpopulation,
		const char *path, const char *file_name_pop,
		const amino_t *primary_seq, const int *nresiduos, 
		const top_global_t *top_global) {
	/*Loads the population file*/
    char chBloco[TAM_BLOCO];
    float phi,psi, omega;
    int number_late, pos;
    float late1,late2,late3,late4;
    FILE *arq_name_pop;
    int res;
    int rot;
    char *fname = path_join_file(path,file_name_pop);
    arq_name_pop = open_file( fname ,fREAD); //path_file_name_pop
    //Checks if population file is empty
    if ( file_is_empty(arq_name_pop) == btrue) {
    	fatal_error("Population file is empty. Check it!!!");
    }
    fgets(chBloco, TAM_BLOCO, arq_name_pop);
    for (int i = 0; i < *innumberpopulation; i++){
    	//Remove Line that contain the ## character
        if (strcmp(chBloco,"##\n") == 0) {
           fgets(chBloco, TAM_BLOCO, arq_name_pop);
        }
        pop[i]->nr_residues = *nresiduos;
	//free(pop[i]->residuo);
	//deAllocateAmino(pop[i]->residuo, *nresiduos);
	//char *c = (char *)malloc(sizeof(char)*5);
	//free(c);
        //pop[i]->residuo = allocateAmino(*nresiduos);
        //Malloc(amino_t, *nresiduos);
        for (res = 0; res < pop[i]->nr_residues; res++ ){
        	sscanf(chBloco,"%f %f %f %i %i %f %f %f %f",&phi,&psi,&omega,
        			&number_late,&pos,&late1,&late2,&late3,&late4);
		 
	        pop[i]->residuo[res].phi         = degree2radians(&phi);
            pop[i]->residuo[res].psi         = degree2radians(&psi);
            pop[i]->residuo[res].omega       = degree2radians(&omega);
            pop[i]->residuo[res].number_late = number_late;
	        pop[i]->residuo[res].pos_late = pos;
		    pop[i]->residuo[res].id =  primary_seq[res].id;
	        strcpy(pop[i]->residuo[res].aminoacido, primary_seq[res].aminoacido);
            pop[i]->residuo[res].late = Malloc(float,pop[i]->residuo[res].number_late);
            for (rot = 0; rot < pop[i]->residuo[res].number_late; rot++){
                if (rot == 0)
                 pop[i]->residuo[res].late[rot] = degree2radians(&late1);
                if (rot == 1)
                 pop[i]->residuo[res].late[rot] = degree2radians(&late2);
                if (rot == 2)
                 pop[i]->residuo[res].late[rot] = degree2radians(&late3);
                if (rot == 3)
                 pop[i]->residuo[res].late[rot] = degree2radians(&late4);
            }
            fgets(chBloco, TAM_BLOCO, arq_name_pop);
        }//end of angle values of residues
      	//Building z_matrix
       	build_z_matrix(pop[i]->z_matrix,top_global);
        //Remove Line that contain the $$ character
        fgets(chBloco, TAM_BLOCO, arq_name_pop);
        trim(chBloco);
        if (is_equal(chBloco,"$$") == btrue){
        	fgets(chBloco, TAM_BLOCO, arq_name_pop);
        }
    }
    free(fname);
    fclose(arq_name_pop);
}


int get_model_number(const char *s){
	/* Returns the number of MODEL presented in PDB file, for example.
	 * In s contains MODEL        1. The function will return 1.
	*/
	 char *aux, *token;
	 int n;

	 aux = Malloc(char, strlen(s)+1);
	 strcpy(aux, s);
	 trim(aux);

	 token = strtok(aux," ");
	 token = strtok(NULL," "); //Skip MODEL

	 n = str2int(token); //Obtaing the number of MODEL

	 free(aux);

	 return n;
}

void set_information_allocate_side_chains(protein *individual, const amino_t *primary_seq, 
	const top_global_t *top_global){
	/* This function set information and allocate side chains based on primary sequence of 
	 * protein.
	*/
	type_aminos_t amino_id;
	char *aux;
    aux = Malloc(char, 4);
	for (int res = 0; res < top_global->numres; res++){
		strcpy(aux, primary_seq[res].idl3);
		amino_id = _get_amino_id_3(aux);
       	individual->residuo[res].pos_late = 0;
    	individual->residuo[res].id =  primary_seq[res].id;
       	strcpy(individual->residuo[res].aminoacido, primary_seq[res].aminoacido);
		if (_has_side_chain(&amino_id) == btrue){
			individual->residuo[res].number_late = _number_side_chains(&amino_id);
			individual->residuo[res].late = Malloc(float, individual->residuo[res].number_late);
			for (int chi = 0; chi < individual->residuo[res].number_late; chi++){
				individual->residuo[res].late[chi] = 0;
			}
		}
	}
	free(aux);
}

void readPopulationFile_PDB(char *chBloco, FILE *arq_name_pop, protein **pop, 
	const int *innumberpopulation, const amino_t *primary_seq, 	
	const top_global_t *top_global){
	/* This function load a population file that is represented by Cartesian (PDB)
	*/
	fatal_error("I have to check the use of the function readPopulationFile_PDB");
	int line_atm;
	int num_model, ind;
	pdb_atom_t *pdb_atoms;

	pdb_atoms = allocate_pdbatom(&top_global->numatom);
	ind = -1;	

	do {		
		if (strncmp(chBloco,"MODEL",5) == 0){
			line_atm = -1;
			num_model = get_model_number(chBloco);
		}else if (strncmp(chBloco,"ATOM",4) == 0){
			line_atm = line_atm + 1;			
			load_pdb_atoms(chBloco, pdb_atoms, &line_atm);
		}else if (strncmp(chBloco,"ENDMDL",6) == 0){ //End of ATOM Section of MODEL
			ind = ind + 1;
			set_information_allocate_side_chains(pop[ind], primary_seq, top_global);
			pdbatoms2protein(pop[ind], pdb_atoms, top_global);			
			//build_z_matrix_from_pdb(pop[ind]->z_matrix, top_global, pdb_atoms);
		}		
	}while ( fgets(chBloco,TAM_BLOCO,arq_name_pop) != NULL);	

	desAllocate_pdbatom(pdb_atoms);	
	//Checking the number of models and number of population
	if (num_model != *innumberpopulation){
		char msg[1024];
		sprintf(msg,"The number of models was %i while the number of individuals is %i. These values have to be equals. Please, check it!", num_model, *innumberpopulation);
		fatal_error(msg);		
	}
}



protein** allocatePopulation(int inPopSize, int nresiduos, int nr_fitness, int numatom){
	 protein **population;
     population = Malloc(protein*,inPopSize);
     for(int i=0;i<inPopSize;i++){
    	population[i] = allocateProtein(nresiduos,nr_fitness, numatom);
 	}
    return population;
}

void deAllocatePopulation(protein** population, int inPopSize){
	//printf("I have to implement deAllocatePopulation \n");
	
    for(int i=0;(i<inPopSize);i++){
    	if (population[i] != NULL){
    		deAllocateProtein(population[i]);
    	}
	}
    free(population);
    
}


type_aminos_t get_amino_id(char c){
/*receives an amino (char) and returns its id*/
	    type_aminos_t amino_id;
		switch (c)   {
			case 'A': //alanina
				amino_id =  aALA;
				break;
			case 'V': //valina
				amino_id =  aVAL;
				break;
			case 'F': //fenilalanina
				amino_id =  aPHE;
				break;
			case 'P': //prolina
				amino_id =  aPRO;
				break;
			case 'L': //leucina
				amino_id =  aLEU;
				break;
			case 'I': //iso leucina
				amino_id =  aILE;
				break;
			case 'R': //argenina
				amino_id =  aARG;
				break;
			case 'D': //acido aspartico
				amino_id =  aASP;
				break;
			case 'E': //acido glutamico
				amino_id =  aGLU;
				break;
			case 'S': //serina
				amino_id =  aSER;
				break;
			case 'T': //treonina
				amino_id =  aTHR;
				break;
			case 'C': //cisteina
				amino_id =  aCYS;
				break;
			case 'N': //asparagina
				amino_id =  aASN;
				break;
			case 'Q': //glutanima
				amino_id =  aGLN;
				break;
			case 'H': //histidina
				amino_id =  aHIS;
				break;
			case 'K': //lisina
				amino_id =  aLYS;
				break;
			case 'Y': //tirosina
				amino_id =  aTYR;
				break;
			case 'M': //metionina
				amino_id =  aMET;
				break;
			case 'W': //triptofano
				amino_id =  aTRP;
				break;
			case 'G': //glicina
				amino_id =  aGLY;
				break;
			default:
				amino_id = aNR;
				char mens [] = "Amino value is not known. Please, check your file.";
				fatal_error(mens);
		}
		return amino_id;
}

void build_pdb_file_name(char *pdb_file_name, const char *aux_name,
		const char *__restrict prefix){
	strcpy(pdb_file_name, prefix);
	strcat(pdb_file_name,aux_name);
	strcat(pdb_file_name,".pdb");
}

void set_population_file_name_with_pop_file_name(char *pop_file_name,
		const int *generation){
	sprintf(pop_file_name,"pop_%d.pop",*generation );
}


void set_do_dssp_percentage(float *percentage_helix, float *percentage_beta,
		float *percentage_total, const char *path_local_execute,
		const char *file_name){
	/* This function reads outputfile of do_dssp program.
	 * It calculates the percentage of secondary structure.
	 * It was developed by Leandro Oliveira Bortot.
	 */
	FILE *input;
	char *file_compute;

	char ss ;
	float helix_count , helix_fraction ;
	float beta_count , beta_fraction ;
	float ss_count , ss_fraction ;
	float n;

	file_compute = path_join_file(path_local_execute,file_name);

	// opens the input file
	input = open_file(file_compute,fREAD);

	ss_count = 0;
	helix_count = 0;
	beta_count = 0;
	n = -1;	// number of residues. -1 because \0 is accounted for

	while( fscanf(input,"%c",&ss) != EOF ){		// reads the whole file char by char

		if( ss == 'H' ) helix_count++;	// alpha-helix
		//if( ss == 'G' ) helix_count++;	// 3/10 helix
		//if( ss == 'I' ) helix_count++;	// pi-helix

		if( ss == 'B' ) beta_count++;	// beta-bridge
		if( ss == 'E' ) beta_count++;	// beta-sheet


		n++;
	}

	ss_count = helix_count + beta_count ;	// total SS

	*percentage_helix = helix_count / n;
	*percentage_beta = beta_count / n;
	*percentage_total = ss_count / n;

	//printf("%.4f\t%.4f\t%.4f", helix_fraction , beta_fraction , ss_fraction );

	fclose(input);
	free(file_compute);

}

/*
 amino_t *load_protein(char *file_name_protein, int *nresiduos, int *numatom) {
		// Build the primary sequence of the protein which is used a lot of functions.
		FILE *arq;
		amino_t * prot_prin;
		arq = open_file(file_name_protein, fREAD);

		char c;
		int i = 0;
		int j;
		int find;
		int fscanfError;
		int n_typesAmino = 1;
		*numatom = 0;
		fscanfError = fscanf(arq, "%d\n", nresiduos);
		prot_prin = Malloc(samino, *nresiduos);
		i = 0;
		fscanfError = fscanf(arq, "%c", &c);
		while (!feof(arq) && ( c != ' ' && c != '\n')) {
			find = search_amino(c, n_typesAmino);
			if (find != -1) {
				aminoacidos[n_typesAmino - 1] = c;
				n_typesAmino++;
			}
			prot_prin[i].aminoacido = c;
			prot_prin[i].phi = 0;
			prot_prin[i].psi = 0;
			switch ( c )   {
				case 'A': //alanina
					prot_prin[i].number_late = 1;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					prot_prin[i].late[0] = 180;
					*numatom += 4;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "ALA");
					break;
				case 'V': //valina
					prot_prin[i].number_late = 1;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 10;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "VAL");
					break;
				case 'F': //fenilalanina
					prot_prin[i].number_late = 2;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 14;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "PHE");
					break;
				case 'P': //prolina
					prot_prin[i].number_late = 2;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 9;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "PRO");
					break;
				case 'L': //leucina
					prot_prin[i].number_late = 2;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++){
						prot_prin[i].late[j] = 180;
					}
					*numatom += 13;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "LEU");
					break;
				case 'I': //iso leucina
					prot_prin[i].number_late = 2;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 13;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "ILE");
					break;
				case 'R': //argenina
					prot_prin[i].number_late = 4;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 18;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "ARG");
					break;
				case 'D': //acido aspartico
					prot_prin[i].number_late = 2;//sp2
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 6;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "ASP");
					break;
				case 'E': //acido glutamico
					prot_prin[i].number_late = 3;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 9;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "GLU");
					break;
				case 'S': //serina
					prot_prin[i].number_late = 1;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 5;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "SER");
					break;
				case 'T': //treonina
					prot_prin[i].number_late = 1;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 8;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "THR");
					break;
				case 'C': //cisteina
					prot_prin[i].number_late = 1;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 4;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "CYS");
					break;
				case 'N': //asparagina
					prot_prin[i].number_late = 2;//sp2
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 8;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "ASN");
					break;
				case 'Q': //glutanima
					prot_prin[i].number_late = 3;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 11;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "GLN");
					break;
				case 'H': //histidina
					prot_prin[i].number_late = 2;//sp2
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 12;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "HIS");
					break;
				case 'K': //lisina
					prot_prin[i].number_late = 4;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 16;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "LYS");
					break;
				case 'Y': //tirosina
					prot_prin[i].number_late = 2;//sp2
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 15;
					prot_prin[i].HP = 0;
					strcpy(prot_prin[i].idl3, "TYR");
					break;
				case 'M': //metionina
					prot_prin[i].number_late = 3;
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 11;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "MET");
					break;
				case 'W': //triptofano
					prot_prin[i].number_late = 2;//sp2
					prot_prin[i].late = Malloc(float, prot_prin[i].number_late);
					for (j = 0; j < prot_prin[i].number_late; j++) {
						prot_prin[i].late[j] = 180;
					}
					*numatom += 18;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "TRP");
					break;
				case 'G': //glicina
					prot_prin[i].number_late = 0;
					prot_prin[i].late = NULL;
					*numatom += 1;
					prot_prin[i].HP = 1;
					strcpy(prot_prin[i].idl3, "GLY");
					break;
			}
			if (i == 0) {
				*numatom += 8;
			} else {
				if (i == (*nresiduos - 1)) {
					*numatom += 7;
				} else {
					if (i < *nresiduos)
						switch (c) {
						case 'P':
							*numatom += 5;
							break;
						default:
							*numatom += 6;
					   }
				}
			 }
			 i++;
			 fscanfError = fscanf(arq,"%c", &c);
	 }
	fclose(arq);
	return prot_prin;
}
 *
 *
 */



/*
 void set_seqtyp(int seqtyp[], const amino_t *prot_prin, int nresiduo){
	for (int i=0; i<nresiduo;i++){
		switch (prot_prin[i].aminoacido)   {
			case 'A': //alanina
				seqtyp[i] = aALA;
				break;
			case 'V': //valina
				seqtyp[i] = aVAL;
				break;
			case 'F': //fenilalanina
				seqtyp[i] = aPHE;
				break;
			case 'P': //prolina
				seqtyp[i] = aPRO;
				break;
			case 'L': //leucina
				seqtyp[i] = aLEU;
				break;
			case 'I': //iso leucina
				seqtyp[i] = aILE;
				break;
			case 'R': //argenina
				seqtyp[i] = aARG;
				break;
			case 'D': //acido aspartico
				seqtyp[i] = aASP;
				break;
			case 'E': //acido glutamico
				seqtyp[i] = aGLU;
				break;
			case 'S': //serina
				seqtyp[i] = aSER;
				break;
			case 'T': //treonina
				seqtyp[i] = aTHR;
				break;
			case 'C': //cisteina
				seqtyp[i] = aCYS;
				break;
			case 'N': //asparagina
				seqtyp[i] = aASN;
				break;
			case 'Q': //glutanima
				seqtyp[i] = aGLN;
				break;
			case 'H': //histidina
				seqtyp[i] = aHIS;
				break;
			case 'K': //lisina
				seqtyp[i] = aLYS;
				break;
			case 'Y': //tirosina
				seqtyp[i] = aTYR;
				break;
			case 'M': //metionina
				seqtyp[i] = aMET;
				break;
			case 'W': //triptofano
				seqtyp[i] = aTRP;
				break;
			case 'G': //glicina
				seqtyp[i] = aGLY;
				break;
		}

	}
}*/
