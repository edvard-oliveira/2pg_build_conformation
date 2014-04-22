#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>

#include "gromacs.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"
#include "consts.h"
#include "string_owner.h"
#include "osutil.h"
#include "parameters_type.h"
#include "pdbio.h"


#define TAM_LINE_ENER 50
#define MAX_ENERGY 9999999999999999999999999.9999
#define MIN_ENERGY -9999999999999999999999999.9999
#define MAX_VALUE 40
#define MAX_VALUE_G_SAS 3

static char *program = NULL; /* program (executable) file name */
static char *filenm1 = NULL; /* data file names */
static char *filenm2 = NULL;
static char *filenm3 = NULL;
static char *filenm4 = NULL;
static char *filenm5 = NULL;
//Options for GROMACS programs
static char *f_step0 = NULL;
static char *prot_gro = NULL;
static char *prot_top = NULL;
static char *prot_tpr = NULL;
static char *prot_trr = NULL;
static char *prot_log = NULL;
static char *confout_gro = NULL;
static char *posre_itp = NULL;
static char *mdout_mdp = NULL;
static char *file_energy_computed_ener_edr = NULL;
static char *prot_sys_trr = NULL;
static char *prot_sys_tpr = NULL;
static char *prot_sys_top = NULL;
static char *prot_sys_gro = NULL;
static char *energy_xvg = NULL;
static char *traj_xtc = NULL;
static char *opt_f = NULL;
static char *opt_s = NULL;
static char *opt_o = NULL;
static char *opt_c = NULL;
static char *opt_ff = NULL;
static char *opt_water = NULL;
static char *opt_none = NULL;
static char *opt_p = NULL;
static char *opt_ignh = NULL;
static char *xvg_1 = NULL;
static char *opt_rerun = NULL;
static char *opt_e = NULL;
static char *opt_g = NULL;
static char *opt_num = NULL;
/* Stores the values of g_sas program 
* g_sas_values[0] Hydrophobic
* g_sas_values[1] Hydrophilic
* g_sas_values[3] Total Area
*/
static double *g_sas_values = NULL;


//It is based on GROMACS version 4.5.3
static option_fitness_gromacs_t option_g_energy_program [] = {
      		                                                  {gmx_potential_ener, "12","Potential"},
      		                                                  {gmx_edw_ener,"8","LJ-14"},
      		                                                  {gmx_elel_ener,"9","Coulomb-14"},
      		                                                  {gmx_hydrophobic,"-1","Hydrophobic"},
      		                                                  {gmx_hydrophilic,"-1","Hydrophilic"},
      		                                                  {gmx_total_area,"-1","Total_Area"},
      		                                                  {gmx_gyrate,"-1","Gyrate"},
      		                                                  {gmx_hbond,"-1","H_Bond"},
      		                                                  {gmx_hbond_main,"-1","H_Bond_Main"},
      		                                                  {gmx_GBSA_Solvatation,"-1","GBSA_Sol"},
      		                                                  {gmx_GB_Polarization,"6","GB-Polarization"},
      		                                                  {gmx_Nonpolar_Sol,"7","Nonpolar-Sol."},
      		                                                  {gmx_stride_total,"-1","Stride_total"},
      		                                                  {gmx_stride_helix,"-1","Stride_helix"},
      		                                                  {gmx_stride_beta,"-1","Stride_beta"}
                                                             };


static inline int run_program(const char *file, char *const argv[])
{
	pid_t pid;
	int status;
	int out, err; /* file descriptors for stdout and stderr */

	pid = fork();

	if (pid == -1) {
		perror("Fork failed to create process");
		return 0;
	} else if (pid == 0) {
		/* child process */
		out = open("/dev/null", O_RDONLY);
		err = open("/dev/null", O_RDONLY);
		dup2(out, 1);
		dup2(err, 2);

		execv(file, argv);
		perror("Execv failed to run program (run_program)");
		_exit(EXIT_FAILURE);
	} else {
		/* parent process */
		waitpid(pid, &status, 0);
	}

	return 1;
}

/* Start a program as in:
 * echo "pipe_msg" | program [args...] */
static inline int run_program_after_pipe(const char *pipe_msg, const char *file,
							char *const argv[])
{
	pid_t pid;
	int status;
	int fd[2], out, err;

	if (pipe(fd) == -1) {
		perror("Failed to create pipe");
		return 0;
	}

	pid = fork();

	if (pid == -1) {
		perror("Fork failed to create process");
		return 0;
	} else if (pid == 0) {
		/* child process */
		/* supress output */
		out = open("/dev/null", O_WRONLY);
		err = open("/dev/null", O_WRONLY);
		dup2(out, 1);
		dup2(err, 2);

		/* read from the pipe */
		close(fd[1]);
		dup2(fd[0], 0);
		execv(file, argv);
		perror("Execv failed to run program");
		_exit(EXIT_FAILURE);
	} else {
		/* parent process */
		close(fd[0]);
		write(fd[1], pipe_msg, strlen(pipe_msg) + 1);
		close(fd[1]);
		waitpid(pid, &status, 0);
	}

	return 1;
}

/*
 * Run nprogs programs with pipes interconecting them and (optionally, if
 * output_file is not NULL) write the output to a file, as in:
 * 
 * prog0 [args0] | prog1 [args1] | ... | progn-1 [argsn-1] > output_file
 *
 * where:
 * argv_list[i] = argsi (argsi is in the format expected by execv(3))
 */
static inline int run_programs_with_pipe(int nprogs, char ***const argv_list,
							const char *output_file)
{
	int oldpipe[2], newpipe[2];
	int status;
	int i;
	int ret_value;
	char **argv;

	if (pipe(oldpipe) == -1) {
		perror("Failed to create pipe");
		return 0;
	}

	ret_value = 1;
	/* create first child */
	switch(fork()) {
		case -1:
			perror("Fork failed to create process");
			return 0;

		case 0:
			/* 1st child process writes to oldpipe, reads from none */
			close(oldpipe[0]);
			dup2(oldpipe[1], 1);

			argv = argv_list[0];
			execv(argv[0], argv);
			perror("Execv failed to run program");
			_exit(EXIT_FAILURE);
	}

	/* create intermediate children */
	/* 1st child process already created, last one to be created later */
	for (i = 1; i < nprogs-1; i++) {
		if (pipe(newpipe) == -1) {
			perror("Failed to create pipe");
			ret_value = 0;
			goto end;
		}

		switch (fork()) {
			case -1:
				perror("Fork failed to create process");
				ret_value = 0;
				goto end;

			case 0:
				/* child reads from oldpipe, writes to newpipe */
				close(oldpipe[1]);
				close(newpipe[0]);
				dup2(oldpipe[0], 0);
				dup2(newpipe[1], 1);

				argv = argv_list[i];
				execv(argv[0], argv);
				perror("Execv failed to run program");
				_exit(EXIT_FAILURE);

			default:
				/* parent process */
				close(oldpipe[0]);
				close(oldpipe[1]);
				oldpipe[0] = newpipe[0];
				oldpipe[1] = newpipe[1];
		}
	}

	/* create last child */
	switch(fork()) {
		case -1:
			perror("Fork failed to create process");
			ret_value = 0;
			goto end;

		case 0:
			/* last child process reads from oldpipe, maybe writes to file */
			close(oldpipe[1]);
			dup2(oldpipe[0], 0);

			if (output_file != NULL) {
				int outfd;

				outfd = open(output_file, O_WRONLY | O_CREAT | O_TRUNC,
						S_IRUSR | S_IWUSR | S_IRGRP);
				if (outfd == -1) {
					perror("Error opening output file");
					_exit(EXIT_FAILURE);
				}
				dup2(outfd, 1);
			}

			argv = argv_list[nprogs-1];
			execv(argv[0], argv);
			perror("Execv failed to run program");
			_exit(EXIT_FAILURE);
	}

	close(oldpipe[0]);
	close(oldpipe[1]);

	end:
	for (; i >= 0; i--) /* wait for all the created children processes */
		wait(&status);

	return ret_value;
}

void initialize_g_sas_values(){
	//Hydrophobic
	g_sas_values[0] = MAX_ENERGY;
	//Hydrophilic
	g_sas_values[1] = MIN_ENERGY;
	//Total Area 
	g_sas_values[2] = MAX_ENERGY;
}

/** Initialize GROMACS execution
* This function must be called exactly once before any other one in this file
 * to obtain memory for internal variables used to interact with Gromacs.
 *
 * Call finish_gromacs_execution when done using Gromacs
 */
void init_gromacs_execution (){

	program = Malloc(char, MAX_COMMAND);
	filenm1 = Malloc(char, MAX_COMMAND);
	filenm2 = Malloc(char, MAX_COMMAND);
	filenm3 = Malloc(char, MAX_COMMAND);
	filenm4 = Malloc(char, MAX_COMMAND);
	filenm5 = Malloc(char, MAX_COMMAND);

	f_step0 = Malloc(char, MAX_FILE_NAME);
	prot_gro = Malloc(char, MAX_FILE_NAME);
	prot_top = Malloc(char, MAX_FILE_NAME);
	prot_tpr = Malloc(char, MAX_FILE_NAME);
	prot_trr = Malloc(char, MAX_FILE_NAME);
	prot_log = Malloc(char, MAX_FILE_NAME);
	confout_gro = Malloc(char, MAX_FILE_NAME);
	posre_itp = Malloc(char, MAX_FILE_NAME);
	mdout_mdp = Malloc(char, MAX_FILE_NAME);
	file_energy_computed_ener_edr = Malloc(char, MAX_FILE_NAME);
	prot_sys_trr = Malloc(char, MAX_FILE_NAME);
	prot_sys_tpr = Malloc(char, MAX_FILE_NAME);
	prot_sys_top = Malloc(char, MAX_FILE_NAME);
	prot_sys_gro = Malloc(char, MAX_FILE_NAME);
	energy_xvg = Malloc(char, MAX_FILE_NAME);
	traj_xtc = Malloc(char, MAX_FILE_NAME);
	xvg_1  = Malloc(char, MAX_FILE_NAME);
	opt_f = Malloc(char,3);
	opt_s = Malloc(char,3);
	opt_o = Malloc(char,3);
	opt_ff  = Malloc(char, 4);
	opt_water = Malloc(char, 7);
	opt_none = Malloc(char, 6);
	opt_p =  Malloc(char, 3);
	opt_ignh = Malloc(char, 7);
	opt_c = Malloc(char, 3);
	opt_rerun  = Malloc(char, 10);
	opt_e = Malloc(char, 3);
	opt_g = Malloc(char, 3);
	opt_num = Malloc(char, 5);
	g_sas_values = Malloc(double, MAX_VALUE_G_SAS);
	
	strcpy(opt_f, "-f");	
	strcpy(opt_o, "-o");
	strcpy(opt_ff, "-ff");
	strcpy(opt_water, "-water");
	strcpy(opt_none, "none");
	strcpy(opt_p, "-p");
	strcpy(opt_ignh, "-ignh");
	strcpy(opt_c, "-c");
	strcpy(opt_f, "-f");
	strcpy(opt_s, "-s");
	strcpy(opt_o, "-o");
	strcpy(f_step0, "step0*.pdb");
	strcpy(prot_gro, "prot.gro");
	strcpy(prot_top, "prot.top");
	strcpy(prot_tpr, "prot.tpr");
	strcpy(prot_trr, "prot.trr");
	strcpy(prot_log, "prot.log");
	strcpy(confout_gro, "confout.gro");
	strcpy(posre_itp, "posre.itp");
	strcpy(mdout_mdp, "mdout.mdp");
	strcpy(file_energy_computed_ener_edr, "file_energy_computed.ener.edr");
	strcpy(prot_sys_trr, "prot_sys.trr");
	strcpy(prot_sys_tpr, "prot_sys.tpr");
	strcpy(prot_sys_top, "prot_sys.top");
	strcpy(prot_sys_gro, "prot_sys.gro");	
	strcpy(energy_xvg, "energy.xvg");
	strcpy(traj_xtc, "traj.xtc");
	strcpy(xvg_1, "*.xvg_1*");
	strcpy(opt_rerun, "-rerun");
	strcpy(opt_e, "-e");	
	strcpy(opt_g, "-g");
	strcpy(opt_num, "-num");
	initialize_g_sas_values();

}

/** Finish GROMACS execution 
* This function is to be called when the Gromacs funcions below will no longer
* be used in order to release memory used by internal variables 
*/
void finish_gromacs_execution(){

	free(program);
	free(filenm1);
	free(filenm2);
	free(filenm3);
	free(filenm4);
	free(filenm5);
	free(f_step0);
	free(prot_gro);
	free(prot_top);
	free(prot_tpr);
	free(prot_trr);
	free(prot_log);
	free(confout_gro);
	free(posre_itp);
	free(mdout_mdp);
	free(file_energy_computed_ener_edr);
	free(prot_sys_trr);
	free(prot_sys_tpr);
	free(prot_sys_top);
	free(prot_sys_gro);
	free(energy_xvg);
	free(traj_xtc);
	free(opt_f);
	free(opt_s);
	free(opt_o);
	free(opt_ff);
	free(opt_water);
	free(opt_none);
	free(opt_p);
	free(opt_ignh);
	free(opt_c);
	free(xvg_1);
	free(opt_rerun);
	free(opt_e);
	free(opt_g);
	free(opt_num);
	free(g_sas_values);
	program = NULL;
	filenm1 = NULL;
	filenm2 = NULL;
	filenm3 = NULL;
	filenm4 = NULL;
	filenm5 = NULL;
}


option_g_energy_t get_option_g_energy_t_from_type_fitness_energy(const type_fitness_energies_t *fit_ener){
	/*Receives type_fitness_energies_t and returns its option_g_energy_t
	 * For example: fit_ener_potential is gmx_potential_ener
	 */
	if (*fit_ener == fit_ener_potential){
		return gmx_potential_ener;
	}else if (*fit_ener == fit_ener_edw){
		return gmx_edw_ener;
	}else if (*fit_ener == fit_ener_ele){
		return gmx_elel_ener;
	}else if (*fit_ener == fit_hydrophobic){
		return gmx_hydrophobic;
	}else if (*fit_ener == fit_hydrophilic){
		return gmx_hydrophilic;
	}else if (*fit_ener == fit_total_area){
		return gmx_total_area;
	}else if (*fit_ener == fit_gyrate){
		return gmx_gyrate;
	}else if (*fit_ener == fit_hbond){
		return gmx_hbond;
	}else if (*fit_ener == fit_hbond_main){
		return gmx_hbond_main;
	}else if (*fit_ener == fit_GBSA_Solvatation){
		return gmx_GBSA_Solvatation;
	}else if (*fit_ener == fit_stride_total){
		return gmx_stride_total;
	}else if (*fit_ener == fit_stride_helix){
		return gmx_stride_helix;
	}else if (*fit_ener == fit_stride_beta){
		return gmx_stride_beta;
	}else{
		fatal_error("Option did not find at get_option_g_energy_t_from_type_fitness_energy function. Please check it! ");
	}
}

/** Calls toremove files used in simulation.
* Simulation means an execution to calculate the objectives by GROMACS.
*/
static void clean_gromacs_simulation(const char *path_local_execute){
	delete_file(path_local_execute, f_step0);
	delete_file(path_local_execute, prot_gro);
	delete_file(path_local_execute, prot_top);
	delete_file(path_local_execute, prot_tpr);
	delete_file(path_local_execute, prot_trr);
	delete_file(path_local_execute, prot_log);	
	delete_file(path_local_execute, confout_gro);
	delete_file(path_local_execute, posre_itp);
	delete_file(path_local_execute, mdout_mdp);
	delete_file(path_local_execute, file_energy_computed_ener_edr);
	delete_file(path_local_execute, prot_sys_trr);
	delete_file(path_local_execute, prot_sys_tpr);
	delete_file(path_local_execute, prot_sys_top);
	delete_file(path_local_execute, prot_sys_gro);
	delete_file(path_local_execute, energy_xvg);
	delete_file(path_local_execute, traj_xtc);
	delete_file(path_local_execute, xvg_1);

}

/** Obtaing the value of objective in Double representation
* This function was created to be a pattern for getting
* the value of objective when it will be in Double
* instead of char.
*/
double get_objective_value(const char *value){
	double aux = MAX_ENERGY;
	aux = str2double(value);
	if ( (aux ==  -1.000000) || (isnan(aux) == btrue)){
		aux = MAX_ENERGY;
	}
	return aux;
}

static void build_pdb_file_name(char *pdb_file_name, const char *aux_name,
		const char *__restrict prefix){
	strcpy(pdb_file_name, prefix);
	strcat(pdb_file_name,aux_name);
	strcat(pdb_file_name,".pdb");
}

/** Creates a tpr file
*/
void build_tpr_file(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs, const char *force_field, const char *mdp_file){
	
	char *pdbfile_aux;
	char *force_field_aux;
	char *mdp_file_aux;

	char *pdb2gmx_args[13];
	char *grompp_args[11];

	pdbfile_aux = Malloc(char, MAX_FILE_NAME);
	force_field_aux = Malloc(char, MAX_FORCE_FIELD_NAME);
	mdp_file_aux = Malloc(char, MAX_FILE_NAME);

	/* pdb2gmx */
	strcpy(program, path_gromacs_programs);
	strcat(program, "pdb2gmx");
	pdb2gmx_args[0] = program;
	//pdb
	strcpy(pdbfile_aux, pdbfile);
	pdb2gmx_args[1] = opt_f;
	strcpy(filenm1, local_execute);	
	strcat(filenm1, pdbfile_aux);
	pdb2gmx_args[2] = filenm1;
	//gro
	pdb2gmx_args[3] = opt_o;
	strcpy(filenm2, local_execute);
	strcat(filenm2, "prot.gro");
	pdb2gmx_args[4] = filenm2;
	//top
	pdb2gmx_args[5] = opt_p;
	strcpy(filenm3, local_execute);
	strcat(filenm3, "prot.top");
	pdb2gmx_args[6] = filenm3;
	//force field
	pdb2gmx_args[7] = opt_ff;
	strcpy(force_field_aux, force_field);
	pdb2gmx_args[8] = force_field_aux;
	//water
	pdb2gmx_args[9] = opt_water;
	pdb2gmx_args[10] = opt_none;
	//Hydrogen
	pdb2gmx_args[11] = opt_ignh;

	pdb2gmx_args[12] = NULL;	

	if (!run_program(program, pdb2gmx_args)){
		fatal_error("Failed to run pdb2gmx at build_tpr_file function \n");
	}

	/* grompp */
	strcpy(program, path_gromacs_programs);
	strcat(program, "grompp");
	grompp_args[0] = program;
	//mdp
	grompp_args[1] = opt_f;
	strcpy(mdp_file_aux, mdp_file);
	strcpy(filenm1, local_execute);
	strcat(filenm1, mdp_file_aux);
	grompp_args[2] = filenm1;
	//top
	grompp_args[3] = opt_p;
	strcpy(filenm2, local_execute);
	strcat(filenm2, "prot.top");
	grompp_args[4] = filenm2;
	//tpŕ
	grompp_args[5] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, "prot.tpr");
	grompp_args[6] = filenm3;
	//gro
	grompp_args[7] = opt_c;
	strcpy(filenm4, local_execute);
	strcat(filenm4, "prot.gro");
	grompp_args[8] = filenm4;

	grompp_args[9] = NULL;
	grompp_args[10] = NULL;

	if (!run_program(program, grompp_args))
		fatal_error("Failed to run grompp at build_tpr_file function \n");

	free(pdbfile_aux);
	free(force_field_aux);
	free(mdp_file_aux);
}

/** Calls mdrun program for miminization
*/
void call_mdrun2minimization(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs){
	char *mdrun_args[13];
	/* mdrun */
	strcpy(program, path_gromacs_programs);
	strcat(program, "mdrun");
	mdrun_args[0] = program;
	//tpŕ
	mdrun_args[1] = opt_s;		
	strcpy(filenm1, local_execute);
	strcat(filenm1, prot_tpr);
	mdrun_args[2] = filenm1;
	//trr
	mdrun_args[3] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, prot_trr);
	mdrun_args[4] = filenm3;
	//energy file
	mdrun_args[5] = opt_e;
	strcpy(filenm4, local_execute);
	strcat(filenm4, file_energy_computed_ener_edr);
	mdrun_args[6] = filenm4;
	//log file
	mdrun_args[7] = opt_g;
	strcpy(filenm5, local_execute);
	strcat(filenm5, prot_log);
	mdrun_args[8] = filenm5;
	//PDB - Last Frame
	mdrun_args[9] = opt_c;
	strcpy(filenm5, local_execute);
	strcat(filenm5, pdbfile);
	mdrun_args[10] = filenm5;

	mdrun_args[11] = NULL;
	mdrun_args[12] = NULL;

	if (!run_program(program, mdrun_args))
		fatal_error("Failed to run mdrun at call_mdrun2minimization function \n");
}



/** Applies minimization process by GROMACS
*/
void minimization_gromacs(pdb_atom_t *pdb_atoms, char *pdbfile_ret, int *numatom_after_min, 
	const input_parameters_t *in_para, const int *numatom){
	char *pdbfile;
	
	pdbfile = Malloc(char, MAX_FILE_NAME);
	strcpy(pdbfile, "conformation_minimization.pdb");
	

	save_pdb_file(in_para->path_local_execute, pdbfile, numatom,
			pdb_atoms, NULL);

    //Building Generic tpr file. It will be used in all GROMACS execution   	    
   	build_tpr_file(pdbfile, in_para->path_local_execute, in_para->path_gromacs_programs, 
	        	in_para->force_field, in_para->mdp_file_min);

   	if (in_para->gromacs_energy_min != ener_min_none){
		//Call mdrun program
		call_mdrun2minimization(pdbfile, in_para->path_local_execute, 
	     		 in_para->path_gromacs_programs);

   	}

	//Call to clean the simulation
	clean_gromacs_simulation(in_para->path_local_execute);

	//Getting information 
	char *path_file = path_join_file(in_para->path_local_execute,pdbfile);
	*numatom_after_min = get_num_atom(path_file);
	strcpy(pdbfile_ret, pdbfile);

	free(path_file);
	free(pdbfile);

}