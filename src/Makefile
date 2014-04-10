#
# --------------------------------------------------
# Author - Rodrigo Antonio Faccioli
# Date   - 22/06/2010
#
# --------------------------------------------------

CXX=g++
OPT=  -c -g
OPT_GSL = -static -Wall -lgsl -lgslcblas -lm -pedantic 
#VALGRIND USE LINE BELOW
#OPT= -c -g -O0

 all:	protpred-Gromacs_pop_initial 
 
protpred-Gromacs_pop_initial:protpred_pop_initial.o LoadConfig.o messages.o protein.o futil.o build_initial_population.o topology.o topologyio.o topologylib.o functions.o populationio.o build_protein_lib.o randomlib.o nerf.o vector_math.o string_owner.o z_matrix.o z_matrix_lib.o z_matrix_io.o pdbatom.o math_owner.o pdbio.o build_protein.o nerfio.o osutil.o random_number_gsl.o load_parameters.o gromacs.o
	$(CXX) -o protpred-Gromacs_pop_initial protpred_pop_initial.o LoadConfig.o messages.o protein.o futil.o build_initial_population.o topology.o topologyio.o topologylib.o functions.o populationio.o build_protein_lib.o randomlib.o nerf.o vector_math.o string_owner.o z_matrix.o z_matrix_lib.o z_matrix_io.o pdbatom.o math_owner.o pdbio.o build_protein.o nerfio.o osutil.o random_number_gsl.o load_parameters.o gromacs.o $(OPT_GSL) 
	

protpred_pop_initial.o:	protpred_pop_initial.c
	$(CXX) $(OPT) protpred_pop_initial.c -I../include

build_initial_population.o:	build_initial_population.c
	$(CXX) $(OPT) build_initial_population.c -I../include

build_protein.o:	build_protein.c
	$(CXX) $(OPT) build_protein.c -I../include
	
build_protein_lib.o:	build_protein_lib.c
	$(CXX) $(OPT) build_protein_lib.c -I../include

populationio.o:	populationio.c
	$(CXX) $(OPT) populationio.c -I../include

LoadConfig.o:	LoadConfig.cpp
	$(CXX) $(OPT) LoadConfig.cpp -I../include

functions.o:	functions.c
	$(CXX) $(OPT) functions.c -I../include
	
topology.o:	topology.c
	$(CXX) $(OPT) topology.c -I../include

topologyio.o:	topologyio.c
	$(CXX) $(OPT) topologyio.c -I../include

topologylib.o:	topologylib.c
	$(CXX) $(OPT) topologylib.c -I../include

pdbio.o:	pdbio.c
	$(CXX) $(OPT) pdbio.c -I../include

protein.o:	protein.c
	$(CXX) $(OPT) protein.c -I../include
	
futil.o:	futil.c
	$(CXX) $(OPT) futil.c -I../include

pdbatom.o:	pdbatom.c
	$(CXX) $(OPT) pdbatom.c -I../include

messages.o:	messages.c
	$(CXX) $(OPT) messages.c -I../include

randomlib.o:	randomlib.c
	$(CXX) $(OPT) randomlib.c -I../include

vector_math.o:	vector_math.c
	$(CXX) $(OPT) vector_math.c -I../include

nerf.o:	nerf.c
	$(CXX) $(OPT) nerf.c -I../include

nerfio.o:	nerfio.c
	$(CXX) $(OPT) nerfio.c -I../include

string_owner.o:	string_owner.c
	$(CXX) $(OPT) string_owner.c -I../include

z_matrix.o:	z_matrix.c
	$(CXX) $(OPT) z_matrix.c -I../include

z_matrix_lib.o:	z_matrix_lib.c
	$(CXX) $(OPT) z_matrix_lib.c -I../include

z_matrix_io.o:	z_matrix_io.c
	$(CXX) $(OPT) z_matrix_io.c -I../include

math_owner.o:	math_owner.c
	$(CXX) $(OPT) math_owner.c -I../include
	
random_number_gsl.o:	random_number_gsl.c
	$(CXX) $(OPT) random_number_gsl.c -I../include

osutil.o:	osutil.c
	$(CXX) $(OPT) osutil.c -I../include

load_parameters.o:	load_parameters.c
	$(CXX) $(OPT) load_parameters.c -I../include
	
fitnessio.o:	fitnessio.c
	$(CXX) $(OPT) fitnessio.c -I../include

gromacs.o:	gromacs.c
	$(CXX) $(OPT) gromacs.c -I../include