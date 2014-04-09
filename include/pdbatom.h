#ifndef OLD_PDBATOM_H
#define OLD_PDBATOM_H

#include "vector_types.h"
#include "enums.h"

#define NUM_ATOM_NAME 5
#define NUM_RES_NAME 4

typedef struct spdbatom{
	type_atoms_t atomid;
    char atmname[NUM_ATOM_NAME];
    char resname[NUM_RES_NAME];
    int resnum;
    own_vector_t coord;
    int atmnumber;
 }pdb_atom_t;

 typedef struct //pdbseqres
  {
    char chain[1];
    char resname[3];
  }pdb_seqres_t;

pdb_atom_t** allocate_Population_pdb(const int *inPopSize, const int *numatom);
void desAllocate_Population_pdb(pdb_atom_t** pdbatoms, const int *inPopSize);
pdb_atom_t * allocate_pdbatom(const int *numatom);
void desAllocate_pdbatom(pdb_atom_t *pdbatoms);
void set_pdb_atom_coordinates(pdb_atom_t *pdbatom,	char *atmname, char *resname,
		const char *chain_name, const int *resnum, 	const float *x,
		const float *y, const float *z, const int *index);
void set_pdb_atom_generic_information(pdb_atom_t *pdbatom,
		char *atmname, char *resname,
		const char *chain_name, const int *resnum,  const int *index);
void get_atom_name_from_atomid(char *atomname, const type_atoms_t *atomid);
type_atoms_t get_atomid_from_atom_name(const  char *__atmname);
const pdb_atom_t * search_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
const pdb_atom_t * search_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname,	const int *num_atom);
const pdb_atom_t * get_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
const pdb_atom_t * get_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname,	const int *num_atom);
int get_number_res_from_atom(const pdb_atom_t *atoms, const int *num_atom);
void get_res_name_from_res_num(char *res_name, const int *num_res,
		const pdb_atom_t *atoms, const int *num_atom);
void renumerate_residue_number(pdb_atom_t *atoms, const int *num_atom);
void rename_atom(pdb_atom_t *atoms, const char *name, const char *name_new,
		const int *res_num, const type_atoms_t *atomid,
		const type_atoms_t *atomid_new, const int *num_atom);
static boolean_t is_residue_number_ok(pdb_atom_t *atoms);
static pdb_atom_t * search_pdb_atom_from_resnum_atomid_alow_change(pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
#endif
