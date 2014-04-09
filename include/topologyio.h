#ifndef OLD_TOPOLOGYIO_H
#define OLD_TOPOLOGYIO_H

#include "topology_types.h"
#include "protein.h"

void _save_topology_file(const char *path, const char *file_name, const top_global_t *top);
amino_t *_load_protein2(const char *file_name_protein, int *n_residues,
		int *numatom, int *num_bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles_type);
void _show_protein_backbone(const protein_backbone_t *prot_back,
		const top_global_t *top);
amino_t *_load_fasta_pdb(const char *file_name_protein, int *n_residues,
		int *numatom, int *num_bond_angles, int *num_protein_side_chains,
		int *num_dihedral_angles_type, boolean_t *has_his);
void _create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein);

static boolean_t _check_pdb_fasta_file(char *line);
#endif
