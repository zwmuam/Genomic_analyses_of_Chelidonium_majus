"""
This script extracts from OrthoFinder results gene trees and corresponding sequence files
for genes of interest provided in a fasta file.
"""
from pathlib import Path
from shutil import copyfile
from typing import Set

from Bio import Phylo
from tqdm import tqdm


def get_fasta_ids(fasta_path: Path):
    """
    Simple fasta parser to get sequence names from a fasta file.
    :param fasta_path: Path to the fasta file.
    :return: A set of sequence seq_ids.
    """
    with fasta_path.open('r') as fasta:
        seq_ids = {line.strip().lstrip('>').split(' ')[0] for line in fasta if line.startswith('>')}
    return seq_ids


def get_trees_from_gene_set(tree_directory: Path,
                            sequence_directory: Path,
                            gene_set: Set[str]) -> Set[Path]:
    """
    Take a list of gene ids and return the tree and sequence files that contain those ids.
    :param tree_directory: The directory containing the trees.
    :param sequence_directory: The directory containing the sequences.
    :param gene_set: A list of gene ids.
    :return: A dictionary with the gene id as key and the tree file as value.
    """
    trees_of_interest = set()
    found_genes = set()
    # select only the tree files
    tree_files = [f for f in tree_directory.iterdir() if f.suffix in {'.nwk', '.txt'}]
    with tqdm(total=len(tree_files)) as bar:
        bar.set_description('searching')
        for tree_file in tree_files:
            # parse the tree file
            tree = Phylo.read(tree_file, 'newick')
            leaves_of_interest = []
            for leaf in tree.get_terminals():
                if leaf.name in gene_set:  # check if the leaf is in the gene set
                    leaves_of_interest.append(leaf.name)
                    found_genes.add(leaf.name)
            if leaves_of_interest:
                bar.set_description(f'{len(found_genes)} found')
                trees_of_interest.add(tree_file)
            bar.update()

    missing_genes = gene_set.difference(found_genes)
    if not trees_of_interest:
        raise ValueError(f'No trees found in {tree_directory.as_posix()} containing any of the provided genes')
    else:
        print(f'Found {len(trees_of_interest)} trees containing {len(found_genes)} genes')
        print(f'Missing genes: {missing_genes} (probably singletons)')

    # find corresponding sequence files
    sequences_of_interest = set()
    for tree in trees_of_interest:
        og_number = tree.stem.split('_')[0]
        sequence_file = sequence_directory.joinpath(f'{og_number}.fa')
        if not sequence_file.exists():
            raise FileNotFoundError(f'{sequence_file.as_posix()} not found')
        else:
            sequences_of_interest.add(sequence_file)

    # find sequence files for singleton genes
    seq_files = [f for f in sequence_directory.iterdir() if f.suffix in {'.fa', }]
    with tqdm(total=len(seq_files)) as bar:
        bar.set_description('searching')
        for seq_file in seq_files:
            seq_genes = get_fasta_ids(seq_file)
            seq_genes = set([f'Chelidonium_majus_{g}' for g in seq_genes])
            found_genes = set()
            for gene in missing_genes:
                if gene in seq_genes:
                    sequences_of_interest.add(seq_file)
                    found_genes.add(gene)
            missing_genes = missing_genes.difference(found_genes)
            bar.set_description(f'{len(missing_genes)} remain')
            if not missing_genes:
                print('All missing genes found among singletons')
                break
            bar.update()
    if missing_genes:
        raise ValueError(f'Missing sequences for: {missing_genes}')
    return trees_of_interest | sequences_of_interest


if __name__ == '__main__':
    # INPUTS
    # set the directory containing the trees
    tree_dir = Path('./OrthoFinder_Results/Resolved_Gene_Trees')
    og_fasta_directory = Path('./OrthoFinder_Results/Orthogroup_Sequences')
    fasta_dir = Path('./C_majus_JBPSLD000000000_families_of_interest')

    # OUTPUT
    out_directory = Path('./Orthogroups_of_interest')
    out_directory.mkdir(exist_ok=True, parents=True)

    for fasta_p in [f for f in fasta_dir.iterdir() if f.suffix in {'.fa', '.faa', '.fasta'}]:
        genes_of_interest = get_fasta_ids(fasta_p)
        files2extract = get_trees_from_gene_set(tree_dir,
                                                og_fasta_directory,
                                                genes_of_interest)
        out_subdir = out_directory.joinpath(fasta_p.stem)
        out_subdir.mkdir(exist_ok=True, parents=True)
        with tqdm(total=len(files2extract)) as bar:
            for efile in files2extract:
                copyfile(efile, out_subdir.joinpath(efile.name))
                bar.set_description(f'{efile.name} copied')
                bar.update()
