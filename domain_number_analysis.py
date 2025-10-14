"""
This script analyzes the number of protein domains per orthogroup and per organism
input is gff with protein annotation generated for OrthoFinder orthologous groups
by InterProScan plugin in Geneious Prime
"""
from collections import defaultdict, Counter
from pathlib import Path
from typing import Dict

import pandas as pd
import plotly.express as px
from loguru import logger


def seq_to_organism_dict(sequence_ids_txt: Path,
                         species_ids_txt: Path) -> Dict[str, str]:
    """
    Create a dictionary mapping sequence IDs to organism names
    based ID files located in WorkingDirectory inside OrthoFinder results
    :param sequence_ids_txt: path to the text file with sequence IDs (SpeciesIDs.txt)
    :param species_ids_txt: path to the text file with species names (SequenceIDs.txt)
    :return: {sequence_id: organism_name, ...}
    """

    species_dict = {}
    with species_ids_txt.open() as sp:
        for line in sp:  # example line: "0: Arabidopsis thaliana.faa"
            species_number, species_name = line.strip().split(': ')
            species_name = species_name.replace('.faa', '')
            species_dict[species_number] = species_name

    out_dict = {}
    with sequence_ids_txt.open() as sq:
        for line in sq:  # example line: "0_0: NP_001030613.1 hypothetical protein 1 [Arabidopsis thaliana]"
            seq_n, seq_id, *rest = line.strip().split(' ')
            species_number = seq_n.split('_')[0]
            out_dict[seq_id] = species_dict[species_number]

    return out_dict


def parse_gff(gff_file: Path,
              organism_dict: Dict[str, str],
              out_dir: Path = None,
              id_filter: str = None,
              database_filter: str = None) -> pd.DataFrame:
    """
    Use Geneious gff with protein domain annotations organism mapping from OrthoFinder
    to create a dataframe with domains annotated with organism and orthologous group information
    :param gff_file: path to the gff file with protein domain annotations
    :param organism_dict: dictionary mapping sequence IDs to organism names
    :param out_dir: path to save the output table
    :param id_filter: if provided, only domains with this ID will be included
    :param database_filter: if provided, only domains from this database will be included
    :return: table with:
                - protein_ID
                - domain_ID
                - domain_name
                - domain_length
                - orthologous_group
                - organism_name
    """
    domain_records = []
    with gff_file.open() as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            else:
                sequence_name, source, feature, start, end, score, strand, frame, attr = line.strip().split('\t')
                try:
                    sequence_id, orthologous_group = sequence_name.split('|')
                except ValueError:
                    sequence_id, orthologous_group = sequence_name, '_SINGLETONS_'
                length = int(end) - int(start)
                attr = dict([a.split('=') for a in attr.split(';')])
                domain_id = attr['InterPro IdX'].split('>')[1].split('<')[0] if 'InterPro IdX' in attr else '-'
                database = attr['Database']
                if id_filter and domain_id != id_filter or database_filter and database != database_filter:
                    continue  # skip this domain if it does not match the filter
                domain_name = attr['InterPro Name'] if 'InterPro Name' in attr else '-'
                domain_dict = {'protein_ID': sequence_id,
                               'domain_length': length,
                               'orthologous_group': orthologous_group,
                               'domain_ID': domain_id,
                               'domain_name': domain_name,
                               'organism': organism_dict[sequence_id]}
                domain_records.append(domain_dict)

    domain_table = pd.DataFrame.from_records(domain_records)
    export_path = out_dir.joinpath('domain_table.xlsx')

    if out_dir is not None:
        domain_table.to_excel(export_path, index=False)

    return domain_table


def domain_number_analysis(domain_table: pd.DataFrame,
                           category: str,
                           out_dir: Path):
    """
    Analyze the number of domains per protein in each category (either organism or orthologous group)
    export the table and create a stacked bar chart
    :param domain_table: dataframe with domains from parse_gff function
    :param category: either 'organism' or 'orthologous_group'
    :param out_dir: path to save both table and plot
    """
    in_records = domain_table.to_dict('records')
    domains_per_protein = defaultdict(lambda: 0)
    for r in in_records:
        domains_per_protein[(r['protein_ID'], r[category])] += 1
    to_count = []
    for (seq, category_value), n_domains in domains_per_protein.items():
        to_count.append((category_value, n_domains))
    counted = Counter(to_count).items()
    domain_count_table = pd.DataFrame([{category: h, 'n_domains': i, 'n_proteins': c} for (h, i), c in counted])
    human_readable_dcounts = domain_count_table.pivot(index=category, columns='n_domains', values='n_proteins')
    human_readable_dcounts.fillna(0, inplace=True)
    table_path = out_dir.joinpath(f'{category}_domain_count_summary.xlsx')
    human_readable_dcounts.to_excel(table_path)
    logger.info(f'Written to {table_path}')

    # sort the dataframe by domain count
    domain_count_table.sort_values(by='n_domains', inplace=True)
    domain_count_table['n_domains'] = domain_count_table['n_domains'].astype(str)
    fig = px.bar(domain_count_table, x=category, y='n_proteins', color='n_domains', title='Domain count',
                 color_discrete_map={'1': 'green',
                                     '2': 'blue',
                                     '3': 'yellow',
                                     '4': 'purple',
                                     '7': 'red'}, )  # this color needs to be adjusted for specific data
    chart_path = out_dir.joinpath(f'{category}_domain_stacked_bar_chart.pdf')
    fig.write_image(chart_path)
    logger.info(f'Written to {chart_path}')


if __name__ == '__main__':
    # INPUT
    protein_gff = Path('./Domains_from_MLP_ogs.gff')
    orthofinder_working_dir = Path('./OrthoFinder_Results/WorkingDirectory')
    orthofinder_species_ids = orthofinder_working_dir.joinpath('SpeciesIDs.txt')
    orthofinder_sequence_ids = orthofinder_working_dir.joinpath('SequenceIDs.txt')
    filtered_id = 'IPR000916'
    filtered_database = 'PFAM'  # IPR000916 id is present in several databases, this filter avoids double counting

    # OUTPUT
    output_dir = Path('./figures_&_tables')
    output_dir.mkdir(exist_ok=True, parents=True)

    # EXECUTION
    sp_dict = seq_to_organism_dict(sequence_ids_txt=orthofinder_sequence_ids,
                                        species_ids_txt=orthofinder_species_ids)

    # parse gff from InterProScan Geneious plugin into a table
    domain_df = parse_gff(gff_file=protein_gff,
                          organism_dict=sp_dict,
                          out_dir=output_dir,
                          id_filter=filtered_id,
                          database_filter=filtered_database)

    # analyse domain counts per orthologous group
    domain_number_analysis(domain_df,
                           category='orthologous_group',
                           out_dir=output_dir)

    # analyse domain counts per organism
    domain_number_analysis(domain_df,
                           category='organism',
                           out_dir=output_dir)
