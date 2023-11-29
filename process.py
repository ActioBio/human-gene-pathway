import csv
import pandas as pd
import gzip
from typing import Dict
from pathlib import Path

GENE_DATA_FILE = "data/input/protein_coding_gene.csv"
PATHWAYCOMMONS_FILE = "data/input/PathwayCommons12.All.hgnc.gmt.gz"
WIKIPATHWAYS_FILE = "data/input/wikipathways-Homo_sapiens.gmt"
NODE_OUTPUT_FILE = "data/output/csv/node_Pathway.csv"
EDGE_OUTPUT_FILE = "data/output/csv/edge_Gene_participatesIn_Pathway.csv"

def load_gene_data_from_csv(file_path):
    gene_ids_set = set()
    gene_symbol_to_id_map = {}

    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            gene_id = row['GeneID']
            gene_symbol = row['Symbol']
            gene_ids_set.add(gene_id)
            gene_symbol_to_id_map[gene_symbol] = gene_id

    return gene_ids_set, gene_symbol_to_id_map

def read_gmt(path: str):
    open_func = gzip.open if path.endswith('.gz') else open

    with open_func(path, 'rt') as read_file:
        return [(row[0], row[1], set(row[2:])) for row in csv.reader(read_file, delimiter='\t') if row]

def parse_description(description: str):
    return dict(item.split(': ', 1) for item in description.split('; ') if ': ' in item)

def process_pathwaycommons_data(path: str, gene_symbol_to_id_map: Dict[str, int]) -> pd.DataFrame:
    data = []

    for url, description, genes in read_gmt(path):
        filtered_genes = sorted({gene_symbol_to_id_map.get(gene) for gene in genes} - {None})

        if any(gene in gene_symbol_to_id_map for gene in genes):
            data.append({
                'identifier': Path(url).stem,
                'name': parse_description(description).get('name'),
                'url': url,
                'genes': filtered_genes
            })

    return pd.DataFrame(data)

def process_wikipathways_data(wikipathways_data, gene_ids_set: set):
    wikipath_df = pd.DataFrame(wikipathways_data, columns=['name', 'description', 'genes'])
    wikipath_df['name'] = wikipath_df['name'].str.split('%').str[0]
    wikipath_df['genes'] = wikipath_df['genes'].apply(lambda genes: sorted({int(gene) for gene in genes} & gene_ids_set))
    
    return wikipath_df[wikipath_df['genes'].astype(bool)]

def create_combined_df(pc_df, wikipath_df):
    wikipath_df = wikipath_df.assign(
        identifier=wikipath_df['description'].str.rsplit('/', n=1).str[-1],
        url=wikipath_df['description']
    )

    combined_df = pd.concat([wikipath_df, pc_df], ignore_index=True).drop_duplicates('genes')
    
    return combined_df[['identifier', 'name', 'url', 'genes']].sort_values('identifier')

def write_to_node_csv(dataframe: pd.DataFrame, filepath: str) -> None:
    dataframe[['identifier', 'name', 'url']].to_csv(filepath, index=False)

def write_to_edge_csv(dataframe: pd.DataFrame, filepath: str) -> None:
    edge_data = []

    for _, row in dataframe.iterrows():
        pathway_id = row['identifier']
        for gene_id in row['genes']:
            edge_data.append({'source_id': gene_id, 'target_id': pathway_id})

    edges_df = pd.DataFrame(edge_data)
    edges_df.to_csv(filepath, index=False)

def main():
    gene_ids_set, gene_symbol_to_id_map = load_gene_data_from_csv(GENE_DATA_FILE)

    pathwaycommons_df = process_pathwaycommons_data(PATHWAYCOMMONS_FILE, gene_symbol_to_id_map)
    wikipathways_data = read_gmt(WIKIPATHWAYS_FILE)
    wikipath_df = process_wikipathways_data(wikipathways_data, gene_ids_set)

    combined_pathway_df = create_combined_df(pathwaycommons_df, wikipath_df)
    
    write_to_node_csv(combined_pathway_df, NODE_OUTPUT_FILE)
    write_to_edge_csv(combined_pathway_df, EDGE_OUTPUT_FILE)

if __name__ == "__main__":
    main()
