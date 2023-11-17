import csv
import pandas as pd
import gzip
from typing import Dict
from pathlib import Path

def read_gmt(path: str):
    open_func = gzip.open if path.endswith('.gz') else open

    with open_func(path, 'rt') as read_file:
        return [(row[0], row[1], set(row[2:])) for row in csv.reader(read_file, delimiter='\t') if row]

def parse_description(description: str):
    return dict(item.split(': ', 1) for item in description.split('; ') if ': ' in item)

def process_pathway_data(path: str, symbol_to_entrez: Dict[str, int]) -> pd.DataFrame:
    data = []

    for url, description, genes in read_gmt(path):
        filtered_genes = {symbol_to_entrez.get(gene) for gene in genes} - {None}

        if any(gene in symbol_to_entrez for gene in genes):
            data.append({
                'identifier': Path(url).stem,
                'name': parse_description(description).get('name'),
                'url': url,
                'genes': filtered_genes
            })

    return pd.DataFrame(data)

def process_wikipathways_data(gmt_data, human_genes: set):
    wikipath_df = pd.DataFrame(gmt_data, columns=['name', 'description', 'genes'])
    wikipath_df['name'] = wikipath_df['name'].str.split('%').str[0]
    wikipath_df['genes'] = wikipath_df['genes'].apply(lambda genes: {int(gene) for gene in genes} & human_genes)
    
    return wikipath_df[wikipath_df['genes'].astype(bool)]

def create_combined_df(pc_df, wikipath_df):
    wikipath_df = wikipath_df.assign(
        identifier=wikipath_df['description'].str.rsplit('/', n=1).str[-1],
        url=wikipath_df['description']
    )

    combined_df = pd.concat([wikipath_df, pc_df], ignore_index=True).drop_duplicates('genes')
    combined_df['n_genes'] = combined_df['genes'].apply(len)
    
    return combined_df[['identifier', 'name', 'url', 'n_genes', 'genes']].sort_values('identifier')

def write_to_node_csv(dataframe: pd.DataFrame, filepath: str) -> None:
    dataframe[['identifier', 'name', 'url']].to_csv(filepath, index=False)

def write_to_edge_csv(dataframe: pd.DataFrame, filepath: str) -> None:
    edge_data = []

    for _, row in dataframe.iterrows():
        pathway_id = row['identifier']
        genes = row['genes'].split('|')

        for gene_id in genes:
            edge_data.append({'source_id': gene_id, 'target_id': pathway_id})

    edges_df = pd.DataFrame(edge_data)
    edges_df.to_csv(filepath, index=False)

def main():
    gene_info_path = 'data/input/Homo_sapiens.gene_info.gz'
    
    with gzip.open(gene_info_path, 'rt') as file:
        gene_df = pd.read_csv(file, delimiter='\t')

    human_genes = set(gene_df['GeneID'])
    symbol_to_entrez = dict(zip(gene_df['Symbol'], gene_df['GeneID']))

    pc_df = process_pathway_data("data/input/PathwayCommons12.All.hgnc.gmt.gz", symbol_to_entrez)
    gmt_data = read_gmt('data/input/wikipathways-Homo_sapiens.gmt')
    wikipath_df = process_wikipathways_data(gmt_data, human_genes)

    pathway_df = create_combined_df(pc_df, wikipath_df)
    
    write_to_node_csv(pathway_df, 'data/output/node_Pathway.csv.gz')
    write_to_edge_csv(pathway_df, 'data/output/edge_Gene_participatesIn_Pathway.csv.gz')

if __name__ == "__main__":
    main()