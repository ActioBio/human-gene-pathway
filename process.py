import csv
import collections
import json
import requests
import pandas as pd

def fetch_data(url, dtype=None):
    return pd.read_table(url, dtype=dtype)

def fetch_json(url):
    return json.loads(requests.get(url).text)

def read_gmt(path):
    with open(path) as read_file:
        reader = csv.reader(read_file, delimiter='\t')
        for row in reader:
            if row:
                yield row[0], row[1], set(row[2:])

def parse_description(description):
    return dict(item.partition(': ')[::2] for item in description.split('; '))

def process_pathways(path, symbol_to_entrez):
    rows = []
    PC_Row = collections.namedtuple('PC_Row', ['identifier', 'name', 'source', 'genes'])
    for identifier, description, genes in read_gmt(path):
        desc_dict = parse_description(description)
        genes = {symbol_to_entrez.get(gene) for gene in genes}
        genes.discard(None)
        if genes:
            rows.append(PC_Row(f'PC12_{len(rows) + 1}', desc_dict.get('name'), desc_dict.get('datasource'), genes))
    return pd.DataFrame(rows)

def process_wikipathways_data(gmt_generator, human_genes):
    wikipath_df = pd.DataFrame(gmt_generator, columns=['name', 'description', 'genes'])
    wikipath_df.name = wikipath_df.name.map(lambda x: x.split('%')[0])
    for genes in wikipath_df.genes:
        genes &= human_genes
    return wikipath_df[wikipath_df.genes.map(bool)]

def create_combined_df(pc_df, wikipath_df, human_coding_genes):
    wikipath_df = pd.DataFrame({
        'identifier': wikipath_df['description'].map(lambda x: x.rsplit('/', 1)[1]),
        'name': wikipath_df['name'],
        'url': wikipath_df['description'],
        'source': 'wikipathways',
        'license': 'CC BY 3.0',
        'genes': wikipath_df.genes
    })
    pathway_df = pd.concat([wikipath_df, pc_df])
    pathway_df = pathway_df[['identifier', 'name', 'url', 'source', 'license', 'genes']]
    pathway_df.genes = pathway_df.genes.map(frozenset)
    pathway_df = pathway_df.drop_duplicates(['genes'])
    pathway_df['coding_genes'] = pathway_df.genes.map(lambda x: x & human_coding_genes)
    pathway_df.insert(3, 'n_genes', pathway_df.genes.map(len))
    pathway_df.insert(4, 'n_coding_genes', pathway_df.coding_genes.map(len))
    return pathway_df.sort_values('identifier')

def write_to_tsv(dataframe, filepath):
    write_df = dataframe.copy()
    join = lambda x: '|'.join(map(str, sorted(x)))
    for column in ['genes', 'coding_genes']:
        write_df[column] = write_df[column].map(join)
    write_df.to_csv(filepath, index=False, sep='\t')

def main():
    # Fetching Entrez Gene data
    entrez_df = fetch_data('https://raw.githubusercontent.com/dhimmel/entrez-gene/a7362748a34211e5df6f2d185bb3246279760546/data/genes-human.tsv', {'GeneID': str})
    human_genes = set(entrez_df.GeneID)
    human_coding_genes = set(entrez_df[entrez_df.type_of_gene == 'protein-coding'].GeneID)

    # Fetching symbol to Entrez Gene mapping
    symbol_to_entrez = fetch_json('https://raw.githubusercontent.com/dhimmel/entrez-gene/a7362748a34211e5df6f2d185bb3246279760546/data/symbol-map.json')
    symbol_to_entrez = {k: str(v) for k, v in symbol_to_entrez.items()}

    # Processing Pathway Commons data
    pc_df = process_pathways("data/input/PathwayCommons12.All.hgnc.gmt", symbol_to_entrez)

    # Processing WikiPathways data
    gmt_generator = read_gmt('data/input/wikipathways-Homo_sapiens.gmt')
    wikipath_df = process_wikipathways_data(gmt_generator, human_genes)

    pathway_df = create_combined_df(pc_df, wikipath_df, human_coding_genes)
    write_to_tsv(pathway_df, 'data/output/pathways.tsv')

if __name__ == "__main__":
    main()