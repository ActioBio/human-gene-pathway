import csv
import collections
import json

import requests
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# read entrez genes
url = 'https://raw.githubusercontent.com/dhimmel/entrez-gene/a7362748a34211e5df6f2d185bb3246279760546/data/genes-human.tsv'
entrez_df = pd.read_table(url, dtype={'GeneID': str})
human_genes = set(entrez_df.GeneID)
human_coding_genes = set(entrez_df[entrez_df.type_of_gene == 'protein-coding'].GeneID)

url = 'https://raw.githubusercontent.com/dhimmel/entrez-gene/a7362748a34211e5df6f2d185bb3246279760546/data/symbol-map.json'
symbol_to_entrez = json.loads(requests.get(url).text)
symbol_to_entrez = {k: str(v) for k, v in symbol_to_entrez.items()}

def read_gmt(path):
    with open(path) as read_file:
        reader = csv.reader(read_file, delimiter='\t')
        for row in reader:
            if not row:
                continue
            identifier = row[0]
            description = row[1]
            genes = set(row[2:])
            yield identifier, description, genes

def parse_description(description):
    description_dict = {}
    for item in description.split('; '):
        key, _, value = item.partition(': ')
        description_dict[key.strip()] = value.strip()
    return description_dict

i = 0
rows = list()
PC_Row = collections.namedtuple('PC_Row', ['identifier', 'name', 'source', 'genes'])
path = "data/input/PathwayCommons12.All.hgnc.gmt"
for identifier, description, genes in read_gmt(path):
    # Process description
    description_dict = parse_description(description)
    name = description_dict.get('name')
    source = description_dict.get('datasource')
    if description_dict.get('organism') != '9606':
        continue

    # Convert genes to Entrez
    genes = {symbol_to_entrez.get(x) for x in genes}
    genes.discard(None)
    if not genes:
        continue
    
    # Add pathway
    i += 1    
    row = PC_Row(
        identifier = 'PC12_{}'.format(i),
        name = name,
        source = source,
        genes = genes,
    )
    rows.append(row)
    
pc_df = pd.DataFrame(rows)
print(pc_df.head(2))

print(pc_df.source.value_counts())

# Assuming pc_df.genes is a Series of sets, we map each set to its length
gene_lengths = pc_df['genes'].map(len)

# Plotting the KDE plot
plt.figure()  # Create a new figure
sns.kdeplot(data=gene_lengths, fill=True)
plt.xlabel('Gene Set Length')
plt.ylabel('Density')
plt.title('Distribution of Gene Set Lengths')
plt.savefig('data/output/gene_set_length_distribution.png')

# Creating the histogram
plt.figure()  # Create a new figure
plt.hist(gene_lengths, bins=np.arange(101))  # 101 to include the right edge of the last bin
plt.xlim(0, 100)
plt.xlabel('Gene Set Length')
plt.ylabel('Frequency')
plt.title('Histogram of Gene Set Lengths')
plt.savefig('data/output/gene_set_length_histogram.png')

gmt_generator = read_gmt('data/input/wikipathways-Homo_sapiens.gmt')
wikipath_df = pd.DataFrame(gmt_generator, columns=['name', 'description', 'genes'])
wikipath_df.name = wikipath_df.name.map(lambda x: x.split('%')[0])
print(len(wikipath_df))

for genes in wikipath_df.genes:
    genes &= human_genes
wikipath_df = wikipath_df[wikipath_df.genes.map(bool)]
print(len(wikipath_df))

print(wikipath_df.head(2))

# Density plot of genes per pathway
plt.figure()  # Create a new figure for the density plot
sns.kdeplot(wikipath_df.genes.map(len), fill=True)
plt.xlabel('Gene Set Length')
plt.ylabel('Density')
plt.title('Density Plot of Genes per Pathway')
plt.savefig('data/output/density_plot_genes_per_pathway.png')  # Save as PNG

# Histogram of genes per pathway
plt.figure()  # Create a new figure for the histogram
plt.hist(list(wikipath_df.genes.map(len)), np.arange(100))
plt.xlim(0, 100)
plt.xlabel('Gene Set Length')
plt.ylabel('Frequency')
plt.title('Histogram of Genes per Pathway')
plt.savefig('data/output/histogram_genes_per_pathway.png')  # Save as PNG

wikipath_df = pd.DataFrame({
    'identifier': wikipath_df['description'].map(lambda x: x.rsplit('/', 1)[1]),
    'name': wikipath_df['name'],
    'url': wikipath_df['description'],
    'source': 'wikipathways',
    'license': 'CC BY 3.0',
    'genes': wikipath_df.genes
})
print(wikipath_df.head(2))

pathway_df = pd.concat([wikipath_df, pc_df])
pathway_df = pathway_df[['identifier', 'name', 'url', 'source', 'license', 'genes']]
print(len(pathway_df))

# Remove duplicate pathways
pathway_df.genes = pathway_df.genes.map(frozenset)
pathway_df = pathway_df.drop_duplicates(['genes'])
print(len(pathway_df))

pathway_df['coding_genes'] = pathway_df.genes.map(lambda x: x & human_coding_genes)

pathway_df.insert(3, 'n_genes', pathway_df.genes.map(len))
pathway_df.insert(4, 'n_coding_genes', pathway_df.coding_genes.map(len))

pathway_df = pathway_df.sort_values('identifier')
print(pathway_df.head())
print(pathway_df.source.value_counts())

# Create a dataframe for writing as a tsv. Multi-element fields are pipe delimited.
write_df = pathway_df.copy()
join = lambda x: '|'.join(map(str, sorted(x)))
for column in 'genes', 'coding_genes':
    write_df[column] = write_df[column].map(join)

write_df.to_csv('data/output/pathways.tsv', index=False, sep='\t')