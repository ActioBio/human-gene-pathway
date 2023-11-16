# Download WikiPathways
wget -O data/input/wikipathways-Homo_sapiens.gmt $(curl -s https://data.wikipathways.org/current/gmt/ | grep -oE 'href="wikipathways-[0-9]+-gmt-Homo_sapiens.gmt"' | awk -F '"' '{print "https://data.wikipathways.org/current/gmt/"$2}' | tail -n 1)

# Download Pathway Commons v12
wget --timestamping --directory-prefix data/input/ https://www.pathwaycommons.org/archives/PC2/v12/PathwayCommons12.All.hgnc.gmt.gz

# Download the Human gene_info.gz file of Entrez Genes
wget --timestamping --directory-prefix data/input/ ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz