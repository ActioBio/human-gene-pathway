# Human Gene Pathway

This repository processes and combines human gene-pathway data from PathwayCommons and WikiPathways. It focuses on the associations between human genes and various biological pathways sourced from these platforms, creating a unified dataset for genomic analysis.

### Execution

```
conda env create -f environment.yml 

conda activate human-gene-pathway

bash download.sh

python process.py
```

### Input

- wikipathways-Homo_sapiens.gmt
  - Downloaded from the WikiPathways database, this file contains curated biological pathways for Homo sapiens (humans). Each entry includes pathway information such as gene sets associated with specific biological processes or diseases.
  - WikiPathways is an open-source platform that offers manually curated biological pathways. 
- PathwayCommons12.All.hgnc.gmt
  - This file from Pathway Commons provides comprehensive pathway data, including gene interactions and pathway information.
  - PathwayCommons is a comprehensive collection of pathways from multiple databases. 
- Homo_sapiens.gene_info.gz
  - The file from NCBI is a compressed archive containing detailed information on genes.

### Output
- node_Pathway.csv.gz
  - This file lists all the pathways with their identifiers, names, and URLs.
- edge_Gene_participatesIn_Pathway.csv.gz
  - This file describes the associations between genes and pathways.

Note: These CSV files are formatted for easy import into Neo4j.

## License

All original content in this repository is released as [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) (public domain). WikiPathways data is [licensed](http://www.wikipathways.org/index.php/WikiPathways:License_Terms) as [CC BY 3.0](http://creativecommons.org/licenses/by/3.0/). Reactome data is licensed as [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/). PID data is in the public domain.
