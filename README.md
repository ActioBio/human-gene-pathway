# Human Gene Pathway

This repository processes and combines human gene-pathway data from PathwayCommons and WikiPathways, creating a unified dataset for genomic analysis.

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
- PathwayCommons12.All.hgnc.gmt
  - This file from Pathway Commons provides comprehensive pathway data, including gene interactions and pathway information. 
- Homo_sapiens.gene_info.gz
  - The file from NCBI is a compressed archive containing detailed information on genes.

### Output
- human-gene-pathway.tsv
  - The file contains combined human gene-pathway data from PathwayCommons and WikiPathways.

## License

All original content in this repository is released as [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) (public domain). WikiPathways data is [licensed](http://www.wikipathways.org/index.php/WikiPathways:License_Terms) as [CC BY 3.0](http://creativecommons.org/licenses/by/3.0/). Reactome data is licensed as [CC BY 4.0](http://creativecommons.org/licenses/by/4.0/). PID data is in the public domain.