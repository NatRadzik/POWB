Final project done for the course "Computetional pipelines in Bioinfromatics".

Project is using data from this publication: https://pmc.ncbi.nlm.nih.gov/articles/PMC5032908/
These are genes from the human X chromosome from 2 populations - from the UK and Nigeria. 

Pipeline outline:
trimming reads (removing adapters and poor-quality sequence ends) -> quality control and generation of a HTML MultiQC report -> indexing -> mapping -> read counting -> visualization using heatmap -> differential gene expression analysis -> pathway enrichment analysis

Project done in Nextflow and R.
