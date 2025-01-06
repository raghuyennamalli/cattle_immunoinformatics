# cattle_immunoinformatics
# Variant Data for "Genetic basis of immunity in Indian cattle as revealed by comparative analysis of Bos genome"

## Description
This repository contains the common variant data identified by both GSAlign and SyRI tool associated with the manuscript titled "Genetic basis of immunity in Indian cattle as revealed by comparative analysis of Bos genome", submitted to Scientific Reports. The data represents genomic variations identified between _Bos indicus_ (Nelore breed) and _Bos taurus_ (Hereford breed) with focusing on identifying significant immune genes marker in indicine (Indian) cattle, as indicine cattle have inherent disease resistance capability in compared to taurine cattle (Exotic) .

## Contents
- `GSAlign SyRI common variants`: The folder contain the variation (such as Insertion, Deletion, and Substitution) data (as CSV file) on chromosome wise, identified between Herefored (Reference) and Nelore breed (Query).
- `Common Deletion variant in immune-related genes, Common Insertion variant in immune-related genes, and Common Substitution variant in immune-related genes`: This CSV file containing the common variation present only in the immune-related genes. (The immune genes identified using InnateDB database and keyword-search)
- `Gene Extraction and Gene Extraction QTL mapping.py`: This Python script is used to extract the genes have common chromosomsal variation and its associated QTLs in Hereford breed
- `SyMAP HitID60`: This CSV file containing the genes have low sequence similarity (HIT:less than 60%) while comparing the whole genome sequence of the Hereford and Nelore breed using SyMAP tool
- `SyMAP Chromosomal1_X variant`:This is the output CSV file obtained from SyMAP tool while comparing the whole genome sequence of the Hereford and Nelore breed    
- `README.md`: Information about the variant data obtained from the pipeline we developed.

## File Information
- **Format:** Variant data in CSV format
- **Data Source:**The whole-genome assembly data were downloaded from the NCBI in both GenBank and FASTA formats for the three breeds, namely Nelore (Accession ID: GCA_000247795.2), Gir (Accession ID: GCA_002933975.1), and Hereford (Accession ID: GCA_002263795.3)
- **Variant identification:** Variant identified using GSAlign, SyRI, and SyMAP
- **Annotation:** Variants were annotated using in-house Python script

- ## Author

**Ragothaman M Yennamalli**  
Senior Assitant Professor, 
Department of Bioinformatics,
School of Chemical and Biotechnolgy,
SASTRA Deemed University 
ragothaman@scbt.sastra.edu
Your ORCID ID:  
GitHub Profile: https://github.com/raghuyennamalli

For any questions or further information, feel free to contact me at ragothaman@scbt.sastra.edu.
