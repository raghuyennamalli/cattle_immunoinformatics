# Data for "Genetic basis of immunity in Indian cattle as revealed by comparative analysis of Bos genome"

Menaka Thambiraja1, Shukkruthi K. Iyengar1$, Brintha Satishkumar1$, Sai Rohith Kavuru1$, Aakanksha Katari1$, Dheer Singh2, Ragothaman M. Yennamalli1*, & Suneel K. Onteru2*
1Department of Bioinformatics, School of Chemical and Biotechnology, SASTRA Deemed to be University, Thanjavur, India. 2Molecular Endocrinology, Functional Genomics and Systems Biology Laboratory, Animal Biochemistry Division, National Dairy Research Institute, Karnal, India. *E-mail: ragothaman@scbt.sastra.edu
$ Equal contribution

## Description
This repository contains the genomic variants, such as single-nucleotide variations (SNV), insertions, and deletions identified by the GATK tool associated with the manuscript titled "Genetic basis of immunity in Indian cattle as revealed by comparative analysis of the Bos genome", submitted to Scientific Reports. The data represent genomic variations identified in _Bos indicus_ (Nelore breed, n=14) and _Bos indicus_ (Gir breed, n=20) with a focus on identifying genomic variation located on the immune-related genes, on the basis that indicine cattle have inherent disease resistance capability compared to taurine cattle (Exotic breed).

## Contents
- `GSAlign SyRI common variants`: The folder contain the variation (such as Insertion, Deletion, and Substitution) data (as CSV file) on chromosome wise, identified between Herefored (Reference) and Nelore breed (Query).
- `Common Deletion variant in immune-related genes, Common Insertion variant in immune-related genes, and Common Substitution variant in immune-related genes`: This CSV file containing the common variation present only in the immune-related genes. (The immune genes identified using InnateDB database and keyword-search)
- `Gene Extraction and Gene Extraction QTL mapping.py`: This Python script is used to extract the genes have common chromosomsal variation in Nelore breed and its associated QTLs in Hereford breed (as QTLs for Nelore breed is not available)
- `SyMAP HitID60`: This CSV file containing the genes have low sequence similarity (HIT:less than 60%) while comparing the whole genome sequence of the Hereford and Nelore breed using SyMAP tool
- `SyMAP Chromosomal1_X variant`: This is the output CSV file obtained from SyMAP tool while comparing the whole genome sequence of the Hereford and Nelore breed    
- `README.md`: Information about the variant data obtained from the pipeline we developed.

## File Information
- **Format:** Variant data in CSV format
- **Data Source:** The whole-genome assembly data were downloaded from the NCBI in both GenBank and FASTA formats for the three breeds, namely Nelore (Accession ID: GCA_000247795.2), Gir (Accession ID: GCA_002933975.1), and Hereford (Accession ID: GCA_002263795.3)
- **Variant identification:** Variant identified using GSAlign, SyRI, and SyMAP
- **Annotation:** Variants were annotated using in-house Python script

**The manuscript has been preprinted:** https://www.biorxiv.org/content/10.1101/2024.12.09.627532v1

- ## Author

**Ragothaman M Yennamalli**  
Senior Assitant Professor, 
Department of Bioinformatics,
School of Chemical and Biotechnolgy,
SASTRA Deemed University 
Email: ragothaman@scbt.sastra.edu
ORCID ID: 0000-0002-3327-1582

For any questions or further information, feel free to contact me at ragothaman@scbt.sastra.edu.
