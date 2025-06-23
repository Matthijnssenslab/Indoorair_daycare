# Indoorair_daycare

Scripts and data to reproduce figures from our indoor air viral metagenomics surveillance study.

## Publication

Preprint available: https://doi.org/10.1101/2024.10.11.24315306

**Citation**: Karatas, M., Geenen, C., Keyaerts, E., Budts, L., Raymenants, J., Eggers, C., Craessaerts, B., André, E., & Matthijnssens, J. (2024). Untargeted viral metagenomics of indoor air as a novel surveillance tool for respiratory, enteric and skin viruses. medRxiv. https://doi.org/10.1101/2024.10.11.24315306

## Repository Structure

```
├── scripts_to_reproduce_figures/    # R scripts for analysis and figure generation
│   ├── Figure 1/                   # Scripts for Figure 1 
│   ├── Figure 2/                   # Scripts for Figure 2 
│   ├── Figure 3/                   # Scripts for Figure 3 
│   ├── Figure 4/                   # Scripts for Figure 4 
│   └── Figure 5/                   # Scripts for Figure 5 
├── input/                          # Input data files
│   ├── csv_files/                  # CSV data files for analysis
│   └── input_phyloseq/             # Phyloseq objects for de-novo assembly analysis
└── output/                         # Generated publication outputs
    ├── PDF_figures/                # Generated publication figures (PDF)
    └── Supplementary_docx/         # Supplementary documents

```

## Abstract

This longitudinal study explores the use of indoor air, in combination with targeted qPCR panels and untargeted viral metagenomics, as a novel virus surveillance tool. Indoor air samples were collected weekly from a daycare center in Leuven, Belgium, over a 12-month period using active indoor air sampling, followed by screening using respiratory and enteric qPCR panels, as well as untargeted viral metagenomics.

Human-associated viruses were detected in 95.2% (40/42) of samples, with MW polyomavirus being the most prevalent at 80.9%. Several other respiratory viruses (e.g., rhinoviruses, RSV-B) and enteric viruses (e.g., rotavirus, astrovirus, adenovirus) were identified, correlating with their known epidemiological circulation patterns.


## Data Availability

Raw data, excluding human reads, is available at the Sequence Read Archive (SRA) with Bioproject number PRJNA1158979. Codes and processed data used to conduct the analyses and produce figures are available in this repository.

## Contact

For questions about the code or data, please open an issue on GitHub or contact: mustafa.karatas@kuleuven.be

## Acknowledgments
Mustafa Karatas is supported by a Research Foundation Flanders (FWO) fundamental research scholarship (number: 11P7I24N). We would like to thank Pieter Wets and Steven Traets for their help in collection of samples. Additionally, we thank the Rotavirus National Reference Center and National Reference Center for Respiratory Pathogens for their support throughout this research.

