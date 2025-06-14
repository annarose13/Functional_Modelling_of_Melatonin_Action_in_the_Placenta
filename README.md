README
================
Anna Davies
2025-06-14

# Functional Modelling of Melatonin Action in the Placenta

This repository contains the code and essential data used in the final
Master’s thesis titled *“Functional Modelling of Melatonin Action in the
Placenta.”* This project investigates melatonin’s role in placental
function using a series of computational modelling techniques.

Please note that only the smaller data files are included here; the full
datasets have been excluded. All genetic data used in this project has
been sourced from publicly available secondary data sources, which are
referenced in the accompanying dissertation.

To ensure proper understanding and execution of the code, it is
recommended to first read the Material and Methods section (Section 2)
of the accompanying dissertation.

<br>

## Repository Structure

- `porcine_data/` - This directory contains the script used to process
  and analyse the porcine bulk RNA sequencing data, as well as
  information on the samples.

- `first_trimester/` - This directory contains the script used to
  process and analyse the first trimester single-cell RNA sequencing
  data.

- `third_trimester/` - This directory contains the script used to
  process and analyse the third trimester single-cell RNA sequencing
  data.

- `hypergraphs/` - This directory contains an example script of how each
  hypergraph was created.

- `README.md` - This file.

- `.gitignore` - git exclusions

<br>

## Getting Started

1.  To use this repository, clone your local machine:

    `git clone https://github.com/annarose13/Functional_Modelling_of_Melatonin_Action_in_the_Placenta`

2.  To run the R scripts (R), open an RStudio project within your
    locally saved repository.

3.  To run the Jupyter Notebooks (Python), ensure that JupyterLab is
    installed (instructions available at:
    <https://github.com/jupyterlab/jupyterlab/blob/main/README.me>)

<br>

### It is recommended to run the scripts in the following order:

1.  `porcine_data/Porcine_Bulkseq_Preprocessing_and_Hypergraph_Analysis.R`

2.  `first_trimester/First_Trimester_scRNAseq_Hypergraph_Analysis.ipynb`

3.  `third_trimester/Third_Trimester_scRNAseq_Analysis.ipynb`

<br>
