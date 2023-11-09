# MARS Streamlit Application

## Overview
This Streamlit application integrates the functionalities of MARS (Microbiome Abundance and Relative Species) and ANT (Automated NCBI Taxonomy) to provide a comprehensive tool for microbiome researchers. It automates the processing of metagenomic datasets, taxonomic name resolution, and mapping to resources like AGORA2, thereby facilitating a streamlined workflow for microbiome data analysis.

## Features
- **MARS**: Automated extraction and normalization of microbial abundances from 16S rRNA and shotgun metagenomic datasets.
- **ANT**: Automated identification and mapping of taxonomic names, including homosynonyms, to databases such as AGORA2.

## Installation
To run the application on your local machine, you will need to have Python installed, along with the necessary libraries. Follow these steps to get started:

1. Clone the repository: `git clone `

2. Navigate to the cloned repository's directory: `cd MARS`

3. Install the required Python packages: `pip install -r requirements.txt`

## Usage
After installation, you can run the application using Streamlit:

The application will open in your default web browser, or you can access it at `http://localhost:8501`.

## Input Data
The application accepts input data in specific formats. For MARS, a feature table and taxonomy table, or an already combined table are required. Please refer to the sample input files provided:
- MARS Input: `SI_table_S1.csv`
- ANT Input: `SI_table_S2.csv`

## Contributing
Contributions to improve the application are welcome. 

## Support
For support, please open an issue in the GitHub repository or contact ines.thiele@universityofgalway.ie.

## Authors
- Tim Hulshof - Initial work on MARS and ANT
- Bram Nap - Contributions to [specific features]
- Filippo Martinelli - Contributions

## How to Cite
If you use the MARS Streamlit application in your research, please cite it as follows: