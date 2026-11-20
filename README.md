# Proteosim Course Package

This repository simulates a proteomics protocol *in silico*. It focuses on protein digestion and LC-MS analysis using Python functions to mimic an experimental workflow. The code is under active development.

## Modules and workflow overview

- input: FASTA files
- output: m/z tables, fragment spectra


- `file_handling.py`: load FASTA file and map protein IDs to sequences. 
- `protein_digestion.py`: simulate the digestion of multiple proteins with a protease of choice, compute sequence coverage
- `liquid_chromatography.py`: predict LC retention times using the `pyteomics` package and select a retention window
- `mass_spectra_simulation.py`: simulate and plot MS1 and MS2 spectra

Please refer to `ms_experiment_final.ipynb` for a detailed tutorial.

## Installation

To install the package, clone the repository and run: 

- `pip install -r requirements.txt`, followed by
- `pip install .`

