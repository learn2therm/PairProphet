# ValidProt

# Background

The following repository is for the project code associated with the two courses: Data Science Methods for Clean Energy Research (ChemE 545) and Software Engineering for Molecular Data Scientists (ChemE 546) at UW.

ValidProt is a project developed by Humood Alanzi, Ryan Francis, Amin Mosallenejad, Logan Roberts, and Chau Vuong.

Purpose: This package is developed to validate functionality between a pair of protein sequences.

# Table of contents

- [Background](#background)
- [Overview](#overview)
- [Requirments](#requirments)
- [Workflow](#workflow)
- [Outputs](#outputs)

# Overview

Protein pair validation is time consuming and resource intensive, given that proteins can be related through many unique functions, both direct and inferred. Many unique softwares specialize in characterizing protein based on a few of these functions. Our pipeline aims to combine different softwares, spanning sequence alignment, structure and folding prediction, and residue conservation into a single pipeline to improve prediction quality and streamline the characterization process.

# Requirments

This package requires a conda environment with external dependencies of biopython and hmmer.
For a more detailed exposition on the modular/importable code, please see [packages](./docs/package_components.md).

# Workflow

1) Retrieve data from Learn2Therm DB
2) Sample from large DB, feature selection, data formatting
3) Family identification with Pfam
4) Family identification with PDB + other softwares
5) Develop model to predict protein pair functionality
6) Optimize scoring methods to improve quality of analysis

# Outputs

Boolean prediction of whether protein pair is functional.
Functionality analysis (ongoing work).
