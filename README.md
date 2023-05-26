# FAFSA (Function annotated from sequence alignment)
Check out our project's documention via [Read The Docs](https://validprot.readthedocs.io/)!

# Background

The following repository is for the project code associated with the two courses: Data Science Methods for Clean Energy Research (ChemE 545) and Software Engineering for Molecular Data Scientists (ChemE 546) at UW.

FAFSA is a project developed by Humood Alanzi, Ryan Francis, Amin Mosallenejad, Logan Roberts, and Chau Vuong.

Purpose: This package is developed to validate functionality between a pair of protein sequences.

# Table of contents

- [Background](#background)
- [Overview](#overview)
- [Requirments](#requirments--installation)
- [Workflow](#workflow)
- [Outputs](#outputs)

# Overview

Protein pair validation is time consuming and resource intensive, given that proteins can be related through many unique functions, both direct and inferred. Many unique softwares specialize in characterizing protein based on a few of these functions. Our pipeline aims to combine different softwares, spanning sequence alignment, structure and folding prediction, and residue conservation into a single pipeline to improve prediction quality and streamline the characterization process.

# Requirments & Installation
To create and activate the environment specified in `environment.yml`
and install the FAFSA package, do the following commands:
```
conda env create --file environment.yml
conda activate validprot
pip install .
```

FAFSA is dependent on Python 3.11.
This package requires a conda environment with external dependencies of biopython and hmmer.
For a more detailed exposition on the modular/importable code, please see [packages](./docs/package_components.md).


# Workflow

1) Retrieve data from Learn2Therm DB
2) Sample from large DB, select features, format data.

    a) sampling notebook found in [notebooks](./notebooks/)
3) Family identification with Pfam

    a) examples of generating outputs from Pfam in [examples](./examples/local_hmmer_example.ipynb)
4) Structure-based family identification with Fatcat + other softwares
5) Develop model to predict protein pair functionality. Include engineered features in addition to features selected in 1).

    a) examples of model development in [notebooks](./notebooks/)
6) Report whether each protein pair is functional along with confidence metrics.

# Outputs

Boolean prediction of whether protein pair is functional.
Text file with confidence statistics.

# Community Guidelines

Our software is open-source. We recommend submission of feature requests and report bugs.