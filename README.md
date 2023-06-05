# PairProphet
Check out our project's documention via [Read The Docs](https://pairprophet.readthedocs.io/)!

# Background

The following repository is for the project code associated with the three courses: Data Science Methods for Clean Energy Research (ChemE 545) Software Engineering for Molecular Data Scientists (ChemE 546), and Molecular Data Science Capstone (ChemE 547) and at UW.

PairProphet is a project developed by Humood Alanzi, Ryan Francis, Amin Mosallenejad, Logan Roberts, and Chau Vuong.

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
and install the PairProphet package, do the following commands:
```
conda env create --file environment.yml
conda activate pairpro
pip install .
```

PairProphet is dependent on Python 3.11.
This package requires a conda environment with external dependencies of biopython and hmmer.
For a more detailed exposition on the modular/importable code, please see [component docs](./docs/component_docs.md).


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
Check out our [code of conduct](./docs/code_of_conduct.md) for more information.
Then, please see our [contributing guidelines](./docs/contributing.md) for more information.