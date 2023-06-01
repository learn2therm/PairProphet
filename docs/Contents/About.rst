Overview
=============

Protein pair validation is time consuming and resource intensive, given that proteins can be related through many unique functions, 
both direct and inferred. Many unique softwares specialize in characterizing protein based on a few of these functions. Our 
pipeline aims to combine different softwares, spanning sequence alignment, structure and folding prediction, and residue conservation 
into a single pipeline to improve prediction quality and streamline the characterization process.

Requirements and Installation
=============

To create and activate the environment specified in `environment.yml` in our 'Github <https://github.com/learn2therm/ValidProt>_:

To install the PairProphet package, use the following commands:

::

    conda env create --file environment.yml

::  
    
    conda activate validprot

::

    pip install .


PairProphet is dependent on Python 3.11.
