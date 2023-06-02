# Scripts Components

Scripts that are executed to take data and create an experimental ML model. This document details each component stage of the whole pipeline at a high-level.

For more detailed information, click [here](../docs/component_docs.md)


## ML training

0. `download_pfam.py`

*Download the pfam HMMs via FTP*
    
    Params: N/A

    Inputs: N/A

    Outputs: Pfam-A.hmm in `data/pfam`

    Metrics: pfam_pulled_timestamp

1. `train_model.py`

*Trains the 'first-pass' ML model*

    Params: N/A

    Inputs: N/A

    Outputs: N/A

    Metrics: Not implemented yet



## Utilizing the model

1. `user_input.py`
*Lorem ipsum dolor sit amet, consectetur adipiscing elit*

    Params: N/A

    Inputs: N/A

    Outputs: N/A

    Metrics: Not implemented yet