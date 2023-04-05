# Scripts Components

Scripts that are executed to take data and create an experimental ML model. This document details each component stage of the whole pipeline at a high-level.

Currently, this document will be rather sparse until we have a functional pipeline.

For more detailed information, click [here](./component_docs.md)

0-2. `preprocessing.py`

*reformats learn2therm database for use in FAFSA pipeline & statistically samples formatted data for model training*
    
    Params: N/A

    Inputs: Learn2ThermDB

    Outputs: - Augmented database file containing relevant proteins and pairs for FAFSA. Keeps original tables intact.
             nb. It is recommended to have a backup copy of the original database as well as 100 GB of free storage and at least 30 GB of available memory. 
             - snakey plots for analysis

    Metrics: N/A


3.`compute_local_hmmer.py` or `HMMER_API.py`

*runs the HMMER algrothim against the pfam database on all protein pairs specified by the user.*

    Params: sampled sequence from the previous component

    Inputs: protein pair amino acid sequences

    Outputs: HMMER domtblout output files for sequences

    Metrics: N/A

4. `compute_pair_function.py` (**ongoing development**)

*parses and filter the results from pfam*

    Params: component 3 outputs

    Inputs: N/A

    Outputs: functional boolean of pair

    Metrics: N/A


5. `find_structures.py`

*acquires structural information*

    Params: component 3 outputs

    Inputs: N/A

    Outputs: functional boolean of pair

    Metrics: N/A


6. `train_val_script*.py`

*trains a machine learning classifier model to predict functional pairs*

    Params: Pandas dataframe containing sampled data from our protein database from Components 0-2 + metric of interest from Component 3.

    Inputs: - Training: Pandas Dataframe containing known protein pairs with Boolean metric of interest.
            - Testing: Pair of amino acid sequences from Component 3.

    Outputs: - Boolean Prediction (Functional pair: True or Functional pair: False)
             - Confidence score on prediction.

    Metrics: N/A