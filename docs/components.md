# Scripts Components

Scripts that are executed to take data and create an experimental ML model. This document details each component stage of the whole pipeline at a high-level.

Currently, this document will be rather sparse until we have a functional pipeline.

For more detailed information, click [here](./component_docs.md)

0. `c0.0_db_formatter.py`

*reformats learn2therm database for use in ValidProt pipeline*
    
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

1. `c1.0_test_input.py`

*Input pipeline for test datasets*
    
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

2. `c2.0_sampling.py`

*statistically samples formatted data for model training*
        - Compute PDF's of population
        - Determine new sampling ratios
        - Figure out how to cover the entire space of feature (make good training material)
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

3.`c3.0_pfam_input.py`

*give AA sequence to pfam, get basic results*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

4. `c3.1_pfam_parse.py`

*parse and filter the results from pfam*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

5.`c3.2_pfam_apply.py`

*apply results to pair data*
    
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

6.`c3.3_pfam_compute.py`

*compute metric of interest*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

<nb. we can have even more steps here re: 3D structure prediction>

7. `c4.0_data_collection.py`

*collect data from various softwares*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

8. `c5.0_relation.py`

*create experimental dataset from outputs of collect (pfam, 3D, etc.) that can be used to formulate ML model*


    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

9. `c6.0_filter.py`

*filter the learn2therm dataset based on this analytical tool*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something