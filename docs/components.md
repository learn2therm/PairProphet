# Scripts Components

Scripts that are executed to take data and create an experimental ML model. This document details each component stage of the whole pipeline

1. `s1.0_db_formatter.py`

*format dataframe per user's needs*
    
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

2. `s2.0_Sample_selection.py`

*strategically sample to best train our model and quere the avialble software tools at adequate speeds*
        - Compute PDF's of population
        - Determine new sampling ratios
        - Figure out how to cover the entire space of feature (make good training material)
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

3.`s3.0_pfam_input.py`

*give AA sequence to pfam, get basic results*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

4. `s3.1_pfam_parse.py`

*parse and filter the results from pfam*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

5.`s3.2_pfam_apply.py`

*apply results to pair data*
    
    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

6.`s3.3_pfam_compute.py`

*compute metric of interest*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

<nb. we can have even more steps here re: 3D structure prediction>

7. `s4.0_data_collection.py`

*collect data from various softwares*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

8. `s5.0_relation.py`

*Create experimental dataset from outputs of collect (pfam, 3D, etc.) that can be used to formulate ML model*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something

9. `s6.0_filter.py`

*filter the learn2therm dataset based on this analytical tool*

    Params: something

    Inputs: something

    Outputs: something

    Metrics: something