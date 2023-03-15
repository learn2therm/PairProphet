# Executables
In this directory, you will find the different sub-modules of FAFSA package:

- Component 0 retrieves data from Learn2Therm database and produces basic data analyses of database transformation via Sankey plots.

- Component 1 formats data for FAFSA pipeline such that it contains only relevant features and target data.

- Component 2 is a sampler for training the Learn2Therm database. Users can seek data subsets that fulfill their parameters and characterization goals.

- Component 3 has two options:

    - `HMMER_API.py` which queries a few amino acid sequeces from a user's database against PFAM and return an interpretable table

    - `compute_local_hmmer.py` upon downloading the pfam-A.hmm db locally, the user can query as much sequences their hardware allows to get HMMER/Pfam functional information

    - `compute_pair_function.py` this script is still in active development, but in the future, it should be able to parse the HMMER outputs, filter them and generate a dataframe if a protein pair are functional or not

- Component 4 is intended to return proteins that have structural similarity to the input by searching the Protein Data Bank (PDB). The current code can search from PDB to return proteins that share sequence homology. This will be leveraged for analyzing structural similarity.

- Component 5 contains a pipeline that cleans the input dataframe, isolates features, and defines target for machine learning. Data is then passed through RandomForestClassifier from scikit-learn.

- Component 6 is a proposed Spring quarter milestone that will include more analysis of functionality between protein pairs.


## Note
Currently these sub-modules/components aren't fully integrated into each other due to temporary technical hurdles. The next version of this software will include integration.

