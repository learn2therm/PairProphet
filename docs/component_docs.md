# Pipeline documentation
This document offers a comprehensive exposition of all the components as well as the subcomponents. For a high-level overview, please refer to [components](./components.md). The purpose, inputs, steps, outputs, metrics, are all described for each pipeline script in this document.


# Table of contents

- [Component 0](#component-0)
- [Component 1](#component-1)
- [Component 2](#component-2)
- [Component 3](#component-3)
- [Component 4](#component-4)
- [Component 5](#component-5)
- [Component 6](#component-6)

# Component 0

    **Params:** 

        **Inputs:** learn2therm protein pair database file

        **Outputs:** Augmented database file containing relevant proteins and pairs for FAFSA. Keeps original tables intact.
                     It is recommended to have a backup copy of the original database as well as 100 GB of free storage and at 
                     least 30 GB of available memory. 

    Component 0 is specifically concerned with the learn2therm database and is not intended for use by users with other data files.
    Additionally, this code need only be executed once to generate all files necessary for downstream processing. Users can input
    their own similarly formatted data for classification starting in component 1. Component 0 can feed both component 1 and 
    component 2. Additionally, users have the option to export intermediate metadata on the flow of information between learn2therm
    and FAFSA as Sankey plots.

    **Packages:** os, time, pandas, numpy, duckdb, plotly, kaleido

### **Subcomponent 1**: 

**Use case**: 

        User has copy of learn2therm database file and generates data for FAFSA. 

**Test**: 

        Test that connection and queries to new FAFSA database are correct. Test connection to old tables as well as new
        within connection object.

### **Subcomponent 2**: 

**Use case**: 

        Produces basic analysis of database transformation via Sankey plots.

**Test**: 

        Test that plot files are created and not empty.

=======
    **Packages:** os, sys, time, pandas, numpy, duckdb, plotly, kaleidoscope

# Component 1

    **Params:** 

        **Inputs:** Protein pair database file or pandas dataframe formatted according to instructions in README.
    
                    Input file should contain at minimum:

                        1. AA sequence for both members of each pair.
                        2. Clear indexing/column names to identify which belongs in meso-thermo groups.
                        3. Most recent metrics for model training (e.g. alignment/coverage scores). Included in docs
                           and c1 docstrings.

        **Outputs:** Pandas dataframe(s) containing data sample of specified size. We hope to include analytics as a
        separate output in the future to inform user decisions about sampling.

    Component 1 formats data for the FAFSA model. The resulting dataframe is pruned to contain only features and 
    target data.

    **Packages:** os, sys, pandas
    
### **Subcomponent 1**: 

**Use case**: 

        User supplies path to database or .csv file and desired sample size. Pulls data into DataFrame using one of several sampling
        methods. 

**Test**: 

        Test assert that DataFrame was generated, compare shape with desired result. Check for obvious errors such as inconsistent data type within columns. Make sure all sampling configurations function as expected.
          

# Component 2

    **Params:** 

        **Inputs:** Protein pair database file or pandas dataframe formatted for FAFSA. Sampling 
                    parameters for training set.

        **Outputs:** DataFrame containing training data. Exported analytics and plots for user reference (optional, not implemented 
                     Winter 2023).

    Component 2 is a sampler for training the learn2therm database. Users can select the desired sampling method, and the function
    can seek data subsets satisfying those parameters. For example, data can be represented with higher information density 
    (fringe data is favored) using the Frank-Wolfe algorithm.

    **Packages:** pandas, numpy
    
### **Subcomponent 1**: 

**Use case**: 

        User supplies large dataset. Returns sampled dataset according to input parameters. Currently supports random sampling
        or Frank-Wolfe D-optimal sampling.

**Test**: 

        Test that sample distributions are significantly different than random. Test that convergence conditions are met and that output dataframe is the correct size.


# Component 3
## Software Component Three: compute_local_hmmer.py or HMMER_API.py
    
    **Params:** sampled sequence data from component 2

    **Inputs:** protein pair amino acid sequences

    **Outputs:** for the two proteins in a pair, you've a target family based on hmmer metrics and a boolean to check if a pair are functional or not according to pfam parsing

    **Metrics:**

    **Packages:** pandas, biopython, HMMER, pfam database (database, optional), unittest

Component 3 aims to use the HMMER algrothim running against the pfam database on all protein pairs specified by the user. In this case, we are using the Learn2thermDB protein pairs.
This component has two options to be ran locally or using the online HMMER server API. The two options have different subcomponents, use-cases and tests. However, for the purposes of this documentation, we will assume that the two options are the same, which is a reasonable assumption to make in this case. This component takes in the sampling data from the upstream component 2. From there, the component starts ensuring that the amino acid sequences supplied by the user are read and written in the appropriate .fasta file type for HMMER to take in. Then, it runs HMMER search against the Pfam database using the hmmscan command. Here, we aim to embarrassingly parallelize to beat the I/O disk-reading of the HMMER algrothim, and make sure that it can take as much sequences as possible. Unfortuntely, this is subject of ongoing work and has not been implemented properly yet. This component will undergo significant further development during the spring. We plan to run HMMER on the newest data sample and we will work on a function that parses the HMMER output, filters, and determines if a pair are functional or not, which will be the input for component 5.

### **Subcomponent 1**: Read sequence

**Use case**: 
        
        User sends their AA sequences of interest. Transform inputs into appropriate fasta seqeunces, and write them as an input file.

**Code**:

        Reads the pair AA sequence returns a list of biopython SeqRecord objects and creates a corresponding input Fasta of them

**Test**: 

        1) Tests if the provided AA sequence is valid (canonical AAs only)
        2) Tests if the input file is created
        3) Tests if the given input from file to original sequence dataframe are empty, invalid, improper type, etc.

### **Subcomponent 2**: Run HMMER against pfam database

**Use case**: 

        Users' input files are provided to HMMER using the hmmscan command to run against the pfam database. Depending on the numbers of the pairs, parallelization will be included here. 

**Code**:

        1) Creates a corresponding HMMER output file from the hmmscan subprocess for the inputs of the user
        --optional (in-development)--
        2) Chunks a large list of sequences based on the class of two pairs
        3) run hmmscan on each chunk using embarrassingly parallel python techinque  

**Test**: 

        Given that we are not using our algrothim here, we are mainly testing if an output file is created


### **Subcomponent 3**: parse and compute functional pair from HMMER outputs (in-development)

**Use case**: 

        parses and filters the various hmmer outputs for the pairs and decides using various methods if the two pairs are functional or not 

**Code**:

        This is still in active-development 

**Test**: 

        Work-in-progress


#### Plan Outline
1. for the first option, explicitly integrate with the upstream and downstream components especially when we have a small dataset
2. for the second approach, parse, filter, and compute metrics of interest for a large dataset

# Component 4 - Acquire structural information
## Software Component Four: find_structures.py

    **Params:**

    **Inputs:** DNA/RNA/protein sequence (FASTA) or Pfam IDs (domain information)

    **Outputs:** DataFrame containing proteins that share sequence or structure similarity to query sequence. Desired parameters include e-value or bitscore and percentage similarity for example.

    **Metrics:**

    **Packages:** pandas, json, pypdb

Component 4 is responsible for searching the RCSB Protein Data Bank (PDB) to identify proteins that share structural or sequence similarity with the query protein. This component first searches for experimentally solved structures, and if not available, then searches for computationally generated structures. The output of this component also includes proteins with sequence similarity to the query protein, allowing for comparison with the results from Pfam (component 3).

This component is in developing phase this quarter as we will only ensure that the pfam functionality works first, so there will be no collecting. Once 3D structure prediction comes in, we formally work on this.

### **Subcomponent 1**: Send request and get results

**Use case**: 

        Searches using single DNA/RNA/protein sequence (FASTA) or Pfam IDs (can be one or multiple Pfam ID(s)) and gets a list of PDB IDs or a more comprehensive lists of result (includes PDB IDs, e-values, bitscores, similarity, etc.) Currently using PyPDB for this.

**Test:** 

        Proper import of PyPDB and expected query formats.

### **Subcomponent 2**: Extract desired information

**Use case**: 

        Flatten the nested dictionary and then filter it. This involves recursively iterating through the nested dictionary to extract all key-value pairs and converting them into a single flat dictionary. Once the nested dictionary is flattened, we can then filter the resulting dictionary to keep only the desired information

**Test:** 

        Absences of dropped parameters.

### **Subcomponent 3**: Format output into DataFrame 

**Use case:** 

        Convert Python dictionary to DataFrame. Manipulate the DataFrame using pandas functions to sort, filter, or group the data.

**Test:** 

        Expected shape and size with correct columns DataFrame.

# Component 5 - Training and Validation
## Software Component Five: train_val_script.py

    **Params:** Pandas dataframe containing sampled data from our protein database from Components 0-2 + metric of interest from Component 3.

    **Inputs:** Training: Pandas Dataframe containing known protein pairs with Boolean metric of interest.
                Testing: Pair of amino acid sequences from Component 3.

    **Outputs:**  Boolean Prediction (Functional pair: True or Functional pair: False)
                  Confidence score on prediction.

    **Metrics:**

    **Packages:** pandas, matplotlib, numpy, scikit-learn (.ensemble, .preprocessing, .model_selection, neighbors, .feature_selection), unittest
    
The purpose of this component is to train a machine learning classifier model on our large dataset of protein pairs. 

For training, the model will take in a pandas dataframe containing protein pairs as examples. The database will be cleaned until 10 quantitative features remain, including alignment metrics, optimal growth temperature, bit scores, sequence length, and alignment coverage. The cleaned dataframe will also include a Boolean target, classifying the protein as either a functional pair or not a functional pair. The model will be trained and tested based on this data.

The upstream Pfam parsing component requires more development before data scraped from Pfam can be used to train and validate the model. Right now, it is trained and validated from a n=50,000 sample from Learn2Therm database. Ongoing work is being done regarding integrating this component with the upstream.

This component will undergo significant further development during the spring. We plan to add functionality such that it can use data scraped from other protein databases than Pfam as features, including tertiary and quarternary structural data. Furthermore, we plan to improve our confidence in predicition and eventually develop a quantitative functionality prediction.


### **Subcomponent 1**: Test for pandas dataframe input

**Use case**: 

        User takes data from component 3 (where data is processed into pandas dataframe) and wants to pass it into relationship component.


**Test**: 
        
        Test asserts that data structure received is a pandas dataframe.


### **Subcomponent 2**: Checks that input data is cleaned properly.

**Use cases**: 

        1) Input data does not include bit score, which we need as an input to our model.

        2) Input data has NaN values.

        3) Input data is missing either a mesophillic or thermophillic protein sequence.
               
**Code**: 

        1) Function that will remove unwanted features from dataframe (protein ID, query ID, 16_s coverage). 

        2) Function that will drop all rows with NaN values. In the future, we may add a functionality that allows for examples with NaN values in certain columns to remain in the dataframe.

        3) Function that verifies two amino acid sequences exist for each example in the dataframe.

**Test**: 

        1) After unwanted columns are removed, assert statement raises error if dataframe has any unwanted features. Separate assert same ensures that all of the features necessary for training are present in the dataframe.

        2) Function that asserts no NaN values exist in the dataframe.

        3) Test than takes an example that is missing one of the sequences and raises an error.


### **Subcomponent 3**: Train the model with sample data. 

**Use case**: 

        Model to predict eventual user input needs to be trained on the protein pair database that we currently have.

**Code**: 

        Function takes in clean pandas dataframe from subcomponent 2. Splits data into dev and test set between the features and target. Trains the data with a RandomForestClassifier model from sklearn. Returns the trained model along with development and test dataframes.
          
          Input: Pandas dataframe
          Output: Model, Pandas dataframe

**Test**: 

        Asserts that input data type is a pandas dataframe and that output has 5 objects (model, dev_X, dev_y, test_X, test_y).


### **Subcomponent 4**: Test the model with sample data. 

**Use case**: 

        Uses test dataframe to evaluate the performance of the model and return a vector of Boolean predictions.

**Code**: 

        Function takes testing dataframe along with the model trained in the previous component. Runs testing dataframe through a RandomForestClassifier and returns a vector of predictions.

              Input: Model, Pandas dataframe
              Output: Numpy array

**Test**: 

        Asserts that input is a pandas dataframe, output is a numpy array, and the length of the output vector is equal to the number of examples.


### **Subcomponent 5**: Plot a confusion matrix of the test results.

**Use case**: 

        User wants to visualize the results of the model testing and receive a confidence score on the prediction.

**Code**: 

        Function takes testing dataframe along with the model trained in the previous component. Plots confusion matrix using sklearn.metrics library. Returns model score.

              Input: Model, Pandas dataframe
              Output: Numpy array, Confusion Matrix
              
**Test**:

        N/A. Looking to eventually develop a test using python image prediction module.

### **Subcomponent 6**: Calculate a 'functionality' metric that is the ultimate output of component five.

Our goal is to develop some functionality score between a pair of proteins that is not just based on our singular classifier model. This will factor in information from multiple softwares as opposed to just Pfam, which is the focus during the first part of the project. This will be built during spring quarter.

**Use case:**

        User wants a reliable and precise metric to assess how closely related their input proteins are.

**Test**: 

        work-in-progress

#### Plan Outline
1. Get data from component 3. This should already be in a pandas dataframe (data prep is included in C3).
2. Clean the data to prepare it for model training and testing.
3. Train and test the model, return scores, plot, and any other necessary indicator of model performance.
4. Input new user data and return a functionality score for the input protein pair.


# Component 6
## Work-in-progress
We do not anticpate to reach this component this Winter quarter, but the Spring quarter! 

