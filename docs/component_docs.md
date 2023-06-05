# Pipeline documentation

This document offers a comprehensive exposition of all the components as well as the subcomponents. The purpose, inputs, steps, outputs, metrics, are all described for all the importable code found in pairpro decribed as components in this document. For a high-level overview of scripts, please refer to [components](../scripts/README.md). The components talked about here constitute the scripts in the scripts folder.

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

        **Inputs:** learn2thermDB protein pair database file

        **Outputs:** Augmented database file containing relevant proteins and pairs for pairpro. Keeps original tables intact.
                     It is recommended to have a backup copy of the original database as well as 100 GB of free storage and at
                     least 30 GB of available memory.

    Component 0 is specifically concerned with the learn2therm database and is not intended for use by users with other data files.
    Additionally, this code need only be executed once to generate all files necessary for downstream processing. Users can input
    their own similarly formatted data for classification starting in component 1. Component 0 can feed both component 1 and
    component 2. Additionally, users have the option to export intermediate metadata on the flow of information between learn2therm
    and pairpro as Sankey plots.

    **Packages:** os, time, pandas, numpy, duckdb, plotly, kaleido

### **Subcomponent 1**:

**Use case**:

        User has copy of learn2therm database file and generates data for pairpro.

**Test**:

        Test that connection and queries to new pairpro database are correct. Test connection to old tables as well as new
        within connection object.

### **Subcomponent 2**:

**Use case**:

        Produces basic analysis of database transformation via Sankey plots.

**Test**:

        Test that plot files are created and not empty.

=======
**Packages:** pandas, numpy, duckdb, plotly, kaleidoscope

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

    Component 1 formats data for the pairpro model. The resulting dataframe is pruned to contain only features and
    target data.

    **Packages:** pandas

### **Subcomponent 1**:

**Use case**:

        User supplies path to database or .csv file and desired sample size. Pulls data into DataFrame using one of several sampling
        methods.

**Test**:

        Test assert that DataFrame was generated, compare shape with desired result. Check for obvious errors such as inconsistent data type within columns. Make sure all sampling configurations function as expected.

# Component 2

    **Params:**

        **Inputs:** Protein pair database file or pandas dataframe formatted for pairpro. Sampling
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

## Software Component Three: hmmer.py

    **Params:** sampled sequence data from component 1

    **Inputs:** protein pair amino acid sequences (API) or their associated PID's (pyhmmer/local)

    **Outputs:** for the two proteins in a pair, you've a target family based on shared accession and a boolean to check if a pair are functional or not according to pfam parsing

    **Metrics:** Not impelemented yet

    **Packages:** pandas, biopython, pyhmmer, httpx, nest-asyncio

Component 3 aims to use the HMMER algrothim running against the pfam database on all protein pairs specified by the user. In this case, we are using the Learn2thermDB protein pairs.
This component has two options to be ran locally or using the online HMMER server API. The two options have different subcomponents, use-cases and tests. Depending on the number of sequences provided by the users, the code will either run locally or using the online HMMER server API. In the future, we aim to make this configurable by the user. The code will run the HMMER program on the user's protein pairs and return the target family and a boolean as well as Jaccard score to check if a pair are functional or not according to pfam parsing.

### API Component Three: hmmer.py

**Use case**:

        The API code comprises the asynchronous implementation of HMMER using the online server API provided by EBI. This solution submits multiple protein sequences for HMMER analysis and processes the server's responses concurrently. It handles request failures through a retry mechanism and parallelizes response processing with a process pool executor.

**Code**:
The code mainly includes four functions that are below.

        1. send_request(): Submits a protein sequence for analysis to the HMMER server.

        *Parameters*:
                semaphore: An asyncio.Semaphore instance. This semaphore controls the number of concurrent requests made to the HMMER API.
                sequence: A string representing the protein sequence to be analyzed.
                client: An instance of httpx.AsyncClient. This is the HTTP client used to send the request.
        *Return*:
                response: The httpx.Response object received from the HMMER API.
        *Exception Handling*:
                Raises an httpx.HTTPStatusError if the HTTP request returned an error status code.
                Raises an httpx.TimeoutException if the request times out.

        2. process_response(): Processes the server's responses and handles request retries. This function processes the response received from the HMMER API. If a request fails, it will attempt to retry the request up to max_retries times.

        *Parameters*:
                semaphore: An asyncio.Semaphore instance. This semaphore controls the number of concurrent
                requests made to the HMMER API.
                sequence: A string representing the protein sequence associated with the response.
                response: An httpx.Response object. This is the response received from the HMMER API.
                client: An instance of httpx.AsyncClient. This is the HTTP client used to send subsequent requests if necessary.
                pid: An integer representing the protein ID associated with the sequence.
                max_retries (optional): An integer specifying the maximum number of retries for failed requests. Defaults to 3.
        *Return*:
                A pandas DataFrame containing the search results for the protein sequence, or None if an error occurred.
        *Exception Handling*:
                Raises a KeyError if the expected key is not found in the response.
                Raises a json.JSONDecodeError if JSON decoding fails.

        3. hmmerscanner(): Submits multiple protein sequences for HMMER analysis and processes the server's responses concurrently.

        *Parameters*:
                df: A pandas DataFrame containing protein sequences.
                k: An integer specifying the number of protein sequences to search.
                max_concurrent_requests: An integer specifying the maximum number of concurrent requests
                to the HMMER API.
                output_path: A string specifying the directory where the output data will be stored.
        *Return*:
                A pandas DataFrame containing the search results for all protein sequences.
        *Exception Handling*:
                Raises a ValueError if the number of sequences exceeds the limit of 1000.

        4. run_hmmerscanner(): A wrapper function that sets up the event loop and initiates the asynchronous HMMER scanning process.

        *Parameters*:
                df: A pandas DataFrame containing protein sequences.
                k: An integer specifying the number of protein sequences to search.
                max_concurrent_requests: An integer specifying the maximum number of concurrent requests to the HMMER API.
        *Return*:
                A pandas DataFrame containing the search results for all protein sequences.
        *Exception Handling*:
                Raises a nest_asyncio.NestingError if the event loop is already running.
                Propagates any exceptions raised by the hmmerscanner() function.

        The send_request() and process_response() functions are utilized within the hmmerscanner() function, which uses asyncio tasks to manage the concurrent execution of these functions.

**Important Points to Note**:

        send_request(): This function submits a protein sequence for analysis to the HMMER server. The sequence is sent as a POST request. A semaphore is used to control the concurrency level of the requests, which avoids overloading the server.

        process_response(): This function processes the server's response and handles retries for failed requests. It retrieves the redirect URL from the response's headers, sends a GET request to the URL, and retrieves the search results from the JSON response. In case of a read timeout, the function retries the request with exponential backoff, up to a maximum number of retries. If the response contains hits, it constructs a pandas DataFrame from the hits and returns it. Otherwise, it returns None.

        hmmerscanner(): This function concurrently submits multiple protein sequences for HMMER analysis and processes the server's responses. It first checks if the number of sequences exceeds the limit of 1000, and if it does, it suggests the user to use the local function and returns an empty DataFrame. Otherwise, it creates tasks for sending requests and processing responses and gathers the results from all tasks. Finally, it constructs a DataFrame containing the search results for all protein sequences and saves it to a CSV file.

        run_hmmerscanner(): This function sets up the event loop and initiates the asynchronous HMMER scanning process. It applies nest_asyncio to allow nested usage of asyncio.run() and then runs the hmmerscanner() function in the event loop.

**Test**:

        The testing phase will involve validating input sequence data, ensuring successful creation of input files, as well as verifying the creation of HMMER output files. We will also validate the correct chunking of sequences and the functioning of the embarrassingly parallel Python technique in the HMMER running subcomponent. The development and testing of the parsing and computation of functional pairs from HMMER outputs subcomponent are ongoing.
        To use the asynchronous HMMER scanner, simply import the code and call the run_hmmerscanner() function with the appropriate arguments. For example:

        """
        df = pd.read_csv('protein_sequences.csv')
        k = 500
        max_concurrent_requests = 50
        run_hmmerscanner(df, k, max_concurrent_requests)
        This will run the HMMER scanner on the first 500 protein sequences from the input DataFrame, with a maximum of 50 concurrent requests to the HMMER server. The results will be returned as a DataFrame.
        """

**Possible Errors**:

        httpx.HTTPStatusError: If the HTTP request returned a status code that denotes an error.
        httpx.TimeoutException: If the request times out.
        KeyError: If expected key is not found in the response.
        json.JSONDecodeError: If JSON decoding fails.
        ValueError: If the number of sequences exceeds the limit of 1000.
        nest_asyncio.NestingError: If the event loop is already running.

### Running HMMER locally via pyhmmer: hmmer.py

This component has two main worker functions:

1. local_hmmer_wrapper(): This function is the main wrapper function that calls the other functions in this component. It queries the constructed database to get sequences only from chunked PID inputs. Then, it converts the query result to a DataFrame, which is then inputted to convert string sequences to pyhmmer digital blocks. Finally, it run HMMER via pyhmmer with the provided sequences, and parses the pyhmmer output and saves it to a CSV file.

2. process_pairs_table(): This function takes the output of the HMMER run and processes it to get the functional pairs. It also saves the output to a CSV file. The processing involves calculating the Jaccard score for each pair via accessions, and then filtering the pairs based on the Jaccard score and the e-value.

**Use case**:

        Utilizes the user's avaiable resources to run HMMER locally. This is useful for users who have a large number of sequences to run HMMER on, and do not want to wait for the HMMER server to process their requests. This is also useful for users who want to run HMMER on sequences that are not available on the HMMER server.

# Component 4

## Software Component Four: structures.py

    **Params:** Pandas Dataframe containing PDB or Uniprot IDs.

    **Inputs:** PDB or Uniprot IDs of proteins of interest.

    **Outputs:** csv file contaning p-values boolean (0 for False pair and 1 for True pair)

    **Metrics:** Not applicable

    **Packages:** pandas, biopython, httpx, numpy

Component 4 is responsible for obtaining protein structures and running structural alignment with FATCAT. This component first searches for experimentally solved structures using PDB IDs, and if not available, then searches for computationally generated structures using Uniprot IDs. The output of this component suggests whether two proteins share structural similarity or not.

**Use case**:
        Uses the provided PDB IDs and Uniprot IDs to compare protein structures. This allows comparison of protein at a secondary and tertiary level to complement the HMMR result which evaluates protein similarity at a primary level (amino acid sequences). 

**Test**:

# Component 5 - Training and Validation

## Software Component Five: train_val_scripts.py

    **Params:** Pandas dataframe containing sampled data from our protein database from sampling component + metric of interest from HMMER/Strucutural components..

    **Inputs:** Training: Pandas Dataframe containing known protein pairs with Boolean metric of interest.
                          Training includes feature engineering using iFeatureOmegaCLI.
                Testing: Pair of amino acid sequences from Component 3.

    **Outputs:**  Boolean Prediction (Functional pair: True or Functional pair: False)
                  Confidence score on prediction.

    **Metrics:**

    **Packages:** pandas, matplotlib, numpy, scikit-learn (.ensemble, .preprocessing, .model_selection, neighbors, .feature_selection), unittest, iFeatureOmegaCLI

The purpose of this component is to train a machine learning classifier model on our large dataset of protein pairs.

For training, the model will take in a pandas dataframe containing protein pairs as examples. The database will be cleaned until 10 quantitative features remain, including alignment metrics, optimal growth temperature, bit scores, sequence length, and alignment coverage. The sequences are then run through a feature generation algorithm, adding a variable-length vector of features to supplement model training. The cleaned dataframe will also include a Boolean target, classifying the protein as either a functional pair or not a functional pair. The model will be trained and tested based on this data.

This component will undergo significant further development to include structural information as a separate target representing Boolean protein pair functionality. The targets from HMMER and structural search are combined to improve classification accuracy.

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

### **Subcomponent 3**: Generate features with iFeatureOmega.

**Use case**:

        User needs an expanded feature space to train a machine learning classifier in order to improve classification accuracy.

**Code**:

        A set of three functions create a pandas dataframe with appended features from iFeatureOmega.

        1) Reads sequences in input pandas dataframe and generates two fasta files for each protein pair.
        2) Takes fasta files from upstream and generates descriptors using methods from iFeatureOmega. Descriptors are generated for each individual sequence, and are represented as floats.
        3) Takes newly generated descriptors and appends them into original dataframe. For each protein pair, embedded function calculates either a ratio or a difference between descriptors.

              Input: Pandas dataframe
              Output: Pandas dataframe with appended features.

**Test**:

        Asserts that input is a pandas dataframe, shape of descriptor dictionary matches number of descriptors in input, and new dataframe has indexing columns removed.

### **Subcomponent 4**: Train the model with sample data.

**Use case**:

        Model to predict eventual user input needs to be trained on the protein pair database that we currently have.

**Code**:

        Function takes in clean pandas dataframe from subcomponent 2. Splits data into dev and test set between the features and target. Trains the data with a RandomForestClassifier model from sklearn. Returns the trained model along with development and test dataframes.

          Input: Pandas dataframe
          Output: Model, Pandas dataframe

**Test**:

        Asserts that input data type is a pandas dataframe and that output has 5 objects (model, dev_X, dev_y, test_X, test_y).

### **Subcomponent 5**: Test performance of many different classifiers on training or testing data.

**Use case**:

        User has basic understanding of ML framework in Python and wants to test performance between multiple classifiers. This subcomponent is optional, and can be accessed as a separate script.

**Code**:

        Function takes cleaned dataframe with appended features. Accesses a predefined dictionary of 9 different classifier models (hyperparameters pre-optimized with Optuna) and runs all models with input data. Returns a .txt file that describes statistical results of each model, as well as a dictionary with each model's accuracy and F1 score.

              Input: Model, Pandas dataframe
              Output: Vector of predictions, dictionary of accuracy scores, evaluationResults.txt

**Test**:

        Asserts that input is a pandas dataframe, output is a numpy array, and the length of the output vector is equal to the number of examples.

### **Subcomponent 6**: Test the model with sample data.

**Use case**:

        Uses test csv to evaluate the performance of the model and return a vector of Boolean predictions. Also creates a .txt file for the user to assess performance as well as a csv with results.

**Code**:

        Function takes testing dataframe along with the model trained in the previous component. Runs testing dataframe through a RandomForestClassifier (default) and returns a vector of predictions, along with a .txt file with scoring metrics. The user also has an option to choose a different classifier model at the command line.

              Input: Model, Pandas dataframe
              Output: Numpy array, evaluationResults.txt

**Test**:

        Asserts that input is a pandas dataframe, output is a numpy array, and the length of the output vector is equal to the number of examples.


#### Plan Outline

1. Get data from component 3. This should already be in a pandas dataframe (data prep is included in C3).
2. Clean the data to prepare it for model training and testing.
3. Train and test the model, return scores, plot, and any other necessary indicator of model performance.
4. Input new user data and return a functionality score for the input protein pair.