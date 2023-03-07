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

        **Outputs:** Augmented database file containing relevant proteins and pairs for ValidProt. Keeps original tables intact.
                     It is recommended to have a backup copy of the original database as well as 100 GB of free storage and at 
                     least 30 GB of available memory. 

    Component 0 is specifically concerned with the learn2therm database and is not intended for use by users with other data files.
    Additionally, this code need only be executed once to generate all files necessary for downstream processing. Users can input
    their own similarly formatted data for classification starting in component 1. Component 0 can feed both component 1 and 
    component 2. Additionally, users have the option to export intermediate metadata on the flow of information between learn2therm
    and ValidProt as Sankey plots.

    **Packages:** os, sys, time, pandas, numpy, duckdb, plotly, kaleidoscope

# Component 1

    **Params:** 

        **Inputs:** Protein pair database file or pandas dataframe formatted according to instructions in README.
    
                    Input file should contain at minimum:

                        1. AA sequence for both members of each pair.
                        2. Clear indexing/column names to identify which belongs in meso-thermo groups.
                        3. Most recent metrics for model training (e.g. alignment/coverage scores). Included in docs
                           and c1 docstrings.

        **Outputs:** Pandas dataframe(s) containing data sample of specified size.

    Component 1 formats data for the ValidProt model. The resulting dataframe is pruned to contain only feature and 
    target data.

    **Packages:** os, sys, pandas, numpy, duckdb
    
### **Subcomponent 1**: 

**Use case**: User supplies path to database or .csv file and desired sample size. Pulls data into DataFrame using one of several sampling
              methods. 

**Test**: Test assert that DataFrame was generated, compare shape with desired result. Check for obvious errors such as inconsistent data type within columns. Make sure all sampling configurations function as expected.

### **Subcomponent 2**: 

**Use case**: User has not supplied E-values for protein pairs, computes on-the-fly for Component 2 processing.

    Subcomponent 2 is under development and will not be release with the initial version of ValidProt. Users will
    need to format and clean their own data.
          

# Component 2

    **Params:** 

        **Inputs:** Protein pair database file or pandas dataframe formatted according to instructions in README. Sampling 
                    parameters for training set.

        **Outputs:** DataFrame containing training data. Exported analytics and plots for user reference (optional).
                     parameters.

    Component 2 is a sampler for training the learn2therm database. Users can select the desired sampling method, and the function
    can seek data subsets satisfying those parameters. For example, rare features can be oversampled or common ones undersampled 
    to give the model broader applicability across protein pairs.

    **Packages:** pandas, numpy, scipy.stats, matplotlib.pyplot, seaborn, math
    
### **Subcomponent 1**: 

**Use case**: Basic statistics and plots associated with Data funneled from Component 1 are exported for user. Includes PDF and other relevant distributions.

**Test**: (If relevant) Test that plots were generated and statistical variables exist.

### **Subcomponent 2**: 

**Use case**: Decides if sample size is sufficient for downstream processing. Re-samples without replacement if needed. Return to Component 1.1

**Test**: (If relevant) Test that sample is representative of mean or matches parameters supplied by user (e.g. specific E-values, etc.)

### **Subcomponent 3**: 

**Use case**: Compiled sample set is sufficient for downstream processing. Collect sequences into DataFrame for pfam analysis. Other information stored in another indexed DataFrame.

**Test**: DataFrame shape matches expected value. Software acknowledges optional parameters from original data input.

# Component 3
## Software Component Three: c3.0_pfam.py
    
    **Params:** sequence data from component 1 and 2

    **Inputs:** protein pair amino acid sequences

    **Outputs:** family identification, pfam E value, and identity scores

    **Metrics:**

    **Packages:** pandas, unittest, biopython, HMMER, pfam database (database, optional)

### **Subcomponent 1**: Send sequence and get result

**Use case**: User sends their AA sequence of interest and gets a basic pfam result

Currently, we are doing two approaches for this:
1. Using interpro's API, and sending HTTPs request, we query pfam
    - pro: minimal amount of packages needed
    - issue: depedent on an online server and potentially slow
2. Downloading HMMER and pfam and running things locally
    - pro: quick
    - issue: requires more packages from users (hence we are thinking of making this optional)

**Test**: check if you get an empty result/call, catch as an exception

```py
assert pfam_in == [], "empty output"
```

### **Subcomponent 2**: parse and filter results

**Use case**: User obtains relevant outputs from pfam

**Test**: work-in-progress

```py
# pseudo-code
def check_result(dataframe):
    
      if no family info in dataframe[]:
          raise Exception "The pair do not have a family or they do not have the same family"
      else:
          pass
          
      if no E value in dataframe[]:
          raise ValueError
      else:
          pass
```

### **Subcompoent 3**: apply pfam to all pairs

**Use case**: User applies pfam to all their desired pair data

**Test**: work-in-progress

### **Subcomponent 4**: compute metric of interst of all pairs

**Use case**: User obtain all relevant metrics of interst for all pairs

**Test**: context-dependent + work-in-progress

#### Plan Outline
1. Have the two approaches for the user to query pfam in Python script code
2. for the first approach, figure out how to do efficent HTTPs requests and get the correct information
3. for the second approach, synergize HMMER and pfam locally
4. parse, filter, and compute metrics of interest for a large dataset

# Component 4
## Work-in-progress
This component won't be done this quarter as we will only ensure that the pfam functionality works first, so there will be no collecting. Once 3D structure prediction comes in, we formally work on this.

# Component 5
## Software Component Five: c5.0_relation.py

    **Params:** 

    **Inputs:** Pandas Dataframe containing Pfam return data. Includes quantitative features (ID, some metric of percent similarlity) and string of amino acid sequence.

    **Outputs:** Quantitative functional similarlity metric.

    **Metrics:**

    **Packages:** pandas, numpy, scipy, seaborn, fuzzywuzzy, unittest

### **Subcomponent 1**: Test for pandas dataframe input (**ALREADY TESTED WITH CODE**)

**Use case**: User takes data from component 4 (where data is processed into pandas dataframe) and wants to pass it into relationship component.

```py

def check_input_type(dataframe):
    tests that input data is a pandas dataframe with assert statement. 
    assert "pandas.core.frame.DataFrame" in str(type(dataframe)) 
    Output should pass unless assert statement fails.

```

**Test**: N/A

### **Subcomponent 2**: Checks that input data is cleaned property (does it have all of the features we need, and are the features we don't need removed).

**Use case**: Input data does not include local E value, which we need as an input to our model.

(**NOT TESTED WITH CODE**)
```py
def check_input_strings(dataframe):
    
      if 'badstring' in dataframe[]:
          dataframe = dataframe.drop(dataframe['string']
      else:
          pass
          
      if 'goodstring' not in dataframe[]:
          raise KeyError
      else:
          pass
```

**Test**: 
1) 
```py
    import unittest
    class TestMissingStrings(unittest.TestCase):
        def test_missing_strings(self):
            try:
                check_input_strings(dataframe)
                self.assertTrue(False)
            except ValueError:
                self.assertTrue(True)
```

(**ALREADY TESTED WITH CODE**)
```py
def check_input_NANs(dataframe):
   #Clean out NAN's
   
   #figure out if dataframe has NaN's
    has_nan = dataframe.isna().any().any()
    nan_rows = dataframe[dataframe.isna().any(axis=1)]

    if has_nan:
        print('Dataframe has {} rows with NaN values!'.format(len(nan_rows)))
    else:
        print("DataFrame does not have any NaN values.")
        
   #Drop rows with NaN's
    dataframe = dataframe.dropna()
    print('Dataframe now has {} rows.'.format(len(dataframe)))
    
    return dataframe
    
    Output: returns dataframe, value count
```

(**NOT TESTED WITH CODE**)
```py
def verify_protein_pair(dataframe, sequence1=dataframe['meso_seq], sequence2=dataframe['thermo_seq']):
    
    Checks that input data has two protein sequences. Will need to generalize this function other data sets 
    to simply make sure two sequences are entered. Code below is for our protein database
    
    assert (dataframe['meso_seq'] in dataframe, 'Dataframe missing mesophillic sequence!')
    assert (dataframe['thermo_seq'] in dataframe, 'Dataframe missing thermophillic sequence!')
    
        
    if len(sequence1) =/ len(sequence2):
        raise ValueError
    else:
        pass
        
    Output: Nothing if test passes
   
```

**Test**: 
1) 
```py
    import unittest
    class TestProteinPairs(unittest.TestCase):
        def test_protein_pair(self):
            try:
                verify_protein_pairs(dataframe)
                self.assertTrue(False)
            except ValueError:
                self.assertTrue(True)
```

### **Subcomponent 3**: Train the model with sample data. (**NOT TESTED WITH CODE**)

**Use case**:
```py
def train_model(dataframe):
    import scipy, numpy
    Split data into dev and test (0.8/0.2 for now)
    Train model (KNN Linear Regression for now)
    Output: Print('Training successful!')
```

**Test**: 

1)
```py 
assert len(dataframe)*0.8 == len(dev_data)
```

### **Subcomponent 4**: Test the model with sample data. (**NOT TESTED WITH CODE**)

**Use case**:
```py
def test_model(dataframe):
    Runs data through model (linear regression (KNN?)
    Output: Returns model_score, confusion matrix, MSE
```
**Test**: Work-in-progress

### **Subcomponent 5**: Run confidence test on model output. (**NOT TESTED WITH CODE**)

**Use case:**

```py
def check_model_confidence(model_score, ci_data):
    Runs a statistical test on model output and compares it to sample
    Output: Returns a confidence score along with the model score
```

**Test**: 
1) Run confidence test on some data for which we know the confidence score assert that the score is correct using numpy.isclose( )

### **Subcomponent 6**: Calculate a 'functionality' metric that is the ultimate output of component five. This will factor in information from multiple software, not just Pfam. This will be built during spring quarter. (**NOT TESTED WITH CODE**)

**Use case:** We need to test that our protein pairs have a near maximal functionality score! This can be used as a basis for eventual user input scores.

```py
def calculate_functionality(model_score, dataframe):
    runs user input data through some mathematical manipulation of their model score and input data
    Output: returns a functionality score, print statement categorizing functionality score
```

**Test**: work-in-progress

#### Plan Outline
1. Get data from component 4. This should already be in a pandas dataframe (data prep is included in C4)
2. Clean the data to prepare it for model training and testing
3. Train and test the model, return scores, MSE, and any other necessary indicator of model performance
4. Run confidence test on model output to determine quality of output
5. Input new user data and return a functionality score for the input protein pair

# Component 6
## Work-in-progress
We do not anticpate to reach this component this Winter quarter, but the Spring quarter! 

