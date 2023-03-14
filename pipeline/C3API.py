import pandas as pd
import requests
import urllib.parse
import time

def hmmerscanner(df: pd.DataFrame):
    """
    ----------------
    Parameters
    ----------------
    df : 
        pandas.core.Dataframe
        a dataframe that has string amino acid sequences.
        
    -----------------    
    Input Data
    -----------------
    We can change them in the function.
    
    1- Meso Sequence:
        `sequence = df['meso_seq'][i]`
        We want to check the Meso sequence and Thermo Sequence separately according to the Humood 
        request, this function is for Meso Sequence for the first step, and we can change 
        that for Thermo sequence.
        And for thermo Sequence we can put `sequence = df['thermo_seq'][i]` or for pairing we can use
        `sequence = df['meso_seq'][i] + df['thermo_seq'][i]`in the function.
        
    2- Range : 
        We can change the range according to our request in the future. Now, for checking the function,
        there is the `range(0,10)` but in the future we can use `range(len(df))`.
        
    ----------------    
    Method
    ----------------
        This method used the API system and JASON to get the information from the HMMER website by 
        sending the sequences as input to scan them and get the information for them.
        
    ----------------    
    Raises
    ----------------
    Exception :
        Raises an exception if the status is pending for too long or 
        if the internet isn't working.
        And if the URL system doesn't wholly answer, it shows the error. 
        
    ---------------    
    Results
    ---------------
        All the families are in the rows, and we have many columns that show
        the information that we need in the future. We can drop some columns
        and keep the needed information.
        
    """
    
    # create an empty DataFrame to store the results
    results_df = pd.DataFrame()
    

    #The loop for checking the seuence - We can change the range for example `range(len(df))`:
    for i in range(0,10): 
        
        # This is for meso_seq; we can change that in the future according to our request. 
        sequence = df['meso_seq'][i]
        
        #url part to request to that and getting response
        url = 'https://www.ebi.ac.uk/Tools/hmmer/search/hmmscan'
        headers = {'Content-Type': 'application/x-www-form-urlencoded', 'Accept': 'application/json'}
        data = {'hmmdb': 'pfam','seq': f'>seq\n{sequence}'}
        data = urllib.parse.urlencode(data).encode('ascii')
        response = requests.post(url, headers=headers, data=data, allow_redirects=False)
        redirect_url = response.headers.get('Location')


        if redirect_url is None:
            # If the server doen't work, show this error.
            print("Error: No redirect URL found in response.")
            
          
        elif redirect_url == 'late':
            #Raises an exception if the status is pending for too long 
            response.raise_for_status()
            time.sleep(180)
            raise IOError("Error notice after 3 minutes.")
                          
        else:
            response2 = requests.get(redirect_url, headers=headers)
            #Put the results in the empty DataFrame 
            results = response2.json()
            hits = results['results']['hits']
            dfff = pd.json_normalize(hits, 'domains', ['acc', 'name', 'score', 'evalue', 'pvalue', 'desc'])
            dfff.insert(0, 'sequence', sequence)
            dfff = dfff.set_index('sequence')
            results_df = pd.concat([results_df, dfff])
            
             
            if redirect_url == 'late':
                # Raises an exception if the status is pending for too long
                response2.raise_for_status()
                time.sleep(180)
                raise IOError("Error notice after 3 minutes.")
            
    return results_df

