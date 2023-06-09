"""
Unit tests for the structures module.
"""
import os
import unittest
import pandas as pd
from pairpro import structures

class TestStructures(unittest.TestCase):
    """
    Test multiple functions in the structures module, including downloading structures and running strutural alignment.
    """
    def setUp(self):
        """
        Set up the test data and directories.
        """
        # Set up any necessary test data
        sample_data_path = './sample_df.csv'
        self.df = pd.read_csv(sample_data_path)
        self.pdb_dir = './pdb_files'
        self.results_file = 'structure_results.csv'

    def test_download_aff(self):
        """
        Providing a valid URL and a filename, the function should download the file and return True.
        """
        # Test the download_aff function
        session = structures.httpx.AsyncClient()
        url = 'https://alphafold.ebi.ac.uk/files/AF-Q5VSL9-F1-model_v4.pdb'
        filename = 'Q5VSL9.pdb'
        success = structures.download_aff(session, url, filename)
        self.assertTrue(success)

    def test_download_af(self):
        """
        Given a row of data, the function should download the pdb file and return True.
        The pdb file should be downloaded in the pdb_dir.
        """
        # Test the download_af function
        row = self.df.iloc[0]
        success = structures.download_af(row, 'meso_pid', self.pdb_dir)
        self.assertTrue(success)
        self.assertTrue(os.path.exists(self.pdb_dir))

    def test_run_download_af_all(self):
        """
        Given a dataframe, the function should download all pdb files and return True.
        """
        # Test the run_download_af_all function
        structures.run_download_af_all(
            self.df, 'meso_pdb', 'meso_pid', self.pdb_dir)
        self.assertTrue(os.path.exists(self.pdb_dir))

    def test_download_structure(self):
        """
        Given a dataframe, the function should download all pdb files and return True.
        Files should be downloaded in the pdb_dir. The pdb_dir should not be empty.
        """
        # Test the download_structure function, pdb files should be downloaded
        # in "pdb_files" directory
        structures.download_structure(
            self.df, 'meso_pdb', 'meso_pid', self.pdb_dir)
        # Check that the pdb_dir is not empty
        files = os.listdir(self.pdb_dir)
        self.assertNotEqual(len(files), 0)

    def test_run_fatcat_dict_job(self):
        """
        Given a dataframe, the function should run FATCAT for all pairs of proteins and return True.
        The output file is in the same directory.
        """
        # "structure_results.csv" is the output file, should be empty because we didn't download structures for thermo species
        structures.run_fatcat_dict_job(
            self.df, self.pdb_dir, self.results_file)
        self.assertTrue(os.path.exists(self.results_file))
