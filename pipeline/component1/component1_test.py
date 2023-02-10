import unittest
import component1
from input import path

# Define a class in which the tests will run
class testComp1(unittest.TestCase):

    # Each method in the class to execute a test
    def test_extension(self):
        self.assertIn(path[-4:], ['.csv','.tsv'])

if __name__ == '__main__':
    unittest.main()
