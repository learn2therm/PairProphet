# need to make more robust in the future
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

setup(
    name="FAFSA",  
    version="2023.3",  # Required
    long_description=long_description, 
    long_description_content_type="text/markdown", 
    url="https://github.com/learn2therm/ValidProt",
    author="Humood Alanzi, Logan Roberts, Ryan Francis, Amin Mosallanejad, Chau Vuong ", 
    author_email="halanzi@uw.edu",
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="proteins, function",
    packages=find_packages(),
    python_requires=">=3.7, <4",
    extras_require={ 
        "test": ["coverage", 'pytest'],
    },

    package_data={},
    data_files=None, 
    entry_points={},
)