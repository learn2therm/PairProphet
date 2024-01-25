# need to make more robust in the future
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / "README.md").read_text(encoding="utf-8")

packages = find_packages()
print(f"Found packages {packages}")

setup(
    name="pairpro", # Required  
    version="2023.6",  # Required
    long_description=long_description, 
    long_description_content_type="text/markdown", 
    url="https://github.com/learn2therm/PairProphet",
    author="Humood Alanzi, Logan Roberts, Ryan Francis, Amin Mosallanejad, Chau Vuong ", 
    author_email="halanzi@uw.edu",
    classifiers=[
        # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 1 - Planning",
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: bioinformatics',
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    keywords="proteins, function",
    packages=packages,
    install_requires=[
        'Click',
    ],
    python_requires=">=3.10, <4",
    extras_require={ 
        "test": ["coverage", 'pytest'],
    },
    package_data={},
    data_files=None, 
    entry_points={},
)
