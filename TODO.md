# Future work for PairProhpet
PairProphet is a software in activate development. We are working on the following features:

## Overall improvements
- [ ] Add more visualization tools for the differents components of the model.
    - [ ] Pfam MSA visualizations
    - [ ] Structure visualizations
- [ ] Adopt an object-oriented approach to the code to allow the user more configuration options.
    - [ ] Add more options to the command line interface.
- [ ] Imlpemenet data version control to allow the user to track the changes in the data and improve model reproducibility.
- [ ] Use the OMA table of orthologs as the main input for the model instead of learn2thermDB.
- [ ] Explore Guassian Processes as an alternative to the current model.
    - [ ] (optional) Explore deep learning models as an alternative to the current model.

## Preprocessing improvements
- [ ] optimize the local alignment step to allow the user to run the model in a cluster.

## Pfam improvements
- [ ] Utilize better HPC techniques to allow the user to run the model in a cluster.

## structure improvements
- [ ] create a python wrapper for FATCAT to allow the user to run the structural alignment from the command line to be more accessible.

## API improvements (optional)
- [ ] Create a REST API to allow the user to run the model from a web interface.