FROM continuumio/miniconda3

WORKDIR /PairProphet

# copy contents of current host directory into container
COPY . /PairProphet/

# install OS dependencies
RUN apt update && apt install -y build-essential

# Create the environment:
#COPY environment.yml .
RUN conda config --add channels conda-forge
RUN conda install mamba

# create the environment
RUN mamba env create -f environment.yml 
#&& \
    #conda init bash && \
    #echo "conda activate myenv" >> ~/.bashrc

# Make RUN commands use the new environment:
SHELL ["conda","run", "-n", "pairpro", "/bin/bash", "-c"]

ENTRYPOINT ["pip", "install" , "."] 

# # set working directory
WORKDIR /PairProphet/scripts/

CMD ["python", "download_pfam.py"]
# CMD ["python", "scripts/train_model.py", "--hmmer", "True"]
 

# CMD ["pwd"]

# ENTRYPOINT ["conda", "run", "-n", "pairpro", "python", "train_val_classification.py"]

# CMD ["python", "train_val_classification.py"] < --- for pairpro/train_val_classification.py
# CMD ["python", "train_model.py", "--blast", "True"] 

# commmand to build the image: docker build -t training_script .

# need init.py to run training script out of container
# The code to run when container is started:
#ENTRYPOINT ["python",  "-c", "train_model.py", "--blast" , "True", "--hmmer", "True", "http.server", "8000"]

# Make RUN commands use the new environment:


