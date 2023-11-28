FROM continuumio/miniconda3

WORKDIR /PairProphet

# copy contents of current host directory into container
COPY . /PairProphet/

# Create the environment:
#COPY environment.yml .
RUN conda config --add channels conda-forge
RUN conda config --add channels defaults
RUN conda env create -f environment.yml

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "pairpro", "/bin/bash", "-c"]

# Demonstrate the environment is activated:
#RUN echo "Make sure flask is installed:"
#RUN python -c "import flask"

# The code to run when container is started:
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "pairpro", "python",  "-c", "import sys; print(sys.path)"]
# "testing/dev-train_model.py", "--blast" , "True", "--hmmer", "True", "http.server", "8000"]





#--------------------------------------------------------------------------------------------------

# # For more information, please refer to https://aka.ms/vscode-docker-python
# # build off alpine image
# # FROM python:3.9-alpine
# FROM continuumio/miniconda3

# # Keeps Python from generating .pyc files in the container
# ENV PYTHONDONTWRITEBYTECODE=1

# # Turns off buffering for easier container logging
# ENV PYTHONUNBUFFERED=1

# # Install pip requirements
# # COPY requirements.txt .
# WORKDIR /
# COPY . .
# # RUN python -m pip install -r requirements.txt
# # RUN conda env create -f environment.yml
# RUN echo "source activate env" > ~/.bashrc
# ENV PATH /opt/conda/envs/env/bin:$PATH

# # Docker file should exist in same directory as script.
# # can either put file outside of scripts/ and pairpro/ or inside of pairpro/,
# # and move training and input scripts inside of pairpro/.
# WORKDIR /scripts
# COPY . /scripts


# # Creates a non-root user with an explicit UID and adds permission to access the /app folder
# # For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
# #RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
# #USER appuser

# # During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
# CMD ["python", "train_model.py"]

# current building command:
# docker build -t training_script -f Dockerfile scripts/