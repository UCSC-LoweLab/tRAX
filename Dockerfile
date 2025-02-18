FROM continuumio/miniconda:latest

# Set working directory and copy TRAX software into docker container
COPY . /opt/trax/
RUN rm -rf /opt/trax/__pycache__
ENV PATH /opt/trax:$PATH

# Install Core TRAX Dependencies
RUN conda env update -n base -f /opt/trax/trax_env.yaml

# Add empty folder for RNA database docker volumes
RUN mkdir /rnadb && \
     chmod -R 777 /rnadb && \
     chmod -R 777 /home

# Add user for the container
RUN useradd -ms /bin/bash jerry
USER jerry
WORKDIR /home/jerry
