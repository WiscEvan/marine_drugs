FROM nextgenusfs/funannotate-slim

# install databases
RUN funannotate setup -i all -b dikarya microsporidia embryophyta metazoa arthropoda vertebrata protists

# required for nextflow
RUN apt-get update \
    && apt install -y procps g++ \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# When image is run, run the code with the environment
SHELL ["/bin/bash", "-c"]
CMD funannotate