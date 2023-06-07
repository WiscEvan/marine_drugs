FROM biobakery/humann:latest
LABEL maintainer="erees@wisc.edu"

RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -f -p $HOME/miniconda \
    && $HOME/miniconda/bin/conda install -c conda-forge -y awscli \
    && rm Miniconda3-latest-Linux-x86_64.sh

RUN rm -rf /usr/local/lib/python3.6/dist-packages/humann/data/{chocophlan_DEMO,uniref_DEMO}

ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
WORKDIR /scratch
ENTRYPOINT [ "/usr/bin/tini", "--" ]
CMD [ "/bin/bash" ]