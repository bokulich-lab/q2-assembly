FROM continuumio/miniconda3:latest AS base

ARG EPOCH
ARG ENVIRONMENT
ARG PLUGIN_NAME

ENV PATH=/opt/conda/envs/${PLUGIN_NAME}-${EPOCH}/bin:$PATH \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    MPLBACKEND=agg \
    UNIFRAC_USE_GPU=N \
    HOME=/home/qiime2 \
    XDG_CONFIG_HOME=/home/qiime2 \
    ENV_NAME=${PLUGIN_NAME}-${EPOCH}

WORKDIR /home/qiime2
COPY environment.yml .

RUN apt-get install -y --no-install-recommends wget procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN conda update -qy conda \
    && conda install -c conda-forge -qy mamba \
    && mamba env create -n ${ENV_NAME} --file environment.yml \
    && mamba clean --all --yes \
    && chmod -R a+rwx /opt/conda

COPY . ./plugin
RUN mamba run -n ${ENV_NAME} pip install ./plugin

RUN /bin/bash -c "source activate ${ENV_NAME}"
ENV CONDA_PREFIX=/opt/conda/envs/${ENV_NAME}/
RUN mamba run -n ${ENV_NAME} qiime dev refresh-cache
RUN echo "source activate ${ENV_NAME}" >> $HOME/.bashrc
RUN echo "source tab-qiime" >> $HOME/.bashrc


FROM base AS test

RUN mamba run -n ${ENV_NAME} pip install pytest pytest-cov coverage parameterized pytest-xdist
CMD mamba run -n ${ENV_NAME} make -f ./plugin/Makefile test-cov

FROM base AS prod

# Important: let any UID modify these directories so that
# `docker run -u UID:GID` works
RUN rm -rf ./plugin
RUN chmod -R a+rwx /home/qiime2
