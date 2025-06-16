FROM continuumio/miniconda3:latest AS base

ARG EPOCH
ARG DISTRO
ARG ENVIRONMENT

ENV PATH=/opt/conda/envs/${DISTRO}-${EPOCH}/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV MPLBACKEND=agg
ENV UNIFRAC_USE_GPU=N
ENV HOME=/home/qiime2
ENV XDG_CONFIG_HOME=/home/qiime2

WORKDIR /home/qiime2

COPY environment.yml environment.yml

RUN conda update -q -y conda \
    && conda install -c conda-forge -q -y wget mamba \
    && apt-get install -y procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget -O ${ENVIRONMENT}-env.yml https://raw.githubusercontent.com/qiime2/distributions/dev/${EPOCH}/${DISTRO}/${ENVIRONMENT}/qiime2-${DISTRO}-ubuntu-latest-conda.yml

RUN mamba env create -n ${DISTRO}-${EPOCH} --file ${ENVIRONMENT}-env.yml \
    && mamba env update -n ${DISTRO}-${EPOCH} --file environment.yml \
    && mamba clean -a -y \
    && chmod -R a+rwx /opt/conda \
    && rm ${ENVIRONMENT}-env.yml

COPY . ./plugin

RUN mamba run -n ${DISTRO}-${EPOCH} pip install ./plugin

RUN /bin/bash -c "source activate ${DISTRO}-${EPOCH}"
ENV CONDA_PREFIX=/opt/conda/envs/${DISTRO}-${EPOCH}/
RUN qiime dev refresh-cache
RUN echo "source activate ${DISTRO}-${EPOCH}" >> $HOME/.bashrc
RUN echo "source tab-qiime" >> $HOME/.bashrc


FROM base AS test

ENV ENV_NAME=${DISTRO}-${EPOCH}

RUN mamba run -n ${ENV_NAME} pip install pytest pytest-cov coverage parameterized pytest-xdist
CMD ["sh", "-c", "mamba run -n ${ENV_NAME} make test-cov"]

FROM base AS prod

# Important: let any UID modify these directories so that
# `docker run -u UID:GID` works
RUN rm -rf ./plugin
RUN chmod -R a+rwx /home/qiime2
