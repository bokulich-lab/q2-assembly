FROM continuumio/miniconda3:latest AS base

ARG EPOCH
ARG DISTRO
ARG ENVIRONMENT

ENV PATH=/opt/conda/envs/${DISTRO}-${EPOCH}/bin:$PATH \
    LC_ALL=C.UTF-8 LANG=C.UTF-8 \
    MPLBACKEND=agg \
    UNIFRAC_USE_GPU=N \
    HOME=/home/qiime2 \
    XDG_CONFIG_HOME=/home/qiime2 \
    ENV_NAME=${DISTRO}-${EPOCH}

WORKDIR /home/qiime2
COPY environment.yml

RUN apt-get install -y --no-install-recommends wget procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN conda update -qy conda \
    && conda install -c conda-forge -qy mamba \
    && mamba env create -n ${DISTRO}-${EPOCH} --file https://raw.githubusercontent.com/qiime2/distributions/dev/${EPOCH}/${DISTRO}/${ENVIRONMENT}/qiime2-${DISTRO}-ubuntu-latest-conda.yml \
    && mamba env update -n ${DISTRO}-${EPOCH} --file environment.yml \
    && mamba clean -all -yes \
    && chmod -R a+rwx /opt/conda

SHELL ["conda", "run", "-n", "${ENV_NAME}", "/bin/bash", "-c"]

COPY . ./plugin
RUN pip install ./plugin

RUN /bin/bash -c "source activate ${DISTRO}-${EPOCH}"
ENV CONDA_PREFIX=/opt/conda/envs/${DISTRO}-${EPOCH}/
RUN qiime dev refresh-cache
RUN echo "source activate ${DISTRO}-${EPOCH}" >> $HOME/.bashrc
RUN echo "source tab-qiime" >> $HOME/.bashrc


FROM base AS test

RUN pip install pytest pytest-cov coverage parameterized pytest-xdist
CMD ["make", "-f", "/plugin/Makefile", "test-cov"]

FROM base AS prod

# Important: let any UID modify these directories so that
# `docker run -u UID:GID` works
RUN rm -rf ./plugin
RUN chmod -R a+rwx /home/qiime2
