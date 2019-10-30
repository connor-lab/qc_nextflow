FROM nfcore/base
LABEL authors="Matt Bull, Amy Gaskin" \
      description="Docker image containing all requirements for connor-lab/nextflow_qc pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get install -y make
ENV PATH /opt/conda/envs/connor-lab-nextflow_qc-0.1/bin:$PATH
