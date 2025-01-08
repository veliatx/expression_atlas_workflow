FROM ubuntu:20.04
LABEL mantainer="Mitch Ledwith (mitch@veliatx.com)"

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        pigz \
        seqtk \
        wget

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.10/sratoolkit.3.0.10-ubuntu64.tar.gz --no-check-certificate && \
    tar -xzvf sratoolkit.3.0.10-ubuntu64.tar.gz

ENV PATH="${PATH}:/sratoolkit.3.0.10-ubuntu64/bin:"
