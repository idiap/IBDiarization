# Build stage
FROM ubuntu:16.04 AS build

LABEL vendor="Idiap Research Institut"
LABEL project="SUMMA"
LABEL maintainer="william.droz@idiap.ch"

# Install dependencies and tools to build IBDiarization
RUN apt-get update -y && apt-get install -y \
  autoconf \
  automake \
  build-essential \
  cmake \
  g++ \
  git \
  libboost-all-dev \
  libomp-dev \
  libatlas-dev \
  libatlas-base-dev \
  python3 \
  sox \
  subversion \
  wget \
  zlib1g-dev

RUN git clone https://github.com/kaldi-asr/kaldi.git && cd kaldi &&  git checkout 5.4 && cd tools && make -j 2 && cd ../src && ./configure && make -j 2

RUN git clone https://github.com/wdroz/IBDiarization.git

# Build IBDiarization
ENV IB_DIARIZATION_HOME /IBDiarization
WORKDIR /IBDiarization/src/diarization/cmake/
RUN cmake . && make -j 4

# Final stage
FROM ubuntu:16.04

# Install only the dependencies to run IBDiarization
RUN apt-get update -y && apt-get install -y \
  ffmpeg \
  libatlas-dev \
  libatlas-base-dev \
  libboost-all-dev \
  libomp-dev \
  python3 \
  python3-pip \
  wget

RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
   /bin/bash ~/miniconda.sh -b -p /opt/conda && \
   rm ~/miniconda.sh && \
   ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
   echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
   echo "conda activate base" >> ~/.bashrc

RUN echo install with conda...
RUN . /opt/conda/etc/profile.d/conda.sh && conda install -c wdroz pydiarization

# Copy IBDiarization from the previous stage
COPY --from=build /IBDiarization /IBDiarization

# Copy Kaldi bins from the previous stage
COPY --from=build /kaldi/src/featbin /featbin

ENV PATH="/featbin:${PATH}"

ENV PATH="/IBDiarization/src/diarization/cmake:${PATH}"

WORKDIR /IBDiarization

# Default result directory
RUN mkdir result.dir

# Test the wrapper
RUN . /opt/conda/etc/profile.d/conda.sh && /opt/conda/bin/python3 -m pydiarization.test_diarization_wrapper

ENV LANG C.UTF-8
ENV PYTHONUNBUFFERED y

RUN /featbin/compute-mfcc-feats --help

ENTRYPOINT [ "/usr/bin/env", "python3", "-u", "/IBDiarization/rabbitmq.py" ]
