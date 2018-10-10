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

RUN git clone https://github.com/kaldi-asr/kaldi.git && cd kaldi &&  git checkout 5.4 && cd tools && make -j 4 && cd ../src && ./configure && make -j 4

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
  python3-pip

RUN pip3 install --upgrade pip 
RUN pip3 install \
  webrtcvad \
  aio-pika==0.21.0 \
  pydiarization

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
RUN python3 -m pydiarization.test_diarization_wrapper

ENV LANG C.UTF-8
ENV PYTHONUNBUFFERED y

RUN /featbin/compute-mfcc-feats --help

ENTRYPOINT [ "/usr/bin/env", "python3", "-u", "/IBDiarization/rabbitmq.py" ]
