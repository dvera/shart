FROM ubuntu:xenial

ENV PATH=/opt/FastQC:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

WORKDIR /root

RUN apt-get update && apt-get install -y \
  gcc \
  g++ \
  wget \
  make \
  zlib1g-dev \
  python-dev \
  python-pip \
  libbz2-dev \
  liblzma-dev \
  unzip \
  r-base-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  openjdk-9-jre \
  perl
