# 
# HPC Base image
# 
# Contents:
#   CUDA version 10.0
#   OPA/infiniband Basic Driver (upstream)
#   GNU compilers (upstream)
#   FFTW version 3.3.8 (default)
#   OpenMPI version 3.1.3 (default)
# 

FROM nvidia/cuda:10.0-devel-centos7 AS devel

RUN yum install -y \
        cuda-samples-10-0 \
        hwloc-libs \
        infinipath-psm \
        libfabric \
        libhfil \
        libibverbs \
        libibverbs-devel \
        libpsm2 \
        libsysfs-devel \
        numactl-libs \
        opa-basic-tools \
        rdma-core \
        ssh && \
    rm -rf /var/cache/yum/*

# GNU compiler
RUN yum install -y \
        gcc \
        gcc-c++ \
        gcc-gfortran && \
    rm -rf /var/cache/yum/*

# FFTW version 3.3.8
RUN yum install -y \
        file \
        make \
        wget && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp ftp://ftp.fftw.org/pub/fftw/fftw-3.3.8.tar.gz && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/fftw-3.3.8.tar.gz -C /var/tmp -z && \
    cd /var/tmp/fftw-3.3.8 &&  CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/fftw --enable-float --enable-sse2 && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/fftw-3.3.8.tar.gz /var/tmp/fftw-3.3.8
ENV LD_LIBRARY_PATH=/usr/local/fftw/lib:$LD_LIBRARY_PATH

# OpenMPI version 3.1.3
RUN yum install -y \
        bzip2 \
        file \
        hwloc \
        make \
        numactl-devel \
        openssh-clients \
        perl \
        tar \
        wget && \
    rm -rf /var/cache/yum/*
RUN mkdir -p /var/tmp && wget -q -nc --no-check-certificate -P /var/tmp https://www.open-mpi.org/software/ompi/v3.1/downloads/openmpi-3.1.3.tar.bz2 && \
    mkdir -p /var/tmp && tar -x -f /var/tmp/openmpi-3.1.3.tar.bz2 -C /var/tmp -j && \
    cd /var/tmp/openmpi-3.1.3 &&  CC=gcc CXX=g++ F77=gfortran F90=gfortran FC=gfortran ./configure --prefix=/usr/local/openmpi --enable-mpi-cxx --with-cuda --with-verbs && \
    make -j$(nproc) && \
    make -j$(nproc) install && \
    rm -rf /var/tmp/openmpi-3.1.3.tar.bz2 /var/tmp/openmpi-3.1.3
ENV LD_LIBRARY_PATH=/usr/local/openmpi/lib:$LD_LIBRARY_PATH \
    PATH=/usr/local/openmpi/bin:$PATH

COPY ./megadock-5.0-alpha-706cb91 /workspace

COPY ./Makefile /workspace/Makefile

RUN cd /workspace && \
    make -j$(nproc)


