"""
HPC Base image

Contents:
  CUDA version 10.0
  OPA/infiniband Basic Driver (upstream)
  GNU compilers (upstream)
  FFTW version 3.3.8 (default)
  OpenMPI version 3.1.3 (default)
"""
# pylint: disable=invalid-name, undefined-variable, used-before-assignment

# base image
devel_image = 'nvidia/cuda:10.0-devel-centos7'

# library version
ompi_version = USERARG.get('ompi', '3.1.3')
fftw_version = USERARG.get('fftw', '3.3.8')

######
# Devel stage
######

Stage0.name = 'devel'

Stage0 += comment(__doc__, reformat=False)

Stage0 += baseimage(image=devel_image, _as='devel')

Stage0 += packages(
    yum=['numactl-libs', 'hwloc-libs', 'libfabric', 'libibverbs', 'infinipath-psm', \
         'opa-basic-tools', 'rdma-core', 'libpsm2', \
         'libhfil', 'libibverbs-devel', 'libsysfs-devel', \
         'cuda-samples-10-0', 'ssh']
)

# GNU compilers
compiler = gnu()
Stage0 += compiler

# FFTW
Stage0 += fftw(
    version=fftw_version, 
    prefix='/usr/local/fftw',
    configure_opts=[
            '--enable-float',
            '--enable-sse2'
            ],
    toolchain=compiler.toolchain
)

# OpenMPI
Stage0 += openmpi(
    version=ompi_version,
    prefix='/usr/local/openmpi',
    cuda=True, 
    infiniband=True,
    configure_opts=[
        '--enable-mpi-cxx'
        ],
    toolchain=compiler.toolchain
)

# MEGADOCK
Stage0 += copy(src='./megadock-5.0-alpha-706cb91', dest='/workspace')
Stage0 += copy(
    src='./Makefile',
    dest='/workspace/Makefile'
)

Stage0 += shell(commands=['cd /workspace', 'make -j$(nproc)'])
