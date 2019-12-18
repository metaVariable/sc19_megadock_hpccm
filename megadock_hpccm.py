"""
HPC Base image

Contents:
  CentOS 7 (default)
  CUDA version 10.0 (default)
  Mellanox OFED version 4.6-1.0.1.1 ('ofed=True')
  Intel OPA driver/library (upstream, 'opa=True')
  GNU compilers (upstream)
  FFTW version 3.3.8 (default)
  OpenMPI version 3.1.3 (default)
"""
# pylint: disable=invalid-name, undefined-variable, used-before-assignment

# userargs
base_image   = USERARG.get('base', 'nvidia/cuda:10.0-devel-centos7')
ompi_version = USERARG.get('ompi', '3.1.3')
fftw_version = USERARG.get('fftw', '3.3.8')
ofed_flag    = USERARG.get('ofed', False)
opa_flag     = USERARG.get('opa', False)

######
# Devel stage
######

# base image
devel_image = base_image

Stage0.name = 'devel'

Stage0 += comment(__doc__, reformat=False)

Stage0 += baseimage(image=devel_image, _as='devel')

# OFED
if ofed_flag:
  Stage0 += mlnx_ofed(version='4.6-1.0.1.1')

# Intel OPA
if opa_flag:
  Stage0 += packages(
    yum=['numactl-libs', 'hwloc-libs', 'libfabric', 'libibverbs', 'infinipath-psm', \
         'opa-basic-tools', 'rdma-core', 'libpsm2', \
         'libhfil', 'libibverbs-devel', 'libsysfs-devel']
  )

# MEGADOCK deps
Stage0 += packages(
    yum=['cuda-samples-10-0', 'ssh']
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

OpenMPI_with_verbs = ofed_flag or opa_flag

# OpenMPI
Stage0 += openmpi(
    version=ompi_version,
    prefix='/usr/local/openmpi',
    cuda=True, 
    infiniband=OpenMPI_with_verbs,
    configure_opts=[
        '--enable-mpi-cxx'
        ],
    toolchain=compiler.toolchain
)

# MEGADOCK
Stage0 += copy(src='./megadock-scfa20', dest='/workspace')
Stage0 += copy(
    src='./Makefile',
    dest='/workspace/Makefile'
)

Stage0 += shell(commands=['cd /workspace', 'make -j$(nproc)'])

