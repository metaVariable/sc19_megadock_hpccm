# MEGADOCK-HPCCM

## Description

MEGADOCK-HPCCM is a HPC container making workflow for [MEGADOCK](https://github.com/akiyamalab/MEGADOCK) application on HPC environment by using [HPCCM (HPC Container Maker)](https://github.com/NVIDIA/hpc-container-maker/) framework. It can generate container specification (recipe) files in [Dockerfile](https://docs.docker.com/engine/reference/builder/) and [Singularity definition](https://sylabs.io/guides/3.3/user-guide/definition_files.html) format from single python code.
The container has necessary GPU, OpenMPI, FFTW, InfiniBand, Intel Omni-Path libraries to perform MEGADOCK application on HPC environment. Users can add user arguments to specify the MPI library version in the container for considering the MPI library compatibility between containers and HPC environments.


## Requirements

- [NVIDIA GPU devices, drivers](https://www.nvidia.com/)
- [HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker/)
- [Docker](https://www.docker.com/) (if you use)
- [Singularity](https://sylabs.io/) (if you use)

## Repository overview
```
.
├── data                                # directory for storing input
│   └── ...                             #   small pdb data for sample docking
├── sample                              # container image recipes used in the poster's experiments
│   ├── Dockerfile                      #   Dockerfile for general environments
│   ├── singularity_ompi-2-1-2_opa.def  #   Singularity definition for TSUBAME3.0
│   └── singularity_ompi-3-1-6_ofed.def #   Singularity definition for ABCI
├── script                              # 
|   └── makeTable.sh                    # script for generating input docking list (table)
├── megadock-scfa20                     # source code of MEGADOCK application
├── megadock_hpccm.py                   # HPCCM recipe
├── Makefile                            # Makefile for image building
└── README.md                           # this document

# The directory will be generated after running scripts
.
├── table                           # directory for storing docking metadata
└── out                             # directory for storing output
```

----

## Quick Link

- [Docker environment](#Docker-environment)
- [Singularity environment](#Singularity-environment)

----

## Docker environment

### Requirements

- pip, python (for HPCCM)
- docker ( > 19.03 )
  - or `nvidia-docker` for gpu support

### Setting up: install HPCCM

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm
```

### Build image: generate Dockerfile and build

``` sh
# generate 'Dockerfile' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format docker > Dockerfile

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format docker --userarg ompi=3.1.3 fftw=3.3.8 > Dockerfile

## Available userargs:
##  ompi=${ompi_version} : version of OpenMPI library
##  fftw=${fftw_version} : version of FFTW library

# build a container image from Dockerfile
docker build . -f Dockerfile -t megadock:hpccm
```

### Run: MEGADOCK calculation with small dataset

```sh
# run with host gpus
docker run --rm -it --gpus all \
  -v `pwd`/data:/data  megadock:hpccm \
  mpirun --allow-run-as-root -n 2 /workspace/megadock-gpu-dp -tb /data/SAMPLE.table
```

### Run: MEGADOCK calculation with ZDOCK Benchmark 5.0

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 5.0)
mkdir -p data
wget https://zlab.umassmed.edu/benchmark/benchmark5.tgz
tar xvzf benchmark5.tgz -C data
rm -f benchmark5.tgz

# create docking table using script (only 100 pairs)
INTERACTIVE=1 TABLE_ITEM_MAX=100 RUNTIME_RELATIVE_ROOT=/ script/makeTable.sh . data/benchmark5/structures/ \*_r_b.pdb \*_l_b.pdb test100pairs

# Note: 
# - unset ${TABLE_ITEM_MAX} variable to unlimit the number of docking calculations (all-to-all)
# - if you need to change the repository root path when runtime, use ${RUNTIME_RELATIVE_ROOT} to modify path in generating the table.

# run
docker run --rm -it --gpus all \
  -v `pwd`/data:/data -v `pwd`/table:/table -v `pwd`/out:/out \
  megadock:hpccm \
    mpirun --allow-run-as-root -n 2 -x OMP_NUM_THREADS=20 \
      /workspace/megadock-gpu-dp -tb /table/test100pairs/test100pairs.table
```

----

## Singularity environment

### Requirements

- pip, python (for HPCCM)
- singularity
  - require `singularity exec` command on HPC system
  - require privilege for `sudo singularity build` or `singularity build --fakeroot` (>= 3.3)

Note: Following commands should be executed on your local environment where you have system privilege.

### Setting up: install HPCCM

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm
```

### Build image: generate Singularity definition file and build

``` sh
# generate 'singularity.def' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format singularity > singularity.def

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format singularity --userarg ompi=3.1.3 fftw=3.3.8 ofed=True > singularity.def

## Available userargs:
##  ompi=${ompi_version} : version of OpenMPI library
##  fftw=${fftw_version} : version of FFTW library
##  ofed=${True|False} : flag for install 'Mellanox OpenFabrics Enterprise Distribution for Linux'
##  opa=${True|False} : flag for install Intel Ompni-Path dependencies

# build a container image from Dockerfile
sudo singularity build megadock-hpccm.sif singularity.def

## or '.simg' format (singularity < 3.2)
sudo singularity build megadock-hpccm.simg singularity.def
```

### Run: MEGADOCK calculation with small dataset

- **Notes:**
  - Following commands should be running on HPC environment (compute-node with gpus).
  - Please replace `${SINGULARITY_IMAGE}` to **path to the container image file** on your environment.
  - **Please read the 'Singularity' section of system manual** which provided by your HPC system. We must add specific options for singularity runtime when using system resources.
    - e.g.) Volume option (`-B XXX`) for mounting system storage, applications, libraries, etc.

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# singularity exec 
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 /workspace/megadock-gpu-dp -tb data/SAMPLE.table
```

### Run: MEGADOCK calculation with ZDOCK Benchmark 5.0

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 5.0)
mkdir -p data
wget https://zlab.umassmed.edu/benchmark/benchmark5.tgz
tar xvzf benchmark5.tgz -C data
rm -f benchmark5.tgz

# create docking table using script (only 100 pairs)
INTERACTIVE=1 TABLE_ITEM_MAX=100 script/makeTable.sh . data/benchmark5/structures/ \*_r_b.pdb \*_l_b.pdb test100pairs

# Note: 
# - unset ${TABLE_ITEM_MAX} variable to unlimit the number of docking calculations (all-to-all)
# - if you need to change file path in compute-node, use ${RUNTIME_RELATIVE_ROOT} to modify path in generating the table.

# ${SINGULARITY_IMAGE}: path to the singularity image file

# singularity exec 
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 -x OMP_NUM_THREADS=20 \
  /workspace/megadock-gpu-dp -tb table/test100pairs/test100pairs.table

# singularity exec (with host MPI library)
mpirun -n 2 -x OMP_NUM_THREADS=20 \
  singularity exec --nv ${SINGULARITY_IMAGE} \
  /workspace/megadock-gpu-dp -tb table/test100pairs/test100pairs.table
```
