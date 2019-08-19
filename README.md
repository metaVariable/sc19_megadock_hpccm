# MEGADOCK-HPCCM

## Description

MEGADOCK-HPCCM is a HPC container making workflow for [MEGADOCK](https://github.com/akiyamalab/MEGADOCK) application on HPC environment by using [HPCCM (HPC Container Maker)](https://github.com/NVIDIA/hpc-container-maker/) framework. It can generate container specification (recipe) files in [Dockerfile](https://docs.docker.com/engine/reference/builder/) and [Singularity definition](https://sylabs.io/guides/3.3/user-guide/definition_files.html) format from single python code.
The container has necessary GPU, OpenMPI, FFTW, InfiniBand, Intel Omni-Path libraries to perform MEGADOCK application on HPC environment. Users can add user arguments to specify the MPI library version in the container for considering the MPI library compatibility between containers and HPC environments.


## Requirements

- [HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker/)
- NVIDIA GPU Device
- [Docker](https://www.docker.com/) (if you use)
- [Singularity](https://sylabs.io/) (if you use)

## Repository overview
```
.
├── data                            # directory for storing input
│   └── ...                         # 
├── sample                          # container image recipes used in the poster's experiments
│   ├── Dockerfile                  #   Dockerfile for general environments
│   ├── singularity_ompi-2-1-3.def  #   Singularity definition for TSUBAME3.0
│   └── singularity_ompi-3-1-3.def  #   Singularity definition for ABCI
├── script                          # 
|   └── makeTable.sh                # script for generating input docking list (table)
├── megadock-5.0-alpha-706cb91      # source code of MEGADOCK application
├── megadock_hpccm.py               # HPCCM recipe
├── Makefile                        # Makefile for image building
└── README.md                       # this document

# The directory will be generated after running scripts
.
├── table                           # directory for storing docking artifacts
└── out                             # directory for storing output
```

## For Singularity environment

**Note: Singularity environment may depend on host libralies, please read the manual on your HPC environment**

- requirements
  - pip, python (for HPCCM)
  - singularity
    - require privilege for `sudo singularity build`

### Installation and setup ( on local environment )

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm
```

### Generate Singularity definition, build Singularity image ( on local environment )
``` sh
# generate 'singularity.def' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format singularity > singularity.def

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format singularity --userarg ompi=3.1.3 fftw=3.3.8 > singularity.def

# build a container image from Dockerfile
sudo singularity build megadock-hpccm.sif singularity.def

# or '.simg' format (singularity < 3.2)
sudo singularity build megadock-hpccm.simg singularity.def

# please copy singularity image to HPC system on yourself ('megadock-hpccm.sif' or 'megadock-hpccm.simg')
```

### Setup and run Singularity container (on HPC environment)

- **Notes:**
  - **Following commands should be running on compute-node.** 
  - Please replace `${SINGULARITY_IMAGE}` to **path to the container image file** on your environment.
  - Please read the system manual about Singularity on your HPC system. We must add options for singularity runtime in general case.
    - e.g.) Volume option (`-B XXX`) for mounting system storage, applications, libraries, etc.

#### Test MEGADOCK calculation with small dataset

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# singularity exec 
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 /workspace/megadock-gpu-dp -tb data/SAMPLE.table
```

#### Run MEGADOCK calculation with ZDOCK Benchmark 1.0

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 1.0)
wget http://zlab.umassmed.edu/zdock/benchmark1.0.tar.gz
tar xvzf benchmark1.0.tar.gz -C data
rm -f benchmark1.0.tar.gz

# generate input docking table for MEGADOCK calculation (all-to-all dockings for ZDOCK benchmark 1.0)
INTERACTIVE=1 TABLE_ITEM_MAX=200 TSV_SIZE=50 \
script/makeTable.sh . data/benchmark1.0/unbound_pdb \*_r.pdb \*_l.pdb 200pairs
# deleting 'TABLE_ITEM_MAX' to prepare all-to-all docking pairs (3481)
# if you need to change file path in compute-node, use 'RUNTIME_RELATIVE_ROOT' to modify PATH in the table.

# singularity exec 
# Note: please replace ${SINGULARITY_IMAGE} to your path to the container image file
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 -x OMP_NUM_THREADS=$(nproc) \
  /workspace/megadock-gpu-dp -tb table/200pairs/200pairs.table

# singularity exec (with host MPI library)
# Note: please read carefully the system manual on your HPC system
mpirun -n 2 -x OMP_NUM_THREADS=$(nproc) \
  singularity exec --nv ${SINGULARITY_IMAGE} \
  /workspace/megadock-gpu-dp -tb table/200pairs/200pairs.table
```

----

## For Docker environment

- requirements
  - pip, python (for HPCCM)
  - docker > 19.03

### Installation and setup

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm
```

### Generate Dockerfile, build Docker image
``` sh
# generate 'Dockerfile' from hpccm recipe
hpccm --recipe megadock_hpccm.py --format docker > Dockerfile

## or adding 'userarg' for specifying library versions
hpccm --recipe megadock_hpccm.py --format docker --userarg ompi=3.1.3 fftw=3.3.8 > Dockerfile

# build a container image from Dockerfile
docker build . -f Dockerfile -t megadock:hpccm
```

### Run Docker container on HPC environment

#### Test MEGADOCK calculation with small dataset

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# run 
docker run --rm -it --gpus all \
  -v `pwd`/data:/data  megadock:hpccm \
  mpirun --allow-run-as-root -n 2 /workspace/megadock-gpu-dp -tb /data/SAMPLE.table
```

#### Run MEGADOCK calculation with ZDOCK Benchmark 1.0

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 1.0)
mkdir -p data
wget http://zlab.umassmed.edu/zdock/benchmark1.0.tar.gz
tar xvzf benchmark1.0.tar.gz -C data
rm -f benchmark1.0.tar.gz

# create docking table using script
INTERACTIVE=1 TABLE_ITEM_MAX=200 TSV_SIZE=50 RUNTIME_RELATIVE_ROOT=/ \
script/makeTable.sh . data/benchmark1.0/unbound_pdb \*_r.pdb \*_l.pdb 200pairs

# deleting 'TABLE_ITEM_MAX' to prepare all-to-all docking pairs (3481)

# run
docker run --rm -it --gpus all \
  -v `pwd`/data:/data -v `pwd`/table:/table -v `pwd`/out:/out \
  megadock:hpccm \
    mpirun --allow-run-as-root -n 2 -x OMP_NUM_THREADS=$(nproc) \
      /workspace/megadock-gpu-dp -tb /table/200pairs/200pairs.table
```
