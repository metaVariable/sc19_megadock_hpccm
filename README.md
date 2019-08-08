# Multiple HPC Environments-Aware Container Image Configuration for Bioinformatics Application

- Authors: 
  - Kento Aoyama(1,2) Hiroki Watanabe(1,2) Ohue Masahito(1)  Yutaka Akiyama(1)

- Affiliations:
  1. Department of Computer Science, School of Computing, Tokyo Institute of Technology, Japan
  2. AIST-Tokyo Tech Real World Big-Data Computation Open Innovative Laboratory (RWBC-OIL), National Institute of Advanced Industrial Science and Technology (AIST), Japan

## Repository overview
```sh
.
├── Makefile                        # Makefile used for container building
├── megadock-5.0-alpha-706cb91      # source code of MEGADOCK application
├── megadock_hpccm.py               # HPCCM recipe for generating Dockerfile and Singularity definition
├── README.md                       # this document
├── sample                          # recipe samples used in the poster's experiments
│   ├── Dockerfile                  #   for General Docker environment
│   ├── singularity_ompi-2-1-3.def  #   for TSUBAME3.0 
│   └── singularity_ompi-3-1-3.def  #   for ABCI
── script                          # scripts for experiments
    └── makeTable.sh                #   for generating input docking list (table)

# Those will be generated
.
├── data                            # directory for storing input
└── out                             # directory for storing output and results
```

----

# MEGADOCK-HPCCM

## Requirements

- [HPC Container Maker](https://github.com/NVIDIA/hpc-container-maker/)
- [Docker](https://www.docker.com/) (if you use)
- [Singularity](https://sylabs.io/) (if you use)

## For Singularity environment

**Note: Singularity environment may depend on host libralies**

- requirements
  - pip, python (for HPCCM)
  - singularity
    - require privilege for `sudo singularity build`

### Setup on your local environment

```sh
# install hpccm
sudo pip install hpccm

# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm
```

### Generate Singularity definition, build singularity image
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

### Run Singularity container on HPC environment

#### Setup on HPC environment

```sh
# clone MEGADOCK-HPCCM repository
git clone https://github.com/metaVariable/sc19_megadock_hpccm.git
cd sc19_megadock_hpccm

# download benchmark dataset (ZDOCK Benchmark 1.0)
mkdir -p data
wget http://zlab.umassmed.edu/zdock/benchmark1.0.tar.gz
tar xvzf benchmark1.0.tar.gz -C data && rm benchmark1.0.tar.gz
```

**Note: following commands should be running on compute-node.** 

- Please read the system manual about singularity on your HPC system. We must add several options for singularity runtime in general case.
  - e.g.) Volume option (`-B`) for mounting system storage, applications, libraries, etc.

```sh
# create input table by script
INTERACTIVE=0 ENABLE_TSV=1 TSV_SIZE=50 \
script/makeTable.sh data/benchmark1.0/unbound_pdb/\*_r.pdb data/benchmark1.0/unbound_pdb/\*_l.pdb test_all2all

# singularity exec (please replace ${SINGULARITY_IMAGE} to your path to singularity image)
singularity exec --nv ${SINGULARITY_IMAGE} \
  mpirun -n 2 -x OMP_NUM_THREADS=$(nproc) \
  /workspace/megadock-gpu-dp -tb `pwd`/out/test_all2all/test_all2all.table

# singularity exec with host MPI library (please read carefully the system manual on your HPC system)
mpirun -n 16 -x OMP_NUM_THREADS=$(nproc) \
  singularity exec --nv ${SINGULARITY_IMAGE} \
  /workspace/megadock-gpu-dp -tb `pwd`/out/test_all2all/test_all2all.table
```

----

## For Docker environment

- requirements
  - pip, python (for HPCCM)
  - docker > 19.03

### Setup

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

### Run Docker container on HPC environment (WIP)

**Note: This section is under development**

```sh
# download benchmark dataset (ZDOCK Benchmark 1.0)
mkdir -p data
wget http://zlab.umassmed.edu/zdock/benchmark1.0.tar.gz
tar xvzf benchmark1.0.tar.gz -C data && rm benchmark1.0.tar.gz

# create input table by script
INTERACTIVE=0 ENABLE_TSV=1 TSV_SIZE=50 \
script/makeTable.sh data/benchmark1.0/unbound_pdb/\*_r.pdb data/benchmark1.0/unbound_pdb/\*_l.pdb test_all2all

# run 
docker run --rm -it --gpus all megadock:hpccm
  mpirun --allow-run-as-root -n 2 megadock-gpu-dp -tb test_all2all.table
```