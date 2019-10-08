```sh
# Sample Recipes

# For general docker environment
hpccm --recipe megadock_hpccm.py --format docker > Dockerfile


# For TSUBAME3.0
# target modules: cuda/8.0.61 openmpi/2.1.2-opa10.9

hpccm --recipe megadock_hpccm.py --format singularity --userarg ompi=2.1.2 fftw=3.3.8 opa=True > singularity_ompi-2-1-2_opa.def


# For ABCI
# target modules: cuda/10.0/10.0.130 openmpi/2.1.6

hpccm --recipe megadock_hpccm.py --format singularity --userarg ompi=2.1.6 fftw=3.3.8 ofed=True > singularity_ompi-2-1-6_ofed.def
```