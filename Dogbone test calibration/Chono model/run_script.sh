#!/bin/bash
#SBATCH --account=p31861  ## YOUR ACCOUNT pXXXX or bXXXX
#SBATCH --partition=normal  ### PARTITION (buyin, short, normal, etc)
#SBATCH --nodes=1 ## how many computers do you need
#SBATCH --ntasks-per-node=1 ## how many cpus or processors do you need on each computer
#SBATCH --time=8:00:00 ## how long does this need to run (remember different partitions have restrictions on this param)
#SBATCH --mem-per-cpu=3G ## how mdogbone2h RAM do you need per CPU (this effects your FairShare score so be careful to not ask for more than you need))
#SBATCH --job-name=CSL_EXP  ## When you run squeue -u NETID this is how you can identify the job
#SBATCH --constraint=quest10

echo "=== SLURM INFO ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $HOSTNAME"
echo "Start time: $(date)"

echo "=== CPU INFO ==="
lscpu
echo "================="

# load modules
module purge
module load singularity

cd /projects/p31861/Users/YuKe/workdir/dogbone2

rm -r build

mkdir -p /projects/p31861/Users/YuKe/workdir/dogbone2/build

pwd

# configure chrono-concrete
singularity exec --pwd /dogbone2/build -B /projects/p31861/Users/YuKe/workdir/dogbone2:/dogbone2 -B /projects/p31861/Users/LaleErol/chrono-concrete:/chrono-concrete /projects/p31861/Users/YuKe/project-chrono-dependencies-with-intel-mkl.sif cmake ../ -G "Ninja" \
        -DCMAKE_BUILD_TYPE=Release \
        -DENABLE_MODULE_POSTPROCESS=TRUE \
        -DENABLE_MODULE_PYTHON=TRUE \
        -DENABLE_MODULE_COSIMULATION=TRUE \
        -DENABLE_MODULE_IRRLICHT=TRUE \
        -DENABLE_MODULE_VEHICLE=FALSE \
        -DENABLE_MODULE_MULTICORE=TRUE \
        -DENABLE_MODULE_OPENGL=TRUE \
        -DENABLE_MODULE_FSI=TRUE \
        -DENABLE_MODULE_SYNCHRONO=TRUE \
        -DENABLE_MODULE_CSHARP=FALSE \
        -DENABLE_MODULE_GPU=TRUE \
        -DENABLE_MODULE_DISTRIBUTED=TRUE \
        -DENABLE_HDF5=TRUE \
        -DCMAKE_C_COMPILER=/usr/bin/gcc \
        -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
        -DCUDA_HOST_COMPILER=/usr/bin/gcc \
        -DPYTHON_EXECUTABLE=/usr/bin/python3 \
        -DEIGEN3_INCLUDE_DIR=/usr/include/eigen3 \
        -DCMAKE_VERBOSE_MAKEFILE=TRUE \
        -DChrono_DIR=/chrono-concrete/build/cmake \
        -DChrono_DATA_DIR=/chrono-concrete/build/data \
        -DCHRONO_CXX_FLAGS=TRUE\
        -DCHORONO_C_FLAGS=TRUE
pwd
# build project chrono-concrete
singularity exec --pwd /dogbone2/build -B /projects/p31861/Users/YuKe/workdir/dogbone2:/dogbone2 -B /projects/p31861/Users/LaleErol/chrono-concrete:/chrono-concrete /projects/p31861/Users/YuKe/project-chrono-dependencies-with-intel-mkl.sif ninja -j 4

# Run job
cd /projects/p31861/Users/YuKe/workdir/dogbone2

singularity exec -B /projects/p31861/Users/YuKe/workdir/dogbone2:/dogbone2 -B /projects/p31861/Users/LaleErol/chrono-concrete:/chrono-concrete /projects/p31861/Users/YuKe/project-chrono-dependencies-with-intel-mkl.sif /dogbone2/build/fea_demo
