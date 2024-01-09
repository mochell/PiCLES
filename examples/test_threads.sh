#!/bin/bash
#PBS -A p93300612
#PBS -N test_threads
#PBS -q main 
#PBS -l walltime=01:00:00 
#PBS -l select=1:ncpus=128:mpiprocs=1:ompthreads=128
#PBS -j oe

module purge
#module load ncarenv/23.09 julia
module load peak-memusage

which julia

#for threads in 2 ; do
for threads in 1 2 4 8 16 32 64 128 ; do
      export JULIA_NUM_THREADS=${threads}
        echo "Num threads: $JULIA_NUM_THREADS"
          peak_memusage julia ./exmpl_homogenous_box_mthread.jl

            echo ""
            done
