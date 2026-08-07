#!/bin/bash
export MSCD_GPU=1
export OMP_NUM_THREADS=1 # To make it purely MPI scaling if they want, wait, by default OpenMP might oversubscribe if we use MPI. The user's run script probably didn't set OMP_NUM_THREADS, but they used --bind-to none.
# Actually, the user's instructions for V5 were:
# "O Open MPI amarra o processo a um núcleo com np baixo, e as 12 threads OpenMP ficaram empilhadas nesse núcleo. Com --bind-to none: pathcut 1,706s"

# If I use `mpirun --use-hwthread-cpus --bind-to none -np $PE`, it means $PE MPI processes, each spawning max threads?
# Let's just run it the way the user did.

rm -f benchmark_gpu_scaling.log
echo "=== GPU V0 vs GPU V2 Scaling Benchmark on 316 atoms ===" | tee -a benchmark_gpu_scaling.log
echo "Date: $(date)" | tee -a benchmark_gpu_scaling.log

for pe in 1 2 4 6 8 12; do
    echo "----------------------------------------" | tee -a benchmark_gpu_scaling.log
    echo "Running PE = $pe (GPU V0)" | tee -a benchmark_gpu_scaling.log
    
    /usr/bin/time -p mpirun --use-hwthread-cpus --bind-to none -np $pe ./baseline/randmscd_gpu.v0 1x2iron.in >> benchmark_gpu_scaling.log 2>&1
    
    echo "Pausing for 40 seconds to cool down..." | tee -a benchmark_gpu_scaling.log
    sleep 40
    
    echo "Running PE = $pe (GPU V2)" | tee -a benchmark_gpu_scaling.log
    /usr/bin/time -p mpirun --use-hwthread-cpus --bind-to none -np $pe ./randmscd_gpu_v2 1x2iron.in >> benchmark_gpu_scaling.log 2>&1
    
    echo "Pausing for 40 seconds to cool down..." | tee -a benchmark_gpu_scaling.log
    sleep 40
done
echo "Benchmark completed." | tee -a benchmark_gpu_scaling.log
