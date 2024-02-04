#!/bin/bash
#SBATCH --job-name="step-1"
#SBATCH -o step1.%A.out
#SBATCH -e step1.%A.err
#SBATCH -p shared
#SBATCH -t 00:01:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G       # memory required, up to 250G on standard nodes.
#SBATCH --gres=tmp:1G  # temporary disk space required on the compute node ($TMPDIR), up to 400G

module load gcc

g++ -fopenmp -O3 -march=native -std=c++0x -fno-math-errno -c -o step-1.o step-1.cpp

g++ -fopenmp -O3 -march=native -std=c++0x -fno-math-errno step-1 NBodySimulation.o step-1.o


./step-1 0.01  10.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0