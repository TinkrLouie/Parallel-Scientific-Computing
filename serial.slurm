#!/bin/bash
#SBATCH --job-name="step-1"
#SBATCH -o results/step1.%A.out
#SBATCH -e results/step1.%A.err
#SBATCH -p shared
#SBATCH -t 00:01:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=1G       # memory required, up to 250G on standard nodes.
#SBATCH --gres=tmp:1G  # temporary disk space required on the compute node ($TMPDIR), up to 400G

./step-3.sh