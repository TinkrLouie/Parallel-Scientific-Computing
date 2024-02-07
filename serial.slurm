#!/bin/bash
#SBATCH --job-name="step-1"
#SBATCH -o step3_slurm.out
#SBATCH -e step3_slurm.err
#SBATCH -p shared
#SBATCH -t 00:01:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=1G       # memory required, up to 250G on standard nodes.
#SBATCH --gres=tmp:1G  # temporary disk space required on the compute node ($TMPDIR), up to 400G

chmod u+x ./step-1.sh
./step-1.sh