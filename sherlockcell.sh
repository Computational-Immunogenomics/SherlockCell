#!/bin/bash
#SBATCH --job-name=nf_sherlock
#SBATCH --partition=franmartinez
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/adrianparrilla@vhio.org/logs/%x_%j.log
#SBATCH --mem=4G

set -euo pipefail

nextflow run main.nf -profile slurm,singularity -resume