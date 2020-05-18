#!/bin/bash

#SBATCH -t 06:00:00
#SBATCH --mem=50G
#SBATCH -o /hpc/pmc_gen/lkester/CNV_calling/stdout/SlurmCNVCallingGATK%A_%a.out
#SBATCH -e /hpc/pmc_gen/lkester/CNV_calling/stdout/SlurmCNVCallingGATK%A_%a.err
