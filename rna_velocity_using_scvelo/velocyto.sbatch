#!/bin/bash
#SBATCH --job-name=velocyto
#SBATCH --output=velocyto-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=80GB
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=user@ufl.edu
#SBATCH --qos=kirst

module load velocyto/0.17

################################
## Parameters set by the user ##

# Note that other options are available for velocyto run10x, see https://velocyto.org/velocyto.py/tutorial/cli.html. For example, if you have a GTF with elements that may be captured by 10x scRNA-seq but that should not be considered, such as transposons or rRNA, you can use the -m.

GTF=/blue/... #path to your GTF
CELLRANGER_OUT=/blue/... #path to your cellranger output folder
N_CORES=6 #number of cores to use
################################

velocyto run10x -vvv --samtools-threads N_CORES CELLRANGER_OUT GTF

