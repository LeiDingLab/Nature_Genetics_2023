#!/bin/bash
#SBATCH -N 1
#SBATCH -p RM-shared
#SBATCH -t 47:59:59
#SBATCH --ntasks-per-node=64

date

~/.conda/envs/CellBender/bin/cellbender remove-background \
	--input ./data/LH002/raw_feature_bc_matrix.h5 \
	--output ./result/20230310_LH002_CellBender_output.h5 \
	--expected-cells 7000 \
	--total-droplets-included 15000 \
	--fpr 0.01 \
	--epochs 150

echo "Completed"
date
