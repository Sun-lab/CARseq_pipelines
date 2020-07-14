#!/bin/bash
mkdir -p sbatch_log
for i in `seq 1 360`; do
        eval 'sbatch -p bat,int -c 12 --mem=20G -t 8:00:00 --output=sbatch_log/hostname_%j.out --error=sbatch_log/hostname_%j.err --wrap="Rscript simulation_step3_test_methods.R '${i}'"'
        sleep 0.1
done
