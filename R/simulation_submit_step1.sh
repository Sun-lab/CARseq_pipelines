#!/bin/bash
mkdir -p sbatch_log
mkdir -p ../data/get_distribution_of_parameters_from_real_data_results
for batch in `seq 1 567`; do
    if [ ! -f ../data/get_distribution_of_parameters_from_real_data_results/estimates_mat${batch}.rds ]; then
        eval 'sbatch -p bat -c 1 --mem=2G -t 200:00:00 --output=sbatch_log/hostname_%j.out --error=sbatch_log/hostname_%j.err --wrap="Rscript step1_get_distribution_of_parameters_from_real_data.R '${batch}'"'
        sleep 0.1
    fi
done
