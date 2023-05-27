#!/bin/bash

source ~/.bashrc
conda activate mashr

strong_file="strong_pairs.MashR_input.txt.gz"
random_file="nominal_pairs.1000000_subset.MashR_input.txt.gz"

Rscript run_MashR.R ${strong_file} ${random_file} 0 ./output_top_paris

Rscript run_MashR.R ${strong_file} ${random_file} 1 ./output_top_paris_across_all