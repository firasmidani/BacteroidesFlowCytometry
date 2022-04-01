#!/bin/bash

# PURPOSE
# this is a driver script for processing flow cytometry data from the Bacteorides experiment in July of 2018

# source Python virtual environment specifci to this project
source /data/davidlab/users/fsm/envs/miniconda3/bin/activate py2env

# for each FCS sample, infer clusters, cluster abundances, and total count of events
python /home/lad44/davidlab/users/fsm/bacteroides_20220307/code/analysis/process_fcs_samples/InferComponents.py $1
