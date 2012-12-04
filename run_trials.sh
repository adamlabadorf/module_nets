#!/bin/bash

RUNID=$RANDOM

for i in {1..20}; do
    ./run_mn_binding.m trials/${RUNID}_trial_easy_${i}.dat easy
done

for i in {1..5}; do
    ./run_mn_binding.m trials/${RUNID}_trial_medium_${i}.dat medium
done
