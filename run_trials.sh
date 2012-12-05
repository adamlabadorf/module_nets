#!/bin/bash

function wait_on_pid {
    echo "waiting for $PID to finish"
    PID=$1
    ps $PID
    STATUS=$?
    while [ $STATUS == 0 ]
    do
        sleep 10
        ps $PID > /dev/null
        STATUS=$?
    done
}

RUNID=$RANDOM

#for i in {1..8}; do
#    ./run_mn_binding.m trials/${RUNID}_trial_easy_${i}.dat easy > /dev/null &
#    PID=$!
#done
#echo $PID
#wait_on_pid $PID

for i in {1..8}; do
    ./run_mn_binding.m trials/${RUNID}_trial_medium_${i}.dat medium > /dev/null &
    PID=$!
done
wait_on_pid $PID

#for i in {1..8}; do
#    ./run_mn_binding.m trials/${RUNID}_trial_hard_${i}.dat hard > /dev/null &
#done
