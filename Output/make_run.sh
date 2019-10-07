#!/bin/bash          
#/* Copyright (c) 2016 Siddhartha Shelton */

rebuild=true
run=true
run_compare=false
compare_only=false
with_manual_input=false
manual_input_only=false
get_flag_list=false

make_checkpoint=false   
run_from_checkpoint=false


if $rebuild
then

    cd nbody_test
    cmake  -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DSEPARATION=OFF -DNBODY_OPENMP=ON    ~/milkywayathome_client/
    make -j 
fi

cd ~/nbody_test/bin
cp ~/shift.txt ./

if $run 
then
    ./milkyway_nbody \
    -f ~/for_developers.lua \
    -z path_to_input_hist \
    -o some_output.out \
    -n 8 -b -P\
    -i 3.95 1.0 0.2 0.2 12 0.2 \
    
fi

if $run_compare
then
    ./milkyway_nbody \
    -f ~/Desktop/research/lua/full_control.lua \
    -o some_output.out \
    -h path_to_input_hist \
    -z path_to_output_hist \
    -n 8 -b  -P -u\
    -i 3.95 1.0 0.2 0.2 12 0.2 \

fi
######################################
if $make_checkpoint
then
    ./milkyway_nbody \
    -f ~/Desktop/research/lua/full_control.lua \
    -z testing_checkpoint.hist \
    -o testing_checkpoint.out \
    -n 8 -b -P  --checkpoint-interval=3 \
    -i --no-clean-checkpoint 3.95 1.0 0.2 0.2 12 0.2 \
    
fi

if $run_from_checkpoint
then
    ./milkyway_nbody \
    -f ~/Desktop/research/lua/full_control.lua \
    -z testing_checkpoint.hist \
    -o testing_checkpoint.out \
    -n 8 -b  --checkpoint=nbody_checkpoint\
     3.95 1.0 0.2 0.2 12 0.2\
     
fi
######################################

#OPTIONS:
#-s -> compare using only emd and cost component
#-S -> use emd, cost, beta dispersion
#-V -> use emd, cost, velocity dispersion
#-D -> use emd, cost, beta dispersion and velocity dispersion
if $compare_only 
then
    ./milkyway_nbody \
    -h testing_checkpoint.hist \
    -S testing_checkpoint1.hist \

fi


# if you run:
if $get_flag_list
then
    ./milkyway_nbody --help
# it wil show you what all the flags mean
fi

if $with_manual_input
then
    ./milkyway_nbody \
    -f ~/Desktop/research/milkywayathome_client/nbody/sample_workunits/for_dev_manual_body.lua \
    -o some_output.out \
    -z path_to_output_hist \
    -n 8 -b  -u --visualizer-bin=~/Desktop/research/nbody_test/bin/milkyway_nbody_graphics\
    -i 4.0 1 0.2 0.5 12 0.5 ~/Desktop/research/disk.out\

fi

if $manual_input_only
then
    ./milkyway_nbody \
    -f ~/Desktop/research/milkywayathome_client/nbody/sample_workunits/manual_body_input.lua \
    -o some_output.out \
    -z path_to_output_hist \
    -n 8 -b  -u\
    -i 4.0 ~/Desktop/research/disk.out\

fi
