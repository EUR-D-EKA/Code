#!/bin/bash          
#/* Copyright (c) 2016 Siddhartha Shelton */
    #rm -r nbody_test
    #mkdir nbody_test
    #cd nbody_test
    #cmake  -DCMAKE_BUILD_TYPE=Release -DNBODY_DEV_OPTIONS=ON -DNBODY_GL=OFF -DBOINC_APPLICATION=OFF -DSEPARATION=OFF -DNBODY_OPENMP=ON    ~/milkywayathome_client/
    #make -j
 



    ../nbody_test/bin/milkyway_nbody \
    -f ./for_developers.lua \
    -o some_output.out \
    -z path_to_output_hist \
    -n 8 -b -P \
    -i 3.95 1.0 0.2 0.2 12 0.2\

# -h path_to_input_hist \
#if you run:
#     ./milkyway_nbody --help\
# it wil show you what all the flags mean
