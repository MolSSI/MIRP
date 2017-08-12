#!/bin/bash

generator/generate_boys_random.py \
           --filename boys_large_random.inp \
           --max-m 100 --power 32 --seed 828843219 \
           --ndigits 32 --ntest 25000

../build/mirp_bin/mirp_create_test \
           --infile boys_large_random.inp \
           --outfile boys_large_random.dat \
           --integral boys --ndigits 101

generator/generate_boys_range.py \
           --filename boys_large_range.inp \
           --max-m 100 --power 32

../build/mirp_bin/mirp_create_test \
           --infile boys_large_range.inp \
           --outfile boys_large_range.dat \
           --integral boys --ndigits 101
