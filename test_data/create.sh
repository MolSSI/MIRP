#!/bin/bash

#######################################
# This is just a script to create the
# test files in this directory. Mostly
# used by Ben to test regeneration
# of the tests after changes
#######################################

set -eu

##################
# Boys function
##################
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


#############################
# Generic 4-center integrals
#############################
generator/generate_integral_single_random.py \
           --filename 4center_single_random_1.inp \
           --seed 64020964 --ncenter 4 \
           --alpha-power 10 --xyz-power 1 --max-am 3 \
           --ntest 2000 --ndigits 20

generator/generate_integral_single_frombasis.py \
           --filename 4center_single_water_sto3g.inp \
           --geo generator/geometry/water.xyz \
           --basis generator/basis/sto-3g.gbs \
           --ncenter 4 --ndigits 20 --ntests 2000 --seed 770225234

generator/generate_integral_random.py \
           --filename 4center_random_1.inp \
           --seed 23833210 --ncenter 4 \
           --max-am 1 --max-nprim 2 --max-ngen 2 \
           --xyz-power 1 --alpha-power 10 --coeff-power 5 \
           --max-z 30 --ntest 1000 --ndigits 20

generator/generate_integral_frombasis.py \
           --filename 4center_water_sto3g.inp \
           --basis generator/basis/sto-3g.gbs \
           --geometry generator/geometry/water.xyz \
           --ndigits 20 --ncenter 4 --ntests 1000 --seed 661232124


##########################
# ERI References
##########################
../build/mirp_bin/mirp_create_test \
           --infile 4center_single_random_1.inp \
           --outfile eri_single_random_1.dat \
           --integral eri_single --ndigits 101

../build/mirp_bin/mirp_create_test \
           --infile 4center_single_water_sto3g.inp \
           --outfile eri_single_water_sto3g.dat \
           --integral eri_single --ndigits 101

../build/mirp_bin/mirp_create_test \
           --infile 4center_random_1.inp \
           --outfile eri_random_1.dat \
           --integral eri --ndigits 101

../build/mirp_bin/mirp_create_test \
           --infile 4center_water_sto3g.inp \
           --outfile eri_water_sto3g.dat \
           --integral eri --ndigits 101
