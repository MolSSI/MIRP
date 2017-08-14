#############################################
# Check the sha sums on the test data       #
# We check for the 'shasum' and 'sha256sum' #
# commands, and then use that on the test   #
# data                                      #
#############################################
find_program(SHASUM_PROG shasum)
if(NOT SHASUM_PROG)
    find_program(SHASUM_PROG sha256sum)
endif()

if(NOT SHASUM_PROG)
    message(STATUS "shasum or sha256sum not found. I wont check the data files")
else()
    add_test(NAME test_data_checksum
             COMMAND ${SHASUM_PROG} -c sha256sums
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/test_data
    )
endif()



######################################
# Creates a test of the boys function
######################################
macro(add_boys_test filename floattype prec extra_m)
    add_test(NAME boys_${filename}_${floattype}_${prec}_extra_${extra_m}
             COMMAND mirp_run_test --integral boys --file ${CMAKE_SOURCE_DIR}/test_data/${filename} --float ${floattype} --prec ${prec} --extra-m ${extra_m}
    )
endmacro()

##############################################
# Creates a test of single 4-center integrals
##############################################
macro(add_integral4_single_test filename integral floattype prec)
    add_test(NAME ${integral}_${filename}_${floattype}_${prec}
             COMMAND mirp_run_test --integral ${integral} --file ${CMAKE_SOURCE_DIR}/test_data/${filename} --float ${floattype} --prec ${prec}
    )
endmacro()

# Testing of examples
add_test(NAME example_simple_eri COMMAND simple_eri)

# Boys function
add_boys_test(boys_large_range.dat double 0 0)
add_boys_test(boys_large_range.dat double 0 10)
add_boys_test(boys_large_range.dat exact 0 0)
add_boys_test(boys_large_range.dat exact 0 10)
add_boys_test(boys_large_range.dat interval 128  0)
add_boys_test(boys_large_range.dat interval 128 10)
add_boys_test(boys_large_range.dat interval 315  0)
add_boys_test(boys_large_range.dat interval 315 10)
add_boys_test(boys_large_random.dat double 0 0)
add_boys_test(boys_large_random.dat double 0 10)
add_boys_test(boys_large_random.dat exact 0 0)
add_boys_test(boys_large_random.dat exact 0 10)
add_boys_test(boys_large_random.dat interval 128  0)
add_boys_test(boys_large_random.dat interval 128 10)
add_boys_test(boys_large_random.dat interval 315  0)
add_boys_test(boys_large_random.dat interval 315 10)

# ERI
add_integral4_single_test(eri_single_random_1.dat eri_single interval 128)
add_integral4_single_test(eri_single_random_1.dat eri_single interval 315)
add_integral4_single_test(eri_single_random_1.dat eri_single exact 0)
add_integral4_single_test(eri_single_water_sto3g.dat eri_single interval 128)
add_integral4_single_test(eri_single_water_sto3g.dat eri_single interval 315)
add_integral4_single_test(eri_single_water_sto3g.dat eri_single exact 0)
