include(CMakeMacros.txt)

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
    message(STATUS "shasum or sha256sum not found. I will not check the data files")
else()
    add_test(NAME tests_checksum
             COMMAND ${SHASUM_PROG} -c sha256sums
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests
    )
endif()


######################################
# Tests the help screen of test progs
######################################
add_test(NAME help_mirp_create_test_1 COMMAND mirp_create_test)
add_test(NAME help_mirp_create_test_2 COMMAND mirp_create_test -h)
add_test(NAME help_mirp_verify_test_1 COMMAND mirp_verify_test)
add_test(NAME help_mirp_verify_test_2 COMMAND mirp_verify_test -h)
add_test(NAME help_mirp_create_reference_1 COMMAND mirp_create_reference)
add_test(NAME help_mirp_create_reference_2 COMMAND mirp_create_reference -h)
add_test(NAME help_mirp_verify_reference_1 COMMAND mirp_verify_reference)
add_test(NAME help_mirp_verify_reference_2 COMMAND mirp_verify_reference -h)

#############################################
# Test failures
# This makes sure the testing infrastructure
# is sane
#############################################
__verify_test_boys(${CMAKE_CURRENT_LIST_DIR}/boys_failure_1.dat interval 332 0 "1 / 1 failed")
__verify_test_boys(${CMAKE_CURRENT_LIST_DIR}/boys_failure_1.dat interval 332 10 "1 / 1 failed")
__verify_test(${CMAKE_CURRENT_LIST_DIR}/gtoeri_failure_1.dat gtoeri interval 332 "1 / 1 failed")


################
# Boys function
################
verify_test_boys(${CMAKE_CURRENT_LIST_DIR}/boys_large_range.dat)
verify_test_boys(${CMAKE_CURRENT_LIST_DIR}/boys_large_random.dat)
create_and_verify_test_boys(${CMAKE_CURRENT_LIST_DIR}/boys_large_random.inp)


############
# ERI
############
verify_test(${CMAKE_CURRENT_LIST_DIR}/gtoeri_single_random_1.dat gtoeri_single)
verify_test(${CMAKE_CURRENT_LIST_DIR}/gtoeri_single_water_sto-3g.dat gtoeri_single)

verify_test(${CMAKE_CURRENT_LIST_DIR}/gtoeri_random_1.dat gtoeri)
verify_test(${CMAKE_CURRENT_LIST_DIR}/gtoeri_water_sto-3g.dat gtoeri)

verify_reference(${CMAKE_CURRENT_LIST_DIR}/gtoeri_water_sto-3g.ref gtoeri)


create_and_verify_test(${CMAKE_CURRENT_LIST_DIR}/4center_single_water_sto-3g.inp gtoeri_single)
create_and_verify_test(${CMAKE_CURRENT_LIST_DIR}/4center_water_sto-3g.inp gtoeri)
create_and_verify_reference(gtoeri)
