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
    add_test(NAME tests_checksum
             COMMAND ${SHASUM_PROG} -c sha256sums
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/tests
    )
endif()



######################################
# Creates a test of the boys function
######################################
macro(add_boys_test filename floattype prec extra_m)
    if(${prec} EQUAL 0)
        add_test(NAME boys_${filename}_${floattype}_extra_${extra_m}
                 COMMAND mirp_verify_test --integral boys
                                          --file ${CMAKE_CURRENT_LIST_DIR}/${filename}
                                          --float ${floattype}
                                          --extra-m ${extra_m}
        )
    else()
        add_test(NAME boys_${filename}_${floattype}_${prec}_extra_${extra_m}
                 COMMAND mirp_verify_test --integral boys
                                          --file ${CMAKE_CURRENT_LIST_DIR}/${filename}
                                          --float ${floattype} --prec ${prec}
                                          --extra-m ${extra_m}
        )
    endif()
endmacro()

##############################################
# Creates a test of general integrals
##############################################
macro(add_integral_test filename integral floattype prec)
    if(${prec} EQUAL 0)
        add_test(NAME ${integral}_${filename}_${floattype}
                 COMMAND mirp_verify_test --integral ${integral} 
                                           --file ${CMAKE_CURRENT_LIST_DIR}/${filename}
                                           --float ${floattype}
        )
    else()
        add_test(NAME ${integral}_${filename}_${floattype}_${prec}
                 COMMAND mirp_verify_test --integral ${integral} 
                                       --file ${CMAKE_CURRENT_LIST_DIR}/${filename}
                                       --float ${floattype} --prec ${prec}
        )
    endif()
endmacro()


##############################################
# Creates a test of reference files
##############################################
macro(add_integral_reftest filename integral)
    add_test(NAME reftest_${integral}_${filename}
             COMMAND mirp_verify_reference --integral ${integral} 
                                           --file ${CMAKE_CURRENT_LIST_DIR}/${filename}
    )
endmacro()


# Boys function
add_boys_test(boys_large_range.dat exact 0 0)
add_boys_test(boys_large_range.dat exact 0 10)
add_boys_test(boys_large_range.dat interval 128  0)
add_boys_test(boys_large_range.dat interval 128 10)
add_boys_test(boys_large_range.dat interval 332  0)
add_boys_test(boys_large_range.dat interval 332 10)
add_boys_test(boys_large_random.dat exact 0 0)
add_boys_test(boys_large_random.dat exact 0 10)
add_boys_test(boys_large_random.dat interval 128  0)
add_boys_test(boys_large_random.dat interval 128 10)
add_boys_test(boys_large_random.dat interval 332  0)
add_boys_test(boys_large_random.dat interval 332 10)

# ERI
add_integral_test(eri_single_random_1.dat gtoeri_single interval 128)
add_integral_test(eri_single_random_1.dat gtoeri_single interval 332)
add_integral_test(eri_single_random_1.dat gtoeri_single exact 0)
add_integral_test(eri_single_water_sto-3g.dat gtoeri_single interval 128)
add_integral_test(eri_single_water_sto-3g.dat gtoeri_single interval 332)
add_integral_test(eri_single_water_sto-3g.dat gtoeri_single exact 0)

add_integral_test(eri_random_1.dat gtoeri interval 128)
add_integral_test(eri_random_1.dat gtoeri interval 332)
add_integral_test(eri_random_1.dat gtoeri exact 0)
add_integral_test(eri_water_sto-3g.dat gtoeri interval 128)
add_integral_test(eri_water_sto-3g.dat gtoeri interval 332)
add_integral_test(eri_water_sto-3g.dat gtoeri exact 0)

add_integral_reftest(eri_water_sto-3g.ref gtoeri)
