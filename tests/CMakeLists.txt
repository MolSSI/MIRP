
macro(add_boys_test filename floattype prec extra_m)
    add_test(NAME boys_${filename}_${floattype}_${prec}_extra_${extra_m}
             COMMAND mirp_run_test --integral boys --filename ${CMAKE_SOURCE_DIR}/test_data/${filename} --float ${floattype} --prec ${prec} --extra-m ${extra_m}
    )
endmacro()

add_test(NAME example_simple_eri COMMAND simple_eri)
add_boys_test(boys_large_range.dat interval 128  0)
add_boys_test(boys_large_range.dat interval 128 10)
add_boys_test(boys_large_range.dat interval 315  0)
add_boys_test(boys_large_range.dat interval 315 10)
add_boys_test(boys_large_random.dat interval 128  0)
add_boys_test(boys_large_random.dat interval 128 10)
add_boys_test(boys_large_random.dat interval 315  0)
add_boys_test(boys_large_random.dat interval 315 10)
