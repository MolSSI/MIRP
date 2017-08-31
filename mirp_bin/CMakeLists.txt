add_library(test_common OBJECT cmdline.cpp
                               testfile_io.cpp
                               reffile_io.cpp
                               read_construct_basis.cpp
                               test_boys.cpp
                               test_common.cpp
                               test_integral4.cpp
                               test_integral4_single.cpp
                               ref_integral4.cpp
)

# Add the include directories to the object library
target_include_directories(test_common PRIVATE
                           $<TARGET_PROPERTY:mirp,INTERFACE_INCLUDE_DIRECTORIES>)

# The executables themselves
add_executable(mirp_run_test         mirp_run_test.cpp         $<TARGET_OBJECTS:test_common>) 
add_executable(mirp_create_test      mirp_create_test.cpp      $<TARGET_OBJECTS:test_common>)
add_executable(mirp_create_reference mirp_create_reference.cpp $<TARGET_OBJECTS:test_common>)
add_executable(mirp_test_reference   mirp_test_reference.cpp   $<TARGET_OBJECTS:test_common>)

# Link these to mirp. The dependency and include directories
# will be included through here as well (they were added as PUBLIC)
target_link_libraries(mirp_run_test         PRIVATE mirp)
target_link_libraries(mirp_create_test      PRIVATE mirp)
target_link_libraries(mirp_create_reference PRIVATE mirp)
target_link_libraries(mirp_test_reference   PRIVATE mirp)
