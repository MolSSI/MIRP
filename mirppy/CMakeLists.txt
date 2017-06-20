find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

add_library(mirppy MODULE 
                   export_boys.c
                   export.c)

set_target_properties(mirppy PROPERTIES PREFIX "")
target_link_libraries(mirppy mirp ${PYTHON_LIBRARIES} arb mpfr gmp m)
