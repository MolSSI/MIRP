
list(APPEND MIRP_FILELIST
               math.c
               gpt.c
               shell.c

               kernels/integral4_wrappers.c

               kernels/boys_d.c
               kernels/boys.c

               kernels/eri_d.c
               kernels/eri.c
)

add_library(mirp SHARED ${MIRP_FILELIST})

target_link_libraries(mirp "${MIRP_DEPS_TARGETS}")

# Set the include directory for installing
set_target_properties(mirp PROPERTIES
                      INTERFACE_INCLUDE_DIRECTORIES "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# Installing the library
install(TARGETS mirp
        EXPORT mirpTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
        
# Install the header files
install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING REGEX ".*\.h$"
)
