
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

target_include_directories(mirp PUBLIC ${MIRP_DEPS_INCLUDE_DIRS})

target_link_libraries(mirp PUBLIC ${MIRP_DEPS_LIBRARIES})


# Installing the library
install(TARGETS mirp
        EXPORT mirpTargets
        ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
        LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
        
# Install the header files
install(DIRECTORY "${CMAKE_CURRENT_LIST_DIR}"
        DESTINATION include
        FILES_MATCHING REGEX ".*\.h$"
)
