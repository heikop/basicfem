##################################################
###                                            ###
###               SUBDIRECTORIES               ###
###                                            ###
##################################################

add_subdirectory (helplib)
add_subdirectory (types)
#add_subdirectory (solver)
#add_subdirectory (assembler)
#add_subdirectory (dataprocessing)

list (APPEND CMAKE_CXX_FLAGS "-std=c++11")
file (GLOB GLOBSOURCES "*.cpp")
#file (GLOB GLOBSOURCES "mpihandler.cpp")
#file (GLOB GLOBSOURCES "global.cpp")
add_library (GLOBLIB ${GLOBSOURCES})
