##################################################
###                                            ###
###                  GENERAL                   ###
###                                            ###
##################################################

cmake_minimum_required (VERSION 2.8)
project (myfem)
# enable testing also in the build dir (not only in test subdir)
include(CTest) # this includes enable_testing() as well


##################################################
###                                            ###
###                  PACKAGES                  ###
###                                            ###
##################################################

# TODO add doxygen(?) and OpenCL or CUDA(?)
find_package (MPI REQUIRED)

include_directories (/opt/openmpi/include)


##################################################
###                                            ###
###               SUBDIRECTORIES               ###
###                                            ###
##################################################

add_subdirectory (src)
add_subdirectory (test)
#add_subdirectory (demos)


##################################################
###                                            ###
###               DOCUMENTATION                ###
###                                            ###
##################################################

# TODO (later)
