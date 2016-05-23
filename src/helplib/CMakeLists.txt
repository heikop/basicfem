##################################################
###                                            ###
###                  COMPILER                  ###
###                                            ###
##################################################

# some flags for C++ code
list (APPEND CMAKE_CXX_FLAGS "-std=c++11")

##################################################
###                                            ###
###                  LIBRARY                   ###
###                                            ###
##################################################

file (GLOB HELPLIBSOURCES "*.cpp")
add_library (HELPLIB ${HELPLIBSOURCES})
