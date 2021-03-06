##################################################
###                                            ###
###                  COMPILER                  ###
###                                            ###
##################################################

# some flags for C++ code
list (APPEND CMAKE_CXX_FLAGS "-std=c++11")


##################################################
###                                            ###
###                    ABCD                    ###
###                                            ###
##################################################

include(CTest) # this includes enable_testing() as well

file (GLOB TESTS "*_test.cpp")
FOREACH(testpath ${TESTS})
    string(REGEX REPLACE "([^ ]*)/test/" "" filename ${testpath})
    string(REPLACE ".cpp" "" testname ${filename})
    add_executable(${testname} ${filename})
    TARGET_LINK_LIBRARIES(${testname} TYPELIB HELPLIB GLOBLIB ${MPI_LIBRARIES})
    add_test(${testname} ${testname}) # for ctest
    add_test(${testname}mpi mpirun -np 4 ${testname}) # for ctest
    #add_test(${testname}mpi mpirun ${testname}) # for ctest
ENDFOREACH(testpath)
