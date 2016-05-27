#ifndef __MPIHANDLER_HPP_
#define __MPIHANDLER_HPP_

#include <mpi.h>
#include <exception>
#include <string>
#include <iostream> //TODO (later) remove this

class MpiHandler
{
public:
    MpiHandler();
    MpiHandler(const MpiHandler&) = delete;
    MpiHandler(MpiHandler&&) = delete;
    ~MpiHandler();
    MpiHandler& operator=(const MpiHandler&) = delete;
    MpiHandler& operator=(MpiHandler&&) = delete;

    size_t get_global_size();
    size_t get_global_rank();
private:
    size_t _size;
    size_t _rank;
    void catch_call(int line, std::string file, const std::exception& ex);
};//class MpiHander

#define MPICALL(X) try { X } \
                   catch (const std::exception& ex) \
                   { catch_call(__LINE__, __FILE__, ex); } 

#endif//__MPIHANDLER_HPP_
