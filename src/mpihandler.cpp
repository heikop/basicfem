#include "mpihandler.hpp"

MpiHandler::MpiHandler():
    _size{0},
    _rank{0}
{
    std::cout << "start mpi" << std::endl;
}

MpiHandler::~MpiHandler()
{
    std::cout << "finish mpi" << std::endl;
}

size_t MpiHandler::get_size()
{
    return _size;
}

size_t MpiHandler::get_rank()
{
    return _rank;
}

void MpiHandler::catch_call(int line, std::string& file, const std::exception& ex)
{
}
