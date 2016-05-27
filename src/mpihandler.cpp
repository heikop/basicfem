#include "mpihandler.hpp"

MpiHandler::MpiHandler():
    _size{0},
    _rank{0}
{
    MPI::Init();
    MPICALL(_size = MPI::COMM_WORLD.Get_size();)
    MPICALL(_rank = MPI::COMM_WORLD.Get_rank();)
}//MpiHandler::MpiHandler()

MpiHandler::~MpiHandler()
{
    MPI::Finalize();
}//MpiHandler::~MpiHandler()

size_t MpiHandler::get_global_size()
{
    return _size;
}//size_t MpiHandler::get_global_size()

size_t MpiHandler::get_global_rank()
{
    return _rank;
}//size_t MpiHandler::get_global_rank()

void MpiHandler::catch_call(int line, std::string file, const std::exception& ex)
{
    //TODO
    std::cout << ex.what() << std::endl;
}//void MpiHandler::catch_call(int line, std::string file, const std::exception& ex)
