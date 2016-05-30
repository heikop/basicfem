#include "include/vector.hpp"

namespace hptypes
{

Vector::Vector(const Vector& other):
    _size_global{other._size_global}, _size_local{other._size_local},
    _firstentrynumber{other._firstentrynumber},
    _data{new double[_size_local]}
{
    for (size_t i{0}; i < _size_local; ++i)
        _data[i] = other._data[i];
    //TODO error handling allocating
}//Vector::Vector(const Vector& other)

Vector::Vector(Vector&& other):
    _size_global{other._size_global}, _size_local{other._size_local},
    _firstentrynumber{other._firstentrynumber},
    _data{other._data}
{
    other._size_global = 0;
    other._size_local = 0;
    other._data = nullptr;
}//Vector::Vector(Vector&& other)

Vector::Vector(const size_t size):
    _size_global{size}, _size_local{size},
    _data{new double[_size_local]}
{
    for (size_t i{0}; i < _size_local; ++i)
        _data[i] = 0.0;
    //TODO error handling allocating
    _size_local = _size_global / __mpi_instance__.get_global_size();
    if (_size_global % __mpi_instance__.get_global_size() > __mpi_instance__.get_global_rank())
    {
        ++_size_local;
        _firstentrynumber = _size_local * __mpi_instance__.get_global_rank();
    }
    else
        _firstentrynumber = _size_local * __mpi_instance__.get_global_rank() + _size_global % __mpi_instance__.get_global_size();
}//Vector::Vector(const size_t size)

Vector& Vector::operator=(const Vector& other)
{
    if (this == &other) return *this;
    if (_data) delete[] _data;
    _size_global = other._size_global;
    _size_local = other._size_local;
    _firstentrynumber = other._firstentrynumber;
    _data = new double[_size_local];
    for (size_t i{0}; i < _size_local; ++i)
        _data[i] = other._data[i];
    return *this;
    //TODO error handling allocating
}//Vector& Vector::operator=(const Vector& other)

Vector& Vector::operator=(Vector&& other)
{
    if (this == &other) return *this;
    if (_data) delete[] _data;
    _size_global = other._size_global;
    _size_local = other._size_local;
    _firstentrynumber = other._firstentrynumber;
    _data = other._data;
    other._size_global = 0;
    other._size_local = 0;
    other._firstentrynumber = 0;
    other._data = nullptr;
    return *this;
}//Vector& Vector::operator=(Vector&& other)

Vector Vector::operator-() const
{
    Vector res(*this);
    for (size_t i{0}; i < _size_local; ++i)
        res._data[i] = - _data[i];
    return res;
}//Vector Vector::operator-() const

Vector& Vector::operator+(const Vector& other)
{
    assert(_size_global == other._size_global);
    assert(_size_local == other._size_local); // redundant
    Vector* res = new Vector(*this);
    *res += other;
    return *res;
}//Vector& Vector::operator+(const Vector& other)

Vector& Vector::operator+=(const Vector& other)
{
    assert(_size_global == other._size_global);
    assert(_size_local == other._size_local); // redundant
    for (size_t i{0}; i < _size_local; i++)
        _data[i] += other._data[i];
    return *this;
}//Vector& Vector::operator+=(const Vector& other)

Vector& Vector::operator-(const Vector& other)
{
    assert(_size_global == other._size_global);
    assert(_size_local == other._size_local); // redundant
    Vector* res = new Vector(*this);
    *res -= other;
    return *res;
}//Vector& Vector::operator-(const Vector& other)

Vector& Vector::operator-=(const Vector& other)
{
    assert(_size_global == other._size_global);
    assert(_size_local == other._size_local); // redundant
    for (size_t i{0}; i < _size_local; i++)
        _data[i] -= other._data[i];
    return *this;
}//Vector& Vector::operator-=(const Vector& other)

double Vector::get_global(const size_t pos) const
{
    assert(pos < _size_global);
    double res{0.0};
    if (pos >= _firstentrynumber && pos < _firstentrynumber + _size_local)
        res = _data[pos - _firstentrynumber];
    double globalres{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&res, &globalres, 1, MPI_DOUBLE, MPI_SUM);) //TODO should work with copying and not adding it up!
    //MPICALL(MPI::COMM_WORLD.Bcast(&res, 1, MPI_DOUBLE, __mpi_instance__.get_global_rank());) // somehow like this, I think
    return globalres;
}//double Vector::get_global(const size_t pos) const

double Vector::get_local(const size_t pos) const
{
    assert(pos < _size_local);
    return _data[pos];
}//double Vector::get_local(const size_t pos) const

void Vector::print_local()
{
    for (size_t i{0}; i < _size_local; ++i)
        std::cout << _data[i] << std::endl;
}//void Vector::print_local()

size_t Vector::get_firstentrynumber() const
{
    return _firstentrynumber;
}//size_t Vector::get_firstentrynumber() const

void Vector::set_global(const size_t pos, const double val)
{
    assert(pos < _size_global);
    if (pos >= _firstentrynumber && pos < _firstentrynumber + _size_local)
        _data[pos - _firstentrynumber] = val;
    MPICALL(MPI::COMM_WORLD.Barrier();) //TODISCUSS necessary?
}//void Vector::set_global(const size_t pos, const double val)

void Vector::set_local(const size_t pos, const double val)
{
    assert(pos < _size_local);
    _data[pos] = val;
}//void Vector::set_local(const size_t pos, const double val)

void Vector::add_global(const size_t pos, const double val)
{
    assert(pos < _size_global);
    if (pos >= _firstentrynumber && pos < _firstentrynumber + _size_local)
        _data[pos - _firstentrynumber] += val;
    MPICALL(MPI::COMM_WORLD.Barrier();) //TODISCUSS necessary?
}//void Vector::add_global(const size_t pos, const double val)

void Vector::add_local(const size_t pos, const double val)
{
    assert(pos < _size_local);
    _data[pos] += val;
}//void Vector::add_local(const size_t pos, const double val)

double Vector::l1norm() const
{
    double localsum{0.0};
    for (size_t i{0}; i < _size_local; ++i)
        localsum += std::abs(_data[i]);
    double globalsum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(globalsum);
}//double Vector::l1norm() const

double Vector::l2norm() const
{
    double localsum{0.0};
    for (size_t i{0}; i < _size_local; ++i)
        localsum += _data[i] * _data[i];
    double globalsum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(globalsum);
}//double Vector::l2norm() const

double Vector::lpnorm(const int p) const
{
    //TODO maybe later
    return lpnorm(static_cast<double>(p));
}//double Vector::lpnorm(const int p) const

double Vector::lpnorm(const double p) const
{
    double localsum{0.0};
    for (size_t i{0}; i < _size_local; ++i)
        localsum += std::pow(std::abs(_data[i]), p);
    double globalsum{0.0};
    MPICALL(MPI::COMM_WORLD.Barrier();)
    MPICALL(MPI::COMM_WORLD.Allreduce(&localsum, &globalsum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::pow(globalsum, 1.0/p);
}//double Vector::lpnorm(const double p) const

double Vector::maxnorm() const
{
    double res{0};
    //TODISCUSS which version is best?
    //for (size_t i{0}; i < _size_local; ++i)
    //    if (std::abs(_data[i]) > res)
    //        res = std::abs(_data[i]);
    //another possibility: std::for_each
    //for (double* ptr{_data}; ptr < _data + _size_local; ++ptr)          //TODISCUSS is this or
    //for (double* ptr{_data}, *upto{_data+_size_local}; ptr < upto; ++ptr) // this version better?
    //    if (std::abs(*ptr) > res)
    //        res = *ptr;
    //return res;
    double localmax{0.0};
    for (size_t i{0}; i < _size_local; ++i)
        if (std::abs(_data[i]) > localmax)
            localmax = std::abs(_data[i]);
    double globalmax{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&localmax, &globalmax, 1, MPI_DOUBLE, MPI_MAX);)
    return globalmax;
}//double Vector::maxnorm() const

}//namespace hptypes
