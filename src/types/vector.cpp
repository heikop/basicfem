#include "include/vector.hpp"

namespace hptypes
{

Vector::Vector(const Vector& other):
    _size_global{other._size_global}, _size_local{other._size_local},
    _data{new double[_size_local]}
{
    for (size_t i{0}; i < _size_local; ++i)
        _data[i] = other._data[i];
    //TODO error handling allocating
}//Vector::Vector(const Vector& other)

Vector::Vector(Vector&& other):
    _size_global{other._size_global}, _size_local{other._size_local},
    _data{other._data}
{
    other._size_global = 0;
    other._size_local = 0;
    other._data = nullptr;
}//Vector::Vector(Vector&& other)

Vector::Vector(size_t size):
    _size_global{size}, _size_local{size},
    _data{new double[_size_local]}
{
    for (size_t i{0}; i < _size_local; ++i)
        _data[i] = 0;
    //TODO error handling allocating
}//Vector::Vector(size_t size)

Vector& Vector::operator=(const Vector& other)
{
    if (this == &other) return *this;
    if (_data) delete[] _data;
    _size_global = other._size_global;
    _size_local = other._size_local;
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
    _data = other._data;
    other._size_global = 0;
    other._size_local = 0;
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
    assert(_size_global = other._size_global && _size_local == other._size_local); //TODO local difference...
    Vector* res = new Vector(*this);
    *res += other;
    return *res;
}//Vector& Vector::operator+(const Vector& other)

Vector& Vector::operator+=(const Vector& other)
{
    assert(_size_global = other._size_global && _size_local == other._size_local); //TODO local difference...
    for (size_t i{0}; i < _size_local; i++)
        _data[i] += other._data[i];
    return *this;
}//Vector& Vector::operator+=(const Vector& other)

Vector& Vector::operator-(const Vector& other)
{
    assert(_size_global = other._size_global && _size_local == other._size_local); //TODO local difference...
    Vector* res = new Vector(*this);
    *res -= other;
    return *res;
}//Vector& Vector::operator-(const Vector& other)

Vector& Vector::operator-=(const Vector& other)
{
    assert(_size_global = other._size_global && _size_local == other._size_local); //TODO local difference...
    for (size_t i{0}; i < _size_local; i++)
        _data[i] -= other._data[i];
    return *this;
}//Vector& Vector::operator-=(const Vector& other)

void Vector::print_local()
{
    for (size_t i{0}; i < _size_local; ++i)
        std::cout << _data[i] << std::endl;
}//void Vector::print_local()

void Vector::set_local(size_t i, double val)
{
    assert(i < _size_local);
    _data[i] = val;
}//void Vector::set_local(size_t i, double val)

double Vector::l2norm()
{
    double res{0};
    for (size_t i{0}; i < _size_local; ++i)
        res += _data[i] * _data[i];
    return std::sqrt(res);
}//double Vector::l2norm()

double Vector::lpnorm(int p)
{
    //TODO maybe later
    return lpnorm(static_cast<double>(p));
}//double Vector::lpnorm(int p)

double Vector::lpnorm(double p)
{
    double res{0};
    for (size_t i{0}; i < _size_local; ++i)
        res += std::pow(std::abs(_data[i]), p);
    return std::pow(res, static_cast<double>(1)/p);
}//double Vector::lpnorm(double p)

double Vector::maxnorm()
{
    double res{0};
    //TODO TODISCUSS which version is best?
    //for (size_t i{0}; i < _size_local; ++i)
    //    if (std::abs(_data[i]) > res)
    //        res = std::abs(_data[i]);
    //another possibility: std::for_each
    //for (double* ptr{_data}; ptr < _data + _size_local; ++ptr)          // TODO TODISCUSS is this or
    for (double* ptr{_data}, *upto{_data+_size_local}; ptr < upto; ++ptr) // this version better?
        if (std::abs(*ptr) > res)
            res = *ptr;
    return res;
}//double Vector::maxnorm()

}//namespace hptypes
