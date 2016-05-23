#ifndef __VECTOR_HPP_
#define __VECTOR_HPP_

#include <cstddef>
#include <cassert>
#include <iostream>

#include "../../helplib/include/helplib.hpp"

#include "matrix.hpp"

namespace hptypes
{

class Vector
{
public:
    Vector() = delete;
    Vector(const Vector& other);
    Vector(Vector&& other);
    Vector(size_t size);
    ~Vector() { if (_data) delete[] _data; }

    Vector& operator=(const Vector& other);
    Vector& operator=(Vector&& other);
    Vector operator-() const;
    Vector& operator+(const Vector& other);
    Vector& operator+=(const Vector& other);
    Vector& operator-(const Vector& other);
    Vector& operator-=(const Vector& other);
    // does not work as I hoped
    //Vector& operator=(const Vector&) = 0;
    //Vector& operator=(Vector&&) = 0;
    //Vector operator+(const Vector&) { throw not_implemented(); }
    //Vector& operator+=(const Vector&) { throw not_implemented(); }
    //Vector operator-(const Vector&) { throw not_implemented(); }
    //Vector& operator-=(const Vector&) { throw not_implemented(); }

    // blocking (deleting) some operators ?
    //bool operator/(const Vector&) = delete;
    //bool operator%(const Vector&) = delete;
    //bool operator!() = delete;
    //bool operator>(const Vector&) = delete;
    //bool operator/=(const Vector&) = delete;
    //bool operator%=(const Vector&) = delete;
    //bool operator<=(const Vector&) = delete;
    //bool operator>=(const Vector&) = delete;
    //bool operator++(const int) = delete;
    //bool operator--(const int) = delete;
    Vector& operator[](size_t) = delete;
    const Vector& operator[](size_t) const = delete;

    size_t get_size_global() const { return _size_global; }
    size_t get_size_local() const { return _size_local; }
    double get_global(size_t i) const { assert(i < _size_global); throw hphelp::not_implemented(); }
    double get_local(size_t i) const { assert(i < _size_local); return _data[i]; }
    // get_datasize_global() is not useful and is not to be implemented
    size_t get_datasize_local() const { return _size_local*sizeof(double) + 2*sizeof(size_t) + sizeof(double*); }
    void print_local();

    void set_global(size_t, double) { throw hphelp::not_implemented(); }
    void set_local(size_t, double);

    double l2norm();
    double lpnorm(int p);
    double lpnorm(double p);
    double maxnorm();

private:
    size_t _size_global;
    size_t _size_local;
    double* _data;
};//class vector

}//namespace hptypes

#endif//__VECTOR_HPP_
