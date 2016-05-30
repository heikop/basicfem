#ifndef __VECTOR_HPP_
#define __VECTOR_HPP_

#include <cstddef>
#include <cassert>
#include <iostream> //TODO remove later

#include "../../mpihandler.hpp"
#include "../../helplib/include/helplib.hpp"

//#include "matrix.hpp"

namespace hptypes
{

class Vector
{
public:
    Vector() = delete;
    Vector(const Vector& other);
    Vector(Vector&& other);
    Vector(const size_t size);
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
    double get_global(const size_t i) const;
    double get_local(const size_t i) const;
    size_t get_datasize_local() const { return _size_local*sizeof(double) + 2*sizeof(size_t) + sizeof(double*); }
    void print_local();
    size_t get_firstentrynumber() const;

    void set_global(const size_t, const double);
    void set_local(const size_t, const double);
    void add_global(const size_t, const double);
    void add_local(const size_t, const double);

    double l1norm() const;
    double l2norm() const;
    double lpnorm(const int p) const;
    double lpnorm(const double p) const;
    double maxnorm() const;

private:
    size_t _size_global;
    size_t _size_local;
    size_t _firstentrynumber;
    double* _data;
};//class vector

}//namespace hptypes

#endif//ifndef __VECTOR_HPP_
