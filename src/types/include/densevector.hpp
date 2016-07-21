#ifndef __DENSEVECTOR_HPP_
#define __DENSEVECTOR_HPP_

#include <cstddef>
#include <cassert>
#include <iostream> //TODO remove later

#include "../../mpihandler.hpp"
#include "../../helplib/include/helplib.hpp"

//#include "matrix.hpp"

namespace hptypes
{

class DenseVector
{
public:
    DenseVector() = delete;
    DenseVector(const DenseVector& other);
    DenseVector(DenseVector&& other);
    DenseVector(const size_t size);
    ~DenseVector() { if (_data) delete[] _data; }

    DenseVector& operator=(const DenseVector& other);
    DenseVector& operator=(DenseVector&& other);
    DenseVector operator-() const;
    DenseVector& operator+(const DenseVector& other);
    DenseVector& operator+=(const DenseVector& other);
    DenseVector& operator-(const DenseVector& other);
    DenseVector& operator-=(const DenseVector& other);
    // does not work as I hoped
    //DenseVector& operator=(const DenseVector&) = 0;
    //DenseVector& operator=(DenseVector&&) = 0;
    //DenseVector operator+(const DenseVector&) { throw not_implemented(); }
    //DenseVector& operator+=(const DenseVector&) { throw not_implemented(); }
    //DenseVector operator-(const DenseVector&) { throw not_implemented(); }
    //DenseVector& operator-=(const DenseVector&) { throw not_implemented(); }

    // blocking (deleting) some operators ?
    //bool operator/(const DenseVector&) = delete;
    //bool operator%(const DenseVector&) = delete;
    //bool operator!() = delete;
    //bool operator>(const DenseVector&) = delete;
    //bool operator/=(const DenseVector&) = delete;
    //bool operator%=(const DenseVector&) = delete;
    //bool operator<=(const DenseVector&) = delete;
    //bool operator>=(const DenseVector&) = delete;
    //bool operator++(const int) = delete;
    //bool operator--(const int) = delete;
    DenseVector& operator[](size_t) = delete;
    const DenseVector& operator[](size_t) const = delete;

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
};//class DenseVector

}//namespace hptypes

#endif//ifndef __DENSEVECTOR_HPP_
