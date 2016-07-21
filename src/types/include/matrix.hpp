#ifndef __MATRIX_HPP_
#define __MATRIX_HPP_

#include <cstddef>
#include <cassert>

#include "../../helplib/include/helplib.hpp"

//#include "vector.hpp"
#include "densevector.hpp"
//#include "sparsevector.hpp"

namespace hptypes
{

template <typename datatype>
class Matrix
{
public:
    virtual ~Matrix() {}

    // does not work as I hoped
    //virtual Matrix& operator= (const Matrix&) = 0;
    //virtual Matrix& operator= (Matrix&&) = 0;
    //virtual Matrix operator+ (const Matrix&) { throw not_implemented(); }
    //virtual Matrix& operator+= (const Matrix&) { throw not_implemented(); }
    //virtual Matrix operator- (const Matrix&) { throw not_implemented(); }
    //virtual Matrix& operator-= (const Matrix&) { throw not_implemented(); }

    virtual size_t get_numrows_global() const = 0;
    virtual size_t get_numrows_local() const = 0;
    virtual size_t get_numcols_global() const = 0;
    virtual size_t get_numcols_local() const = 0;
    virtual bool isquadratic() const = 0;
    virtual datatype get_global(const size_t, const size_t) const = 0;
    virtual datatype get_local(const size_t, const size_t) const = 0;
    // get_datasize_global() is not useful and is not to be implemented
    virtual size_t get_datasize_local() const { throw hphelp::not_implemented(); }
    virtual void print_local() const { throw hphelp::not_implemented(); }

    // blocking (deleting) some operators ?
    //bool operator/ (const Matrix&) = delete;
    //bool operator% (const Matrix&) = delete;
    //bool operator! () = delete;
    //bool operator> (const Matrix&) = delete;
    //bool operator/= (const Matrix&) = delete;
    //bool operator%= (const Matrix&) = delete;
    //bool operator<= (const Matrix&) = delete;
    //bool operator>= (const Matrix&) = delete;
    //bool operator++ (const int) = delete;
    //bool operator-- (const int) = delete;
    Matrix& operator[] (size_t) = delete;
    const Matrix& operator[] (size_t) const = delete;

    virtual void set_global(const size_t, const size_t, datatype) = 0;
    virtual void set_local(const size_t, const size_t, datatype) = 0;

    virtual bool issymmetric() const { throw hphelp::not_implemented(); }

    virtual double norm_1() const = 0;
    virtual double norm_2() const = 0;
    virtual double norm_inf() const = 0;

    virtual Matrix& get_transpose() const { throw hphelp::not_implemented(); }
    virtual void transpose() { throw hphelp::not_implemented(); }
    virtual Matrix& get_mat_add(const Matrix& other) const { throw hphelp::not_implemented(); }
    virtual void mat_add(const Matrix& other) { throw hphelp::not_implemented(); }
    virtual Matrix& get_mat_sub(const Matrix& other) const { throw hphelp::not_implemented(); }
    virtual void mat_sub(const Matrix& other) { throw hphelp::not_implemented(); }
    virtual Matrix& get_scal_mul(const double scal) const = 0;
    virtual void scal_mul(const double scal) = 0;
//    virtual void pow(unsigned int exp);
    virtual Matrix& get_mat_mul(const Matrix& other) const { throw hphelp::not_implemented(); }
    virtual void mat_mul(const Matrix& other) { throw hphelp::not_implemented(); }
    virtual Matrix& get_inverse() const { throw hphelp::not_implemented(); }
    virtual void invert() { throw hphelp::not_implemented(); }

    virtual DenseVector& get_vec_mul(const DenseVector& vec) const { throw hphelp::not_implemented(); } //TODO later this should be "= 0"
    virtual DenseVector& get_pre_vec_mul(const DenseVector& vec) const { throw hphelp::not_implemented(); } //TODO maybe change later to "= 0"
};//class Matrix

}//namespace hptypes

#endif//ifndef __MATRIX_HPP_
