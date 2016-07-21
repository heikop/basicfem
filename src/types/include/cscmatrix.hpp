#ifndef __CSCMATRIX_HPP_
#define __CSCMATRIX_HPP_

#include "matrix.hpp"

#include <vector>
#include <iostream> //TODO remove later

#include "../../mpihandler.hpp"
#include "../../helplib/include/helplib.hpp"

namespace hptypes
{

// this matrix has always unique entries (one data per matrixentry)
// and it is always sorted (rows and columns)
class CscMatrix : public Matrix<double>
{
public:
    CscMatrix() = delete;
    CscMatrix(const CscMatrix& other);
    CscMatrix(CscMatrix&& other);
    CscMatrix(const size_t numrows, const size_t numcols);
    ~CscMatrix();

    CscMatrix& operator=(const CscMatrix&);
    CscMatrix& operator=(CscMatrix&&);
    //CscMatrix operator+ (const CscMatrix&) { throw not_implemented(); }
    //CscMatrix& operator+= (const CscMatrix&) { throw not_implemented(); }
    //CscMatrix operator- (const CscMatrix&) { throw not_implemented(); }
    //CscMatrix& operator-= (const CscMatrix&) { throw not_implemented(); }
    bool operator==(const CscMatrix&);
    bool operator!=(const CscMatrix&);

    void free_unused_space();

    size_t get_numrows_global() const { return _numrows_global; }
    size_t get_numrows_local() const { return _numrows_local; }
    size_t get_numcols_global() const { return _numcols_global; }
    size_t get_numcols_local() const { return _numcols_local; }
    bool isquadratic() const { return _numrows_global == _numcols_global; }
    double get_global(const size_t, const size_t) const;
    double get_local(const size_t, const size_t) const;
    //size_t get_datasize_local() const { throw hphelp::not_implemented(); }
    //void print_local() { throw hphelp::not_implemented(); }

    void set_global(const size_t, const size_t, const double);
    void set_local(const size_t, const size_t, const double);

    //bool issymmetric() const { throw hphelp::not_implemented(); }

    double norm_1() const;
    double norm_2() const;
    double norm_inf() const;

    //CscMatrix& get_transpose() const { throw hphelp::not_implemented(); }
    //void transpose() { throw hphelp::not_implemented(); }
    //CscMatrix& get_mat_add(CscMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_add(CscMatrix& other) { throw hphelp::not_implemented(); }
    //CscMatrix& get_mat_sub(CscMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_sub(CscMatrix& other) { throw hphelp::not_implemented(); }
    CscMatrix& get_scal_mul(const double scal) const;
    void scal_mul(const double scal);
//    void pow(unsigned int exp);
    //CscMatrix& get_mat_mul(CscMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_mul(CscMatrix& other) { throw hphelp::not_implemented(); }
    //CscMatrix& get_inverse() const { throw hphelp::not_implemented(); }
    //void invert() { throw hphelp::not_implemented(); }

    DenseVector& get_vec_mul(const DenseVector& vec) const;
    DenseVector& get_pre_vec_mul(const DenseVector& vec) const;

private:
    size_t _numrows_global, _numcols_global;
    size_t _numrows_local, _numcols_local;
    size_t _firstrownum_globalcount;
    std::vector<double> _data;
    std::vector<size_t> _rowindex;
    std::vector<size_t> _firstcolentry;
};//class CscMatrix

}//namespace hptypes

#endif//ifndef __CSCMATRIX_HPP_
