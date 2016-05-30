#ifndef __DENSEMATRIX_HPP_
#define __DENSEMATRIX_HPP_

#include "matrix.hpp"
#include "../../helplib/include/helplib.hpp"

#include <iostream>
#include <algorithm>

namespace hptypes
{

class DenseMatrix : public Matrix<double>
{
public:
    DenseMatrix() = delete;
    DenseMatrix(const DenseMatrix& other);
    DenseMatrix(DenseMatrix&& other);
    DenseMatrix(size_t numrows, size_t numcols);
    ~DenseMatrix() { delete[] _data; delete[] _rowptrs; }

    DenseMatrix& operator= (const DenseMatrix& other);
    DenseMatrix& operator= (DenseMatrix&& other);

    size_t get_numrows_global() const { return _numrows_global; }
    size_t get_numrows_local() const { return _numrows_local; }
    size_t get_numcols_global() const { return _numcols_global; }
    size_t get_numcols_local() const { return _numcols_local; }
    bool isquadratic() const { return _numrows_global == _numcols_global; }
    double get_global(const size_t row, const size_t col) const { assert(row < _numrows_global && col < _numcols_global); throw hphelp::not_implemented(); }
    double get_local(const size_t row, const size_t col) const { assert(row < _numrows_local && col < _numcols_local); return _rowptrs[row][col]; }
    size_t get_datasize_local() const { return _numrows_local*_numcols_local*sizeof(double) + 4*sizeof(size_t)+ sizeof(double**) + sizeof(double*); }
    void print_local() const;

    void set_global(const size_t row, const size_t col, const double val) { throw hphelp::not_implemented(); }
    void set_local(const size_t row, const size_t col, const double val) { _rowptrs[row][col] = val; }

    bool issymmetric() const;

    double norm_1() const;
    double norm_2() const;
    double norm_inf() const;

    DenseMatrix& get_transpose() const;
    void transpose();
    DenseMatrix& get_mat_add(const DenseMatrix& other) const;
    void mat_add(const DenseMatrix& other);
    DenseMatrix& get_mat_sub(const DenseMatrix& other) const;
    void mat_sub(const DenseMatrix& other);
    DenseMatrix& get_scal_mul(const double scal) const;
    void scal_mul(const double scal);
//    void pow(unsigned int exp);
    DenseMatrix& get_mat_mul(const DenseMatrix& other) const;
    void mat_mul(const DenseMatrix& other);
    DenseMatrix& get_inverse() const;
    void invert();

private:
    size_t _numrows_global, _numcols_global;
    size_t _numrows_local, _numcols_local;
    double** _rowptrs;
    double* _data;
};//class DenseMatrix

}//namespace hptypes

#endif//ifndef __DENSEMATRIX_HPP_
