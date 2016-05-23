#ifndef __CSRMATRIX_HPP_
#define __CSRMATRIX_HPP_

#include "matrix.hpp"
#include "../../helplib/include/helplib.hpp"

namespace hptypes
{

class CsrMatrix : public Matrix<double>
{
public:
    CsrMatrix() = delete;
    CsrMatrix(const CsrMatrix& other);
    CsrMatrix(CsrMatrix&& other);
    CsrMatrix(size_t numrows, size_t numcols);
    ~CsrMatrix();

    CsrMatrix& operator= (const CsrMatrix&);
    CsrMatrix& operator= (CsrMatrix&&);
    //CsrMatrix operator+ (const CsrMatrix&) { throw not_implemented(); }
    //CsrMatrix& operator+= (const CsrMatrix&) { throw not_implemented(); }
    //CsrMatrix operator- (const CsrMatrix&) { throw not_implemented(); }
    //CsrMatrix& operator-= (const CsrMatrix&) { throw not_implemented(); }

    size_t get_numrows_global() const { return _numrows_global; }
    size_t get_numrows_local() const { return _numrows_local; }
    size_t get_numcols_global() const { return _numcols_global; }
    size_t get_numcols_local() const { return _numcols_local; }
    bool isquadratic() const { return _numrows_global == _numcols_global; }
    double get_global(size_t, size_t) const;
    double get_local(size_t, size_t) const;
    //size_t get_datasize_local() const { throw hphelp::not_implemented(); }
    //void print_local() { throw hphelp::not_implemented(); }

    void set_global(size_t, size_t, double);
    void set_local(size_t, size_t, double);

    //bool issymmetric() const { throw hphelp::not_implemented(); }

    double norm_1() const;
    double norm_2() const;
    double norm_inf() const;

    //CsrMatrix& get_transpose() const { throw hphelp::not_implemented(); }
    //void transpose() { throw hphelp::not_implemented(); }
    //CsrMatrix& get_mat_add(CsrMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_add(CsrMatrix& other) { throw hphelp::not_implemented(); }
    //CsrMatrix& get_mat_sub(CsrMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_sub(CsrMatrix& other) { throw hphelp::not_implemented(); }
    CsrMatrix& get_scal_mul(double scal) const;
    void scal_mul(double scal);
//    void pow(unsigned int exp);
    //CsrMatrix& get_mat_mul(CsrMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_mul(CsrMatrix& other) { throw hphelp::not_implemented(); }
    //CsrMatrix& get_inverse() const { throw hphelp::not_implemented(); }
    //void invert() { throw hphelp::not_implemented(); }

private:
    size_t _numrows_global, _numcols_global;
    size_t _numrows_local, _numcols_local;
};//class CsrMatrix

}//namespace hptypes

#endif//__CSRMATRIX_HPP_
