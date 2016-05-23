#ifndef __LILMATRIX_HPP_
#define __LILMATRIX_HPP_

#include "matrix.hpp"
#include "../../helplib/include/helplib.hpp"

namespace hptypes
{

class LilMatrix : public Matrix<double>
{
public:
    LilMatrix() = delete;
    LilMatrix(const LilMatrix& other);
    LilMatrix(LilMatrix&& other);
    LilMatrix(size_t numrows, size_t numcols);
    ~LilMatrix();

    LilMatrix& operator= (const LilMatrix&);
    LilMatrix& operator= (LilMatrix&&);
    //LilMatrix operator+ (const LilMatrix&) { throw not_implemented(); }
    //LilMatrix& operator+= (const LilMatrix&) { throw not_implemented(); }
    //LilMatrix operator- (const LilMatrix&) { throw not_implemented(); }
    //LilMatrix& operator-= (const LilMatrix&) { throw not_implemented(); }

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

    //LilMatrix& get_transpose() const { throw hphelp::not_implemented(); }
    //void transpose() { throw hphelp::not_implemented(); }
    //LilMatrix& get_mat_add(LilMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_add(LilMatrix& other) { throw hphelp::not_implemented(); }
    //LilMatrix& get_mat_sub(LilMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_sub(LilMatrix& other) { throw hphelp::not_implemented(); }
    LilMatrix& get_scal_mul(double scal) const;
    void scal_mul(double scal);
//    void pow(unsigned int exp);
    //LilMatrix& get_mat_mul(LilMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_mul(LilMatrix& other) { throw hphelp::not_implemented(); }
    //LilMatrix& get_inverse() const { throw hphelp::not_implemented(); }
    //void invert() { throw hphelp::not_implemented(); }

private:
    size_t _numrows_global, _numcols_global;
    size_t _numrows_local, _numcols_local;
};//class LilMatrix

}//namespace hptypes

#endif//ifndef __LILMATRIX_HPP_
