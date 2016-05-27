#ifndef __COOMATRIX_HPP_
#define __COOMATRIX_HPP_

#include "matrix.hpp"

#include <vector>
#include <iostream> //TODO remove later

#include "../../mpihandler.hpp"
#include "../../helplib/include/helplib.hpp"

namespace hptypes
{

class CooMatrix : public Matrix<double>
{
public:
    CooMatrix() = delete;
    CooMatrix(const CooMatrix& other);
    CooMatrix(CooMatrix&& other);
    CooMatrix(size_t numrows, size_t numcols);
    ~CooMatrix();

    CooMatrix& operator= (const CooMatrix&);
    CooMatrix& operator= (CooMatrix&&);
    //CooMatrix operator+ (const CooMatrix&) { throw not_implemented(); }
    //CooMatrix& operator+= (const CooMatrix&) { throw not_implemented(); }
    //CooMatrix operator- (const CooMatrix&) { throw not_implemented(); }
    //CooMatrix& operator-= (const CooMatrix&) { throw not_implemented(); }

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

    //CooMatrix& get_transpose() const { throw hphelp::not_implemented(); }
    //void transpose() { throw hphelp::not_implemented(); }
    //CooMatrix& get_mat_add(CooMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_add(CooMatrix& other) { throw hphelp::not_implemented(); }
    //CooMatrix& get_mat_sub(CooMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_sub(CooMatrix& other) { throw hphelp::not_implemented(); }
    CooMatrix& get_scal_mul(double scal) const;
    void scal_mul(double scal);
//    void pow(unsigned int exp);
    //CooMatrix& get_mat_mul(CooMatrix& other) const { throw hphelp::not_implemented(); }
    //void mat_mul(CooMatrix& other) { throw hphelp::not_implemented(); }
    //CooMatrix& get_inverse() const { throw hphelp::not_implemented(); }
    //void invert() { throw hphelp::not_implemented(); }

private:
    bool _uniquemapping; // indicator whether there is max one _data entry per matrix entry or not
    size_t _numrows_global, _numcols_global;
    size_t _numrows_local, _numcols_local;
    size_t _firstrownumber;
    std::vector<double> _data;
    std::vector<size_t> _row;
    std::vector<size_t> _col;
};//class CooMatrix

}//namespace hptypes

#endif//ifndef __COOMATRIX_HPP_
