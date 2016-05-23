#include "include/cscmatrix.hpp"

namespace hptypes
{

CscMatrix::CscMatrix(const CscMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

CscMatrix::CscMatrix(CscMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

CscMatrix::CscMatrix(size_t numrows, size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{numrows}, _numcols_local{numcols} // TODO
{
    //TODO
}

CscMatrix::~CscMatrix()
{
    //TODO
}

CscMatrix& CscMatrix::operator=(const CscMatrix& other)
{
    if (this == &other) return *this;
    //TODO
}

CscMatrix& CscMatrix::operator=(CscMatrix&& other)
{
    if (this == &other) return *this;
    //TODO
}

double CscMatrix::get_global(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

double CscMatrix::get_local(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

void CscMatrix::set_global(size_t row, size_t col, double val)
{
    //TODO
}

void CscMatrix::set_local(size_t row, size_t col, double val)
{
    //TODO
}

double CscMatrix::norm_1() const
{
    //TODO
    return 0.0;
}

double CscMatrix::norm_2() const
{
    //TODO
    return 0.0;
}

double CscMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}

CscMatrix& CscMatrix::get_scal_mul(double scal) const
{
    //TODO
    CscMatrix* res = new CscMatrix(*this);
    return *res;
}

void CscMatrix::scal_mul(double scal)
{
    //TODO
}

}//namespace hptypes
