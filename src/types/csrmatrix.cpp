#include "include/csrmatrix.hpp"

namespace hptypes
{

CsrMatrix::CsrMatrix(const CsrMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

CsrMatrix::CsrMatrix(CsrMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

CsrMatrix::CsrMatrix(size_t numrows, size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{numrows}, _numcols_local{numcols} // TODO
{
    //TODO
}

CsrMatrix::~CsrMatrix()
{
    //TODO
}

CsrMatrix& CsrMatrix::operator=(const CsrMatrix& other)
{
    if (this == &other) return *this;
    //TODO
}

CsrMatrix& CsrMatrix::operator=(CsrMatrix&& other)
{
    if (this == &other) return *this;
    //TODO
}

double CsrMatrix::get_global(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

double CsrMatrix::get_local(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

void CsrMatrix::set_global(size_t row, size_t col, double val)
{
    //TODO
}

void CsrMatrix::set_local(size_t row, size_t col, double val)
{
    //TODO
}

double CsrMatrix::norm_1() const
{
    //TODO
    return 0.0;
}

double CsrMatrix::norm_2() const
{
    //TODO
    return 0.0;
}

double CsrMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}

CsrMatrix& CsrMatrix::get_scal_mul(double scal) const
{
    //TODO
    CsrMatrix* res = new CsrMatrix(*this);
    return *res;
}

void CsrMatrix::scal_mul(double scal)
{
    //TODO
}

}//namespace hptypes
