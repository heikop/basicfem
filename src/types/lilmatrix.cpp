#include "include/lilmatrix.hpp"

namespace hptypes
{

LilMatrix::LilMatrix(const LilMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

LilMatrix::LilMatrix(LilMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

LilMatrix::LilMatrix(size_t numrows, size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{numrows}, _numcols_local{numcols} // TODO
{
    //TODO
}

LilMatrix::~LilMatrix()
{
    //TODO
}

LilMatrix& LilMatrix::operator=(const LilMatrix& other)
{
    if (this == &other) return *this;
    //TODO
}

LilMatrix& LilMatrix::operator=(LilMatrix&& other)
{
    if (this == &other) return *this;
    //TODO
}

double LilMatrix::get_global(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

double LilMatrix::get_local(size_t row, size_t col) const
{
    //TODO
    return 0.0;
}

void LilMatrix::set_global(size_t row, size_t col, double val)
{
    //TODO
}

void LilMatrix::set_local(size_t row, size_t col, double val)
{
    //TODO
}

double LilMatrix::norm_1() const
{
    //TODO
    return 0.0;
}

double LilMatrix::norm_2() const
{
    //TODO
    return 0.0;
}

double LilMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}

LilMatrix& LilMatrix::get_scal_mul(double scal) const
{
    //TODO
    LilMatrix* res = new LilMatrix(*this);
    return *res;
}

void LilMatrix::scal_mul(double scal)
{
    //TODO
}

}//namespace hptypes
