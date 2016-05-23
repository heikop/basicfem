#include "include/coomatrix.hpp"

namespace hptypes
{

CooMatrix::CooMatrix(const CooMatrix& other):
    _uniquemapping{other._uniquemapping},
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    _data = other._data;
    _row = other._row;
    _col = other._col;
}

CooMatrix::CooMatrix(CooMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local}
{
    //TODO
}

CooMatrix::CooMatrix(size_t numrows, size_t numcols):
    _uniquemapping{true},
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{numrows}, _numcols_local{numcols} // TODO
{
    std::vector<double>().swap(_data);
    std::vector<size_t>().swap(_row);
    std::vector<size_t>().swap(_col);
}

CooMatrix::~CooMatrix()
{
    //TODO
}

CooMatrix& CooMatrix::operator=(const CooMatrix& other)
{
    if (this == &other) return *this;
    _uniquemapping = other._uniquemapping;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _data = other._data;
    _row = other._row;
    _col = other._col;
}

CooMatrix& CooMatrix::operator=(CooMatrix&& other)
{
    if (this == &other) return *this;
    //TODO
}

double CooMatrix::get_global(size_t row, size_t col) const
{
    return get_local(row, col);
    //TODO
}

double CooMatrix::get_local(size_t row, size_t col) const
{
    if (_uniquemapping)
    {
        for (size_t i{0}; i < _row.size(); ++i)
            if (_row[i] == row && _col[i] == col)
                return _data[i];
        return 0.0;
    }
    else
    {
        double res{0.0};
        for (size_t i{0}; i < _row.size(); ++i)
            if (_row[i] == row && _col[i] == col)
                res += _data[i];
        return res;
    }
}

void CooMatrix::set_global(size_t row, size_t col, double val)
{
    set_local(row, col, val);
    //TODO
}

void CooMatrix::set_local(size_t row, size_t col, double val)
{
    if (_uniquemapping)
    {
        for (size_t i{0}; i < _row.size(); ++i)
            if (_row[i] == row && _col[i] == col)
            {
                _data[i] = val;
                return;
            }
        _data.push_back(val);
        _row.push_back(row);
        _col.push_back(col);
    }
    else
    {
        for (size_t i{0}; i < _row.size(); ++i)
            if (_row[i] == row && _col[i] == col)
            {
                _data.erase(_data.begin() + i);
                _row.erase(_row.begin() + i);
                _col.erase(_col.begin() + i);
            }
        _data.push_back(val);
        _row.push_back(row);
        _col.push_back(col);
    }
}

double CooMatrix::norm_1() const
{
    //TODO
    return 0.0;
}

double CooMatrix::norm_2() const
{
    //TODO
    return 0.0;
}

double CooMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}

CooMatrix& CooMatrix::get_scal_mul(double scal) const
{
    //TODO
    CooMatrix* res = new CooMatrix(*this);
    return *res;
}

void CooMatrix::scal_mul(double scal)
{
    //TODO
}

}//namespace hptypes
