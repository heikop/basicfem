#include "include/coomatrix.hpp"

namespace hptypes
{

CooMatrix::CooMatrix(const CooMatrix& other):
    _uniquemapping{other._uniquemapping},
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownumber{other._firstrownumber}
{
    try
    {
        _data = other._data;
        _row = other._row;
        _col = other._col;
    }
    catch (...)
    {
        //TODO
    }
}//CooMatrix::CooMatrix(const CooMatrix& other)

CooMatrix::CooMatrix(CooMatrix&& other):
    _uniquemapping{other._uniquemapping},
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownumber{other._firstrownumber},
    _data{std::move(other._data)},
    _row{std::move(other._row)}, _col{std::move(other._col)}
{
}//CooMatrix::CooMatrix(CooMatrix&& other):

CooMatrix::CooMatrix(size_t numrows, size_t numcols):
    _uniquemapping{true},
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{0}, _numcols_local{numcols}
{
    std::vector<double>().swap(_data);
    std::vector<size_t>().swap(_row);
    std::vector<size_t>().swap(_col);
    _numrows_local = _numrows_global / __mpi_instance__.get_global_size();
    if (_numrows_global % __mpi_instance__.get_global_size() > __mpi_instance__.get_global_rank())
    {
        ++_numrows_local;
        _firstrownumber = _numrows_local * __mpi_instance__.get_global_rank();
    }
    else
        _firstrownumber = _numrows_local * __mpi_instance__.get_global_rank() + _numrows_global % __mpi_instance__.get_global_size();
}//CooMatrix::CooMatrix(size_t numrows, size_t numcols):

CooMatrix::~CooMatrix()
{
}//CooMatrix::~CooMatrix()

CooMatrix& CooMatrix::operator=(const CooMatrix& other)
{
    if (this == &other) return *this;
    _uniquemapping = other._uniquemapping;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownumber = other._firstrownumber;
    _data = other._data;
    _row = other._row;
    _col = other._col;
}//CooMatrix& CooMatrix::operator=(const CooMatrix& other)

CooMatrix& CooMatrix::operator=(CooMatrix&& other)
{
    if (this == &other) return *this;
    std::vector<double>().swap(_data);
    std::vector<size_t>().swap(_row);
    std::vector<size_t>().swap(_col);
    _uniquemapping = other._uniquemapping;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownumber = other._firstrownumber;
    _data = std::move(other._data);
    _row = std::move(other._row);
    _col = std::move(other._col);
    other._uniquemapping = true;
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownumber = 0;
}//CooMatrix& CooMatrix::operator=(CooMatrix&& other)

double CooMatrix::get_global(size_t row, size_t col) const
{
    assert(row < _numrows_global && col < _numcols_global);
    double val{0.0};
    if (row >= _firstrownumber && row < _firstrownumber + _numrows_local)
        val = get_local(row - _firstrownumber, col);
    double val_global{0.0};
    //MPICALL(MPI::COMM_WORLD.Barrier();)
    MPICALL(MPI::COMM_WORLD.Allreduce(&val, &val_global, 1, MPI_DOUBLE, MPI_SUM);) //TODO should work with copying and not adding it up!
    //MPICALL(MPI::COMM_WORLD.Bcast(&val, 1, MPI_DOUBLE, __mpi_instance__.get_global_rank());) // somehow like this, I think
    return val_global;
}//double CooMatrix::get_global(size_t row, size_t col) const

double CooMatrix::get_local(size_t row, size_t col) const
{
    assert(row < _numrows_local && col < _numcols_local);
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
}//double CooMatrix::get_local(size_t row, size_t col) const

void CooMatrix::set_global(size_t row, size_t col, double val)
{
    assert(row < _numrows_global && col < _numcols_global);
    if (row >= _firstrownumber && row < _firstrownumber + _numrows_local)
        set_local(row - _firstrownumber, col, val);
}//void CooMatrix::set_global(size_t row, size_t col, double val)

void CooMatrix::set_local(size_t row, size_t col, double val)
{
    assert(row < _numrows_local && col < _numcols_local);
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
}//void CooMatrix::set_local(size_t row, size_t col, double val)

double CooMatrix::norm_1() const
{
    //TODO
    return 0.0;
}//double CooMatrix::norm_1() const

double CooMatrix::norm_2() const
{
    double local_squaresum{0.0};

    if (_uniquemapping)
    {
        for (auto val : _data)
            local_squaresum += val*val;
    }
    else
    {
        //TODO
    }

    double global_squaresum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&local_squaresum, &global_squaresum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(global_squaresum);
}//double CooMatrix::norm_2() const

double CooMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}//double CooMatrix::norm_inf() const

CooMatrix& CooMatrix::get_scal_mul(double scal) const
{
    CooMatrix* res = new CooMatrix(*this);
    if (scal == 1.0) return *res;
    for (auto& val : res->_data)
        val *= scal;
    return *res;
}//CooMatrix& CooMatrix::get_scal_mul(double scal) const

void CooMatrix::scal_mul(double scal)
{
    if (scal != 1.0)
        for (auto& val : _data)
            val *= scal;
}//void CooMatrix::scal_mul(double scal)

}//namespace hptypes
