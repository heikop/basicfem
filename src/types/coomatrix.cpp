#include "include/coomatrix.hpp"

namespace hptypes
{

CooMatrix::CooMatrix(const CooMatrix& other):
    _uniquemapping{other._uniquemapping}, _sorted{other._sorted},
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
    _uniquemapping{other._uniquemapping}, _sorted{other._sorted},
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownumber{other._firstrownumber},
    _data{std::move(other._data)},
    _row{std::move(other._row)}, _col{std::move(other._col)}
{
    other._uniquemapping = true;
    other._sorted = true;
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
}//CooMatrix::CooMatrix(CooMatrix&& other):

CooMatrix::CooMatrix(const size_t numrows, const size_t numcols):
    _uniquemapping{true}, _sorted{true},
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
    _sorted = other._sorted;
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
    _sorted = other._sorted;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownumber = other._firstrownumber;
    _data = std::move(other._data);
    _row = std::move(other._row);
    _col = std::move(other._col);
    other._uniquemapping = true;
    other._sorted = true;
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownumber = 0;
}//CooMatrix& CooMatrix::operator=(CooMatrix&& other)

bool CooMatrix::operator==(const CooMatrix& other)
{
    if (*this == other) return true;
    if (_numrows_global != other._numcols_global
        || _numcols_global != other._numcols_global)
        return false;
    for (size_t i{0}; i < _numrows_local; ++i)
        for (size_t j{0}; j < _numcols_local; ++j)
            if (get_local(i, j) != other.get_local(i, j))
                return false;
    return true;
}//bool CooMatrix::operator==(const CooMatrix& other)

bool CooMatrix::operator!=(const CooMatrix& other)
{
    return !operator==(other);
}//bool CooMatrix::operator!=(const CooMatrix& other)

void CooMatrix::make_unique()
{
    if (_uniquemapping) return;
    if (_sorted)
    {
        for (size_t i{0}; i < _data.size(); ++i)
            for (size_t j{i+1}; j < _data.size(); ++j)
            {
                if (_row[i] == _row[j] && _col[i] == _col[j])
                {
                    _data[i] += _data[j];
                    _data.erase(_data.begin() + i);
                    _row.erase(_row.begin() + i);
                    _col.erase(_col.begin() + i);
                }
                else
                    j == _data.size();
            }
    }
    else
    {
        for (size_t i{0}; i < _data.size(); ++i)
            for (size_t j{i+1}; j < _data.size(); ++j)
                if (_row[i] == _row[j] && _col[i] == _col[j])
                {
                    _data[i] += _data[j];
                    _data.erase(_data.begin() + i);
                    _row.erase(_row.begin() + i);
                    _col.erase(_col.begin() + i);
                }
    }
    _uniquemapping = true;
}

void CooMatrix::sort()
{
    //TODO
    if (_uniquemapping)
    {
    }
    else
    {
    }
    _sorted = true;
}

void CooMatrix::make_unique_and_sort()
{
    if (_uniquemapping)
        sort();
    else if (_sorted)
        make_unique();
    else
    {
        //TODO this probably could be done better
        sort();
        make_unique();
    }
    _uniquemapping = true;
    _sorted = true;
}

void CooMatrix::free_unused_space()
{
    _data.shrink_to_fit();
    _row.shrink_to_fit();
    _col.shrink_to_fit();
}//void CooMatrix::free_unused_space()

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
    if (_sorted)
    {
        if (_uniquemapping)
        {
        //TODO implement better finding algorithm 
        //this is just a copy of the unsorted case
            for (size_t i{0}; i < _row.size(); ++i)
                if (_row[i] == row && _col[i] == col)
                    return _data[i];
            return 0.0;
        //    size_t row_lowerbound{0}, row_upperbound{_numrows_local};
        //    while (row_lowerbound < row_upperbound)
        //    {
        //        if (
        //    }
        }
        else
        {
        //TODO implement better finding algorithm 
        //this is just a copy of the unsorted case
            double res{0.0};
            for (size_t i{0}; i < _row.size(); ++i)
                if (_row[i] == row && _col[i] == col)
                    res += _data[i];
            return res;
        }
    }
    else
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
}//double CooMatrix::get_local(size_t row, size_t col) const

void CooMatrix::set_global(size_t row, size_t col, double val)
{
    assert(row < _numrows_global && col < _numcols_global);
    if (row >= _firstrownumber && row < _firstrownumber + _numrows_local)
        set_local(row - _firstrownumber, col, val);
    MPICALL(MPI::COMM_WORLD.Barrier();) //TODISCUSS necessary?
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
        _sorted = false;
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
        _sorted = false;
    }
}//void CooMatrix::set_local(size_t row, size_t col, double val)

double CooMatrix::norm_1() const
{
    //TODO this is not the columnsum norm it is ment to be...
    if (!_uniquemapping)
    {
        //TODO without not possible
        assert(false);
    }
    std::vector<double> localsums(_numcols_global, 0.0);
    for (size_t i{0}; i < _data.size(); ++i)
        localsums[_col[i]] += std::abs(_data[i]);
    std::vector<double> globalsums(_numcols_global);
    MPICALL(MPI::COMM_WORLD.Allreduce(localsums.data(), globalsums.data(), localsums.size(), MPI_DOUBLE, MPI_SUM);)
    double res{0.0};
    for (const auto x : globalsums)
        if (x > res)
            res = x;
    return res;
}//double CooMatrix::norm_1() const

double CooMatrix::norm_2() const
{
    if (!_uniquemapping)
    {
        //TODO without not possible
        assert(false);
    }
    double local_squaresum{0.0};
    for (const auto val : _data)
        local_squaresum += val*val;
    double global_squaresum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&local_squaresum, &global_squaresum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(global_squaresum);
}//double CooMatrix::norm_2() const

double CooMatrix::norm_inf() const
{
    if (!_uniquemapping)
    {
        //TODO without not possible
        assert(false);
    }
    std::vector<double> rowsums(_numrows_local, 0.0);
    for (size_t i{0}; i < _data.size(); ++i)
        rowsums[_row[i]] += std::abs(_data[i]);
    double localnorm{0.0};
    for (const auto x : rowsums)
        if (x > localnorm)
            localnorm = x;
    double globalnorm{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&localnorm, &globalnorm, 1, MPI_DOUBLE, MPI_MAX);)
    return globalnorm;
}//double CooMatrix::norm_inf() const

CooMatrix& CooMatrix::get_scal_mul(const double scal) const
{
    CooMatrix* res = new CooMatrix(*this);
    if (scal == 1.0) return *res;
    for (auto& val : res->_data)
        val *= scal;
    return *res;
}//CooMatrix& CooMatrix::get_scal_mul(const double scal) const

void CooMatrix::scal_mul(const double scal)
{
    if (scal != 1.0)
        for (auto& val : _data)
            val *= scal;
}//void CooMatrix::scal_mul(const double scal)

DenseVector& CooMatrix::get_vec_mul(const DenseVector& vec) const
{
    //TODISCUSS is it faster to copy the full vector to every node? (-> to avoid get_global)
    assert(_numcols_global == vec.get_size_global());
    DenseVector* res = new DenseVector(vec.get_size_global());
    for (size_t i{0}; i < _data.size(); ++i)
        res->add_local(_row[i] - _firstrownumber, _data[i] * vec.get_global(_col[i]));
    return *res;
}//DenseVector& CooMatrix::get_vec_mul(const DenseVector& vec) const

DenseVector& CooMatrix::get_pre_vec_mul(const DenseVector& vec) const
{
    assert(_numrows_global == vec.get_size_global());
    std::vector<double> localres(_numcols_global, 0.0);
    for (size_t i{0}; i < _data.size(); ++i)
        localres[_col[i]] += vec.get_local(_row[i]) * _data[i];
    std::vector<double> globalres(_numcols_global);
    MPICALL(MPI::COMM_WORLD.Allreduce(localres.data(), globalres.data(), _numcols_global, MPI_DOUBLE, MPI_SUM);)
    DenseVector* res = new DenseVector(vec.get_size_global()); //TODO use std::vector constructor for Vectcor (when available)
    for (size_t locind{0}, globind{res->get_firstentrynumber()}; locind < res->get_size_local(); ++locind, ++globind)
        res->set_local(locind, globalres[globind]);
    return *res;
}//DenseVector& CooMatrix::get_pre_vec_mul(const DenseVector& vec) const

}//namespace hptypes
