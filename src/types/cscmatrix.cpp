#include "include/cscmatrix.hpp"

namespace hptypes
{

CscMatrix::CscMatrix(const CscMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownum_globalcount{other._firstrownum_globalcount}
{
    _data = other._data;
    _rowindex = other._rowindex;
    _firstcolentry = other._firstcolentry;
}//CscMatrix::CscMatrix(const CscMatrix& other)

CscMatrix::CscMatrix(CscMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownum_globalcount{other._firstrownum_globalcount},
    _data{std::move(other._data)},
    _rowindex{std::move(other._rowindex)},
    _firstcolentry{std::move(other._firstcolentry)}
{
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownum_globalcount = 0;
}//CscMatrix::CscMatrix(CscMatrix&& other)

CscMatrix::CscMatrix(const size_t numrows, const size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{0}, _numcols_local{numcols},
    _data{std::vector<double>(0)},
    _rowindex{std::vector<size_t>(0)}
{
    _firstcolentry = std::vector<size_t>(_numcols_local + 1, 0);
    _numrows_local = _numrows_global / __mpi_instance__.get_global_size();
    if (_numrows_global % __mpi_instance__.get_global_size() > __mpi_instance__.get_global_rank())
    {
        ++_numrows_local;
        _firstrownum_globalcount = _numrows_local * __mpi_instance__.get_global_rank();
    }
    else
        _firstrownum_globalcount = _numrows_local * __mpi_instance__.get_global_rank() + _numrows_global % __mpi_instance__.get_global_size();
}//CscMatrix::CscMatrix(const size_t numrows, const size_t numcols)

CscMatrix::~CscMatrix()
{
}//CscMatrix::~CscMatrix()

CscMatrix& CscMatrix::operator=(const CscMatrix& other)
{
    if (this == &other) return *this;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownum_globalcount = other._firstrownum_globalcount;
    _data = other._data;
    _rowindex = other._rowindex;
    _firstcolentry = other._firstcolentry;
}//CscMatrix& CscMatrix::operator=(const CscMatrix& other)

CscMatrix& CscMatrix::operator=(CscMatrix&& other)
{
    if (this == &other) return *this;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownum_globalcount = other._firstrownum_globalcount;
    _data = std::move(other._data);
    _rowindex = std::move(other._rowindex);
    _firstcolentry = std::move(other._firstcolentry);
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownum_globalcount = 0;
}//CscMatrix& CscMatrix::operator=(CscMatrix&& other)

bool CscMatrix::operator==(const CscMatrix& other)
{
    if (*this == other) return true;
    if (_numrows_global != other._numcols_global
        || _numcols_global != other._numcols_global)
        return false;
    for (size_t i{0}; i < _data.size(); ++i)
        if (_data[i] != other._data[i])
            return false;
    return true;
}//bool CscMatrix::operator==(const CscMatrix& other)

bool CscMatrix::operator!=(const CscMatrix& other)
{
    return !operator==(other);
}//bool CscMatrix::operator!=(const CscMatrix& other)

void CscMatrix::free_unused_space()
{
    _data.shrink_to_fit();
    _rowindex.shrink_to_fit();
}//void CscMatrix::free_unused_space()

double CscMatrix::get_global(const size_t row, const size_t col) const
{
    assert(row < _numrows_global && col < _numcols_global);
    double val{0.0};
    if (row >= _firstrownum_globalcount && row < _firstrownum_globalcount + _numrows_local)
        val = get_local(row - _firstrownum_globalcount, col);
    double val_global{0.0};
    //MPICALL(MPI::COMM_WORLD.Barrier();)
    MPICALL(MPI::COMM_WORLD.Allreduce(&val, &val_global, 1, MPI_DOUBLE, MPI_SUM);) //TODO should work with copying and not adding it up!
    //MPICALL(MPI::COMM_WORLD.Bcast(&val, 1, MPI_DOUBLE, __mpi_instance__.get_global_rank());) // somehow like this, I think
    return val_global;
}//double CscMatrix::get_global(const size_t row, const size_t col) const

double CscMatrix::get_local(const size_t row, const size_t col) const
{
    assert(row < _numrows_local && col < _numcols_local);
    for (size_t datalocation{_firstcolentry[col]}; datalocation < _firstcolentry[col+1]; ++datalocation)
    {
        if (_rowindex[datalocation] == row)
            return _data[datalocation];
        else if (_rowindex[datalocation] > row)
            return 0.0;
    }
    return 0.0;
}//double CscMatrix::get_local(const size_t row, const size_t col) const

void CscMatrix::set_global(const size_t row, const size_t col, const double val)
{
    assert(row < _numrows_global && col < _numcols_global);
    if (row >= _firstrownum_globalcount && row < _firstrownum_globalcount + _numrows_local)
        set_local(row - _firstrownum_globalcount, col, val);
    MPICALL(MPI::COMM_WORLD.Barrier();) //TODISCUSS necessary?
}//void CscMatrix::set_global(const size_t row, const size_t col, const double val)

void CscMatrix::set_local(const size_t row, const size_t col, const double val)
{
    assert(row < _numrows_local && col < _numcols_local);
    size_t pos_to_insert{_firstcolentry[col]};
    bool new_entry{true};
    while (pos_to_insert < _firstcolentry[col+1])
    {
        if (_rowindex[pos_to_insert] < row)
            ++pos_to_insert;
        else if (_rowindex[pos_to_insert] == row)
        {
            _data[pos_to_insert] = val;
            return;
        }
    }
    _rowindex.insert(_rowindex.begin() + pos_to_insert, row);
    _data.insert(_data.begin() + pos_to_insert, val);
    for (size_t i{col + 1}; i <= _numcols_local; ++i)
        ++_firstcolentry[i]; //TOCHECK does this work?
        //_firstcolentry[i] += 1; // otherwise do this
}//void CscMatrix::set_local(const size_t row, const size_t col, const double val)

double CscMatrix::norm_1() const
{
    std::vector<double> localsums(_numcols_global, 0.0);
    for (size_t j{0}; j < _firstcolentry.size() - 1; ++j)
        for (size_t i{_firstcolentry[i]}; i < +_firstcolentry[i+1]; ++i)
            localsums[j] += std::abs(_data[i]);
    std::vector<double> globalsums(_numcols_global);
    MPICALL(MPI::COMM_WORLD.Allreduce(localsums.data(), globalsums.data(), localsums.size(), MPI_DOUBLE, MPI_SUM);)
    double res{0.0};
    for (const auto x : globalsums)
        if (x > res)
            res = x;
    return res;
}//double CscMatrix::norm_1() const

double CscMatrix::norm_2() const
{
    double local_squaresum{0.0};
    for (const auto val : _data)
        local_squaresum += val*val;
    double global_squaresum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&local_squaresum, &global_squaresum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(global_squaresum);
}//double CscMatrix::norm_2() const

double CscMatrix::norm_inf() const
{
    std::vector<double> rowsums(_numrows_local, 0.0);
    for (size_t i{0}; i < _data.size(); ++i)
        rowsums[_rowindex[i]] += _data[i];
    double localmax{0.0};
    for (const auto x : rowsums)
        if (x > localmax)
            localmax = x;
    double globalmax{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&localmax, &globalmax, 1, MPI_DOUBLE, MPI_MAX);)
    return globalmax;
}//double CscMatrix::norm_inf() const

CscMatrix& CscMatrix::get_scal_mul(const double scal) const
{
    CscMatrix* res = new CscMatrix(*this);
    if (scal == 1.0) return *res;
    for (auto& val : res->_data)
        val *= scal;
    return *res;
}//CscMatrix& CscMatrix::get_scal_mul(const double scal) const

void CscMatrix::scal_mul(const double scal)
{
    if (scal != 1.0)
        for (auto& val : _data)
            val *= scal;
}//void CscMatrix::scal_mul(const double scal)

DenseVector& CscMatrix::get_vec_mul(const DenseVector& vec) const
{
    assert(vec.get_size_global() == _numcols_global);
    DenseVector* res = new DenseVector(_numrows_global);
    for (size_t j{0}; j < _numcols_local; ++j)
        for (size_t i{_firstcolentry[j]}; i < _firstcolentry[j+1]; ++i)
            res->add_local(_rowindex[i], _data[i] * vec.get_local(j));
    return *res;
}//DenseVector& CscMatrix::get_vec_mul(const DenseVector& vec) const

DenseVector& CscMatrix::get_pre_vec_mul(const DenseVector& vec) const
{
    assert(vec.get_size_global() == _numrows_global);
    DenseVector* res = new DenseVector(_numcols_global);
    for (size_t j{0}; j < _numcols_local; ++j)
        for (size_t i{_firstcolentry[j]}; i < _firstcolentry[j+1]; ++i)
            res->add_local(j, vec.get_local(_rowindex[i]) * _data[i]);
    return *res;
}//DenseVector& CscMatrix::get_pre_vec_mul(const DenseVector& vec) const

}//namespace hptypes
