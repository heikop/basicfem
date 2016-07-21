#include "include/csrmatrix.hpp"

namespace hptypes
{

CsrMatrix::CsrMatrix(const CsrMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownum_globalcount{other._firstrownum_globalcount}
{
    _data = other._data;
    _colindex = other._colindex;
    _firstrowentry = other._firstrowentry;
}//CsrMatrix::CsrMatrix(const CsrMatrix& other)

CsrMatrix::CsrMatrix(CsrMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _firstrownum_globalcount{other._firstrownum_globalcount},
    _data{std::move(other._data)},
    _colindex{std::move(other._colindex)},
    _firstrowentry{std::move(other._firstrowentry)}
{
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownum_globalcount = 0;
}//CsrMatrix::CsrMatrix(CsrMatrix&& other)

CsrMatrix::CsrMatrix(const size_t numrows, const size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{0}, _numcols_local{numcols},
    _data{std::vector<double>(0)},
    _colindex{std::vector<size_t>(0)}
{
    _numrows_local = _numrows_global / __mpi_instance__.get_global_size();
    if (_numrows_global % __mpi_instance__.get_global_size() > __mpi_instance__.get_global_rank())
    {
        ++_numrows_local;
        _firstrownum_globalcount = _numrows_local * __mpi_instance__.get_global_rank();
    }
    else
        _firstrownum_globalcount = _numrows_local * __mpi_instance__.get_global_rank() + _numrows_global % __mpi_instance__.get_global_size();
    _firstrowentry = std::vector<size_t>(_numrows_local + 1, 0);
}//CsrMatrix::CsrMatrix(const size_t numrows, const size_t numcols)

CsrMatrix::~CsrMatrix()
{
}//CsrMatrix::~CsrMatrix()

CsrMatrix& CsrMatrix::operator=(const CsrMatrix& other)
{
    if (this == &other) return *this;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownum_globalcount = other._firstrownum_globalcount;
    _data = other._data;
    _colindex = other._colindex;
    _firstrowentry = other._firstrowentry;
}//CsrMatrix& CsrMatrix::operator=(const CsrMatrix& other)

CsrMatrix& CsrMatrix::operator=(CsrMatrix&& other)
{
    if (this == &other) return *this;
    _numrows_global = other._numrows_global;
    _numcols_global = other._numcols_global;
    _numrows_local = other._numrows_local;
    _numcols_local = other._numcols_local;
    _firstrownum_globalcount = other._firstrownum_globalcount;
    _data = other._data;
    _colindex = other._colindex;
    _firstrowentry = other._firstrowentry;
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._firstrownum_globalcount = 0;
}//CsrMatrix& CsrMatrix::operator=(CsrMatrix&& other)

bool CsrMatrix::operator==(const CsrMatrix& other)
{
    if (*this == other) return true;
    if (_numrows_global != other._numcols_global
        || _numcols_global != other._numcols_global)
        return false;
    for (size_t i{0}; i < _data.size(); ++i)
        if (_data[i] != other._data[i])
            return false;
    return true;
}//bool CsrMatrix::operator==(const CsrMatrix& other)

bool CsrMatrix::operator!=(const CsrMatrix& other)
{
    return !operator==(other);
}//bool CsrMatrix::operator!=(const CsrMatrix& other)

void CsrMatrix::free_unused_space()
{
    _data.shrink_to_fit();
    _colindex.shrink_to_fit();
}

double CsrMatrix::get_global(const size_t row, const size_t col) const
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
}//double CsrMatrix::get_global(const size_t row, const size_t col) const

double CsrMatrix::get_local(const size_t row, const size_t col) const
{
    assert(row < _numrows_local && col < _numcols_local);
    for (size_t datalocation{_firstrowentry[row]}; datalocation < _firstrowentry[row+1]; ++datalocation)
    {
        if (_colindex[datalocation] == col)
            return _data[datalocation];
        else if (_colindex[datalocation] > col)
            return 0.0;
    }
    return 0.0;
}//double CsrMatrix::get_local(const size_t row, const size_t col) const

void CsrMatrix::set_global(const size_t row, const size_t col, const double val)
{
    assert(row < _numrows_global && col < _numcols_global);
    if (row >= _firstrownum_globalcount && row < _firstrownum_globalcount + _numrows_local)
        set_local(row - _firstrownum_globalcount, col, val);
    MPICALL(MPI::COMM_WORLD.Barrier();) //TODISCUSS necessary?
}//void CsrMatrix::set_global(const size_t row, const size_t col, const double val)

void CsrMatrix::set_local(const size_t row, const size_t col, const double val)
{
    assert(row < _numrows_local && col < _numcols_local);
    size_t pos_to_insert{_firstrowentry[row]};
    bool new_entry{true};
    while (pos_to_insert < _firstrowentry[row+1])
    {
        if (_colindex[pos_to_insert] < col)
            ++pos_to_insert;
        else if (_colindex[pos_to_insert] == col)
        {
            _data[pos_to_insert] = val;
            return;
        }
    }
    _colindex.insert(_colindex.begin() + pos_to_insert, col);
    _data.insert(_data.begin() + pos_to_insert, val);
    for (size_t i{row + 1}; i <= _numrows_local; ++i)
        ++_firstrowentry[i]; //TOCHECK does this work?
        //_firstrowentry[i] += 1; // otherwise do this
}//void CsrMatrix::set_local(const size_t row, const size_t col, const double val)

double CsrMatrix::norm_1() const
{
    //TODO
    return 0.0;
}//double CsrMatrix::norm_1() const

double CsrMatrix::norm_2() const
{
    double local_squaresum{0.0};
    for (const auto val : _data)
        local_squaresum += val*val;
    double global_squaresum{0.0};
    MPICALL(MPI::COMM_WORLD.Allreduce(&local_squaresum, &global_squaresum, 1, MPI_DOUBLE, MPI_SUM);)
    return std::sqrt(global_squaresum);
}//double CsrMatrix::norm_2() const

double CsrMatrix::norm_inf() const
{
    //TODO
    return 0.0;
}//double CsrMatrix::norm_inf() const

CsrMatrix& CsrMatrix::get_scal_mul(const double scal) const
{
    CsrMatrix* res = new CsrMatrix(*this);
    if (scal == 1.0) return *res;
    for (auto& val : res->_data)
        val *= scal;
    return *res;
}//CsrMatrix& CsrMatrix::get_scal_mul(const double scal) const

void CsrMatrix::scal_mul(const double scal)
{
    if (scal != 1.0)
        for (auto& val : _data)
            val *= scal;
}//void CsrMatrix::scal_mul(const double scal)

DenseVector& CsrMatrix::get_vec_mul(const DenseVector& vec) const
{
    assert(vec.get_size_global() == _numcols_global);
    std::cout << _numrows_global << std::endl;
    std::cout << _numcols_global << std::endl;
    DenseVector* res = new DenseVector(_numrows_global);
    for (size_t i{0}; i < _numrows_local; ++i)
    {
        double tmp{0.0};
        for (size_t j{_firstrowentry[i]}; j < _firstrowentry[j+1]; ++j)
        {
            tmp += _data[j] * vec.get_global(_colindex[j]);
            std::cout << _colindex[j] << std::endl;
        }
        res->add_local(i, tmp);
    }
    std::cout << "haha" << std::endl;
    return *res;
}//DenseVector& CsrMatrix::get_vec_mul(const DenseVector& vec) const

DenseVector& CsrMatrix::get_pre_vec_mul(const DenseVector& vec) const
{
    assert(vec.get_size_global() == _numrows_global);
    DenseVector* res = new DenseVector(_numrows_global);
    //TODO
    return *res;
}//DenseVector& CsrMatrix::get_pre_vec_mul(const DenseVector& vec) const

}//namespace hptypes
