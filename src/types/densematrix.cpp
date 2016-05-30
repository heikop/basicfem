#include "include/densematrix.hpp"

namespace hptypes
{

DenseMatrix::DenseMatrix(const DenseMatrix& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _rowptrs{new double*[_numrows_local]}, _data{new double[_numrows_local*_numcols_local]}
{
    for (size_t i{0}; i < _numrows_local; ++i)
        _rowptrs[i] = _data + i*_numcols_local;
    for (size_t i{0}; i < _numrows_local*_numcols_local; ++i)
        _data[i] = other._data[i];
    //TODO error handling allocating
}//DenseMatrix::DenseMatrix(const DenseMatrix& other)

DenseMatrix::DenseMatrix(DenseMatrix&& other):
    _numrows_global{other._numrows_global}, _numcols_global{other._numcols_global},
    _numrows_local{other._numrows_local}, _numcols_local{other._numcols_local},
    _rowptrs{other._rowptrs}, _data{other._data}
{
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    other._rowptrs = nullptr;
    other._data = nullptr;
}//DenseMatrix::DenseMatrix(DenseMatrix&& other)

DenseMatrix::DenseMatrix(size_t numrows, size_t numcols):
    _numrows_global{numrows}, _numcols_global{numcols},
    _numrows_local{numrows}, _numcols_local{numcols}, // TODO
    _rowptrs{new double*[_numrows_local]}, _data{new double[_numrows_local*_numcols_local]}
{
    for (size_t i{0}; i < _numrows_local; ++i)
        _rowptrs[i] = _data + i*_numcols_local;
    for (size_t i{0}; i < _numrows_local*_numcols_local; ++i)
        _data[i] = 0.0;
    //TODO error handling allocating
}//DenseMatrix::DenseMatrix(size_t numrows, size_t numcols)

DenseMatrix& DenseMatrix::operator= (const DenseMatrix& other)
{
    if (this == &other) return *this;
    if (_data) delete[] _data;
    if (_rowptrs) delete[] _rowptrs;
    _numrows_global = other.get_numrows_global();
    _numcols_global = other.get_numcols_global();
    _numrows_local = other.get_numrows_local();
    _numcols_local = other.get_numcols_local();
    _rowptrs = new double*[_numrows_local];
    _data = new double[_numrows_local*_numcols_local];
    for (size_t i{0}; i < _numrows_local; ++i)
        _rowptrs[i] = _data + i*_numcols_local;
    for (size_t i{0}; i < _numrows_local*_numcols_local; ++i)
        _data[i] = other._data[i];
    return *this;
    //TODO error handling allocating
}//DenseMatrix& DenseMatrix::operator= (const DenseMatrix& other)

DenseMatrix& DenseMatrix::operator= (DenseMatrix&& other)
{
    if (this == &other) return *this;
    if (_data) delete[] _data;
    if (_rowptrs) delete[] _rowptrs;
    _numrows_global = other.get_numrows_global();
    _numcols_global = other.get_numcols_global();
    _numrows_local = other.get_numrows_local();
    _numcols_local = other.get_numcols_local();
    _rowptrs = other._rowptrs;
    _data = other._data;
    other._rowptrs = nullptr;
    other._data = nullptr;
    other._numrows_global = 0;
    other._numcols_global = 0;
    other._numrows_local = 0;
    other._numcols_local = 0;
    return *this;
    //TODO error handling allocating
}//DenseMatrix& DenseMatrix::operator= (DenseMatrix&& other)

void DenseMatrix::print_local() const
{
    for (size_t i{0}; i < _numrows_local; ++i)
    {
        for (size_t j{0}; j < _numcols_local; ++j)
            std::cout << _rowptrs[i][j] << " ";
        std::cout << std::endl;
    }
}//void DenseMatrix::print_local() const

bool DenseMatrix::issymmetric() const
{
    if (!isquadratic()) return false;
    for (size_t i{0}; i < _numrows_local; ++i)
        for (size_t j{i+1}; j < _numcols_local; j++)
            if (_data[i*_numcols_local + j] != _data[j*_numcols_local + i])
                return false;
    return true;
}//bool DenseMatrix::issymmetric() const

double DenseMatrix::norm_1() const
{
    double res{0};
    for (size_t col{0}; col < _numcols_local; ++col)
    {
        double tmp{0};
        for (double** rowptr{_rowptrs}, ** upto{_rowptrs + _numrows_local}; rowptr < upto; ++rowptr)
            tmp += *(*rowptr + col);
        if (tmp > res) res = tmp;
    }
    return res;
}//double DenseMatrix::norm_1() const

double DenseMatrix::norm_2() const
{
    double res{0};
    for (double* data{_data}, * upto{_data + _numrows_local*_numcols_local}; data < upto; ++data)
        res += (*data) * (*data);
    return res;
}//double DenseMatrix::norm_2() const

double DenseMatrix::norm_inf() const
{
    double res{0};
    for (double** rowptr{_rowptrs}, ** upto{_rowptrs + _numrows_local}; rowptr < upto; ++rowptr)
    {
        double tmp{0};
        for (double* dataptr{*rowptr}, * upto{*rowptr + _numcols_local}; dataptr < upto; ++dataptr)
            tmp += *dataptr;
        if (tmp > res) res = tmp;
    }
    return res;
}//double DenseMatrix::norm_inf() const

DenseMatrix& DenseMatrix::get_transpose() const
{
    assert(_data);
    DenseMatrix* transposed = new DenseMatrix(_numcols_global, _numrows_global);
    size_t colnum{0};
    for (double* dataptr{_data}, * upto{_data + _numrows_local*_numcols_local}; dataptr < upto;)
    {
        for (double** rowptr{transposed->_rowptrs}, **upto{transposed->_rowptrs + transposed->_numrows_local}; rowptr < upto; ++rowptr, ++dataptr)
            *(*rowptr + colnum) = *dataptr;
        ++colnum;
    }
    return *transposed;
}//DenseMatrix& DenseMatrix::get_transpose() const

void DenseMatrix::transpose()
{
//    if (!_data) return;
//    DenseMatrix transposed(_numcols_global, _numrows_global);
//    size_t colnum{0};
//    for (double* dataptr{_data}, *upto{_data + _numrows_local*_numcols_local}; dataptr < upto;)
//    {
//        for (double** rowptr{transposed._rowptrs}, **upto{transposed._rowptrs + transposed._numrows_local}; rowptr < upto; ++rowptr, ++dataptr)
//            *(*rowptr + colnum) = *dataptr;
//        ++colnum;
//    }
//    *this = std::move(transposed);

    //*this = std::move(get_transpose()); //TODO TOTHINK what is better?
    *this = get_transpose();
}//void DenseMatrix::transpose()

DenseMatrix& DenseMatrix::get_mat_add(const DenseMatrix& other) const
{
    DenseMatrix* res = new DenseMatrix(*this);
    assert(_numrows_global == other._numrows_global && _numcols_global == other._numcols_global);
    for (double* data{res->_data}, * otherdata{other._data}, * upto{res->_data + _numrows_local*_numcols_local};
                                                                        data < upto; ++data, ++otherdata)
        *data += *otherdata;
    return *res;
}//DenseMatrix& DenseMatrix::get_mat_add(const DenseMatrix& other) const

void DenseMatrix::mat_add(const DenseMatrix& other)
{
    *this = get_mat_add(other);
}//void DenseMatrix::mat_add(const DenseMatrix& other)

DenseMatrix& DenseMatrix::get_mat_sub(const DenseMatrix& other) const
{
    DenseMatrix* res = new DenseMatrix(*this);
    assert(_numrows_global == other._numrows_global && _numcols_global == other._numcols_global);
    for (double* data{res->_data}, * otherdata{other._data}, * upto{res->_data + _numrows_local*_numcols_local};
                                                                        data < upto; ++data, ++otherdata)
        *data -= *otherdata;
    return *res;
}//DenseMatrix& DenseMatrix::get_mat_sub(const DenseMatrix& other) const

void DenseMatrix::mat_sub(const DenseMatrix& other)
{
    *this = get_mat_sub(other);
}//void DenseMatrix::mat_sub(const DenseMatrix& other)

DenseMatrix& DenseMatrix::get_scal_mul(const double scalar) const
{
    DenseMatrix* res = new DenseMatrix(*this);
    res->scal_mul(scalar);
    return *res;
}//DenseMatrix& DenseMatrix::get_scal_mul(const double scalar) const

void DenseMatrix::scal_mul(const double scalar)
{
    for (double* data{_data}, * upto{_data + _numrows_local*_numcols_local}; data < upto; ++data)
        *data *= scalar;
}//void DenseMatrix::scal_mul(const double scalar)

//void DenseMatrix::pow(const unsigned int exp)
//{
//}void DenseMatrix::pow(const unsigned int exp)

DenseMatrix& DenseMatrix::get_mat_mul(const DenseMatrix& other) const
{
    assert(_numcols_global == other._numrows_global);
    DenseMatrix* res = new DenseMatrix(_numrows_local, _numrows_local);
    double* resdata{res->_data};
    for (size_t i{0}; i < res->_numrows_local; ++i)
    {
        for (size_t j{0}; j < res->_numcols_local; ++j)
        {
            double* rowval{*(_rowptrs + i)};
            double* colval{other._data + j};
            double tmp{0};
            for (size_t l{0}; l < _numcols_local; ++l)
            {
                tmp += (*rowval) * (*colval);
                ++rowval;
                colval += other._numcols_local;
            }
            *resdata = tmp;
            ++resdata;
        }
    }
    return *res;
}//DenseMatrix& DenseMatrix::get_mat_mul(const DenseMatrix& other) const

void DenseMatrix::mat_mul(const DenseMatrix& other)
{
//    if (&other == this) { pow(2); return; }
    *this = get_mat_mul(other);
}//void DenseMatrix::mat_mul(const DenseMatrix& other)

DenseMatrix& DenseMatrix::get_inverse() const
{
    assert(isquadratic());
    DenseMatrix* inverse = new DenseMatrix(get_transpose());
    inverse->scal_mul(1 / this->norm_inf() / inverse->norm_inf()); // = inf * 1 norm. inf norm is easier to compute -> faster
    DenseMatrix two_idmat(_numrows_local, _numcols_local);
    DenseMatrix idmat(_numrows_local, _numcols_local);
    for (size_t i{0}; i < two_idmat._numrows_local; ++i)
    {
        idmat.set_local(i, i, 1.0);
        two_idmat.set_local(i, i, 2.0);
    }
    for (size_t i{0}; std::abs((this->get_mat_mul(*inverse)).get_mat_sub(idmat).norm_2()) > 1.0e-12 && i < 10000; ++i)
    {
        inverse->mat_mul(two_idmat.get_mat_sub(this->get_mat_mul(*inverse)));
        std::cout << "after " << i << "iterations: " << std::abs((this->get_mat_mul(*inverse)).get_mat_sub(idmat).norm_2()) << std::endl;
    }
    return *inverse;
}//DenseMatrix& DenseMatrix::get_inverse() const

void DenseMatrix::invert()
{
    //*this = std::move(get_inverse()); //TODO TOTHINK what is better?
    *this = get_inverse();
}//void DenseMatrix::invert()

}//namespace hptypes
