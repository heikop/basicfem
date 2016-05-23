#ifndef __SPARSEMATRIXTESTS_HPP_
#define __SPARSEMATRIXTESTS_HPP_

#include <stdexcept>
#include <vector>

#include "../../src/helplib/include/helplib.hpp"

template <typename mattype>
void sparsetest_basics(mattype* mata, mattype* matb, mattype* matc, size_t num_tests=50)
{
    hphelp::activate_random();
    for (size_t num_checks{0}; num_checks < num_tests; ++num_checks)
    {
        // 0-initialization
        const size_t numrows{hphelp::randomsize(10e1, 10e2)};
        const size_t numcols{hphelp::randomsize(10e1, 10e2)};
        //hptypes::CooMatrix mata(numrows, numcols);
        mata = new mattype(numrows, numcols);
        if (mata->get_numrows_global() != numrows)
            throw std::length_error("initialized with wrong row size");
        if (mata->get_numcols_global() != numcols)
            throw std::length_error("initialized with wrong col size");
        //TODO check sum over all local sizes
        for (size_t i{0}; i < mata->get_numrows_local(); ++i)
            for (size_t j{0}; j < mata->get_numcols_local(); ++j)
                if (mata->get_local(i, j) != 0.0)
                    throw std::invalid_argument("not zero-initialized/ get_local not correct");

        // set values - global
        const size_t numglobvals{hphelp::randomsize(10 * (numrows + numcols), numrows * numcols / 10)};
        const std::vector<size_t> rowposglob{hphelp::randompositions(numglobvals, numcols)};
        const std::vector<size_t> colposglob{hphelp::randompositions(numglobvals, numrows)};
        const std::vector<double> globalvals{hphelp::randomvalues(numglobvals, 8, -4)};
        for (size_t i{0}; i < numglobvals; ++i)
        {
            mata->set_global(rowposglob[i], colposglob[i], globalvals[i]);
            if (mata->get_global(rowposglob[i], colposglob[i]) != globalvals[i])
                throw std::invalid_argument("set_global/ get_global not correct");
        }
        // set values - local
        const size_t numlocvals{hphelp::randomsize(10 * (numrows + numcols), numrows * numcols / 10)};
        const std::vector<size_t> rowposloc{hphelp::randompositions(numlocvals, numcols)};
        const std::vector<size_t> colposloc{hphelp::randompositions(numlocvals, numrows)};
        const std::vector<double> localvals{hphelp::randomvalues(numlocvals, 8, -4)};
        for (size_t i{0}; i < numlocvals; ++i)
        {
            mata->set_local(rowposloc[i], colposloc[i], localvals[i]);
            if (mata->get_local(rowposloc[i], colposloc[i]) != localvals[i])
                throw std::invalid_argument("set_local/ get_local not correct");
        }

        // copy-constructor
        //hptypes::CooMatrix matb(mata);
        matb = new mattype(*mata);
        if (matb->get_numrows_global() != mata->get_numrows_global())
            throw std::length_error("copy-constructor: wrong global row size");
        if (matb->get_numcols_global() != mata->get_numcols_global())
            throw std::length_error("copy-constructor: wrong global col size");
        if (matb->get_numrows_local() != mata->get_numrows_local())
            throw std::length_error("copy-constructor: wrong local row size");
        if (matb->get_numcols_local() != mata->get_numcols_local())
            throw std::length_error("copy-constructor: wrong local col size");
        for (size_t i{0}; i < matb->get_numrows_local(); ++i)
            for (size_t j{0}; j < matb->get_numcols_local(); ++j)
                if (matb->get_local(i, j) != mata->get_local(i, j))
                    throw std::invalid_argument("copy-constructor: value wrong assigned");

        // assignment-operator
        //hptypes::CooMatrix matc(0, 0);
        matc = new mattype(0, 0);
        *matc = *mata;
        if (matc->get_numrows_global() != mata->get_numrows_global())
            throw std::length_error("assignment-operator: wrong global row size");
        if (matc->get_numcols_global() != mata->get_numcols_global())
            throw std::length_error("assignment-operator: wrong global col size");
        if (matc->get_numrows_local() != mata->get_numrows_local())
            throw std::length_error("assignment-operator: wrong local row size");
        if (matc->get_numcols_local() != mata->get_numcols_local())
            throw std::length_error("assignment-operator: wrong local col size");
        for (size_t i{0}; i < matc->get_numrows_local(); ++i)
            for (size_t j{0}; j < matc->get_numcols_local(); ++j)
                if (matc->get_local(i, j) != mata->get_local(i, j))
                    throw std::invalid_argument("assignment-operator: value wrong assigned");

        // 1 norm
        //TODO

        // 2 norm
        //TODO

        // inf norm
        //TODO
        delete mata;
        delete matb;
        delete matc;
    }
}

#endif//ifndef __SPARSEMATRIXTESTS_HPP_
