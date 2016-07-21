#ifndef __VECMULTEST_HPP_
#define __VECMULTEST_HPP_

#include <stdexcept>
#include <vector>
#include <mpi.h>

#include "../../src/helplib/include/helplib.hpp"
#include "../../src/mpihandler.hpp"

template <typename mattype, typename vectype>
void vecmultest(mattype* mat, vectype* vec, size_t num_tests=1)
//void vecmultest(mattype* mat, vectype* vec, size_t num_tests=50)
{
    hphelp::activate_random();
    for (size_t num_checks{0}; num_checks < num_tests; ++num_checks)
    {
        // matrix initialization
        const size_t numrows{hphelp::randomsize(1e1, 1e2)};
        const size_t numcols{hphelp::randomsize(1e1, 1e2)};
        const std::vector<double> matvals{hphelp::randomvalues(numrows*numcols, 8, -4)};
        mat = new mattype(numrows, numcols);
        for (size_t i{0}; i < numrows; ++i)
            for (size_t j{0}; j < numcols; ++j)
                mat->set_global(i, j, matvals[i*numcols + j]);

        // matrix * vector
        vec = new vectype(numcols);
        const std::vector<double> vecvals{hphelp::randomvalues(std::max(numrows, numcols), 4, -2)};

        //hptypes::DenseVector resa{mat->get_vec_mul(*vec)};
        vectype resa{mat->get_vec_mul(*vec)};

        delete mat;
    }
}//void vecmultest(mattype* mat, vectype* vec, size_t num_tests=50)


#endif//#ifndef __VECMULTEST_HPP_
