#include "../src/types/include/densematrix.hpp"
#include "../src/types/include/densevector.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/basicmatrixtests.hpp"
#include "include/vecmultest.hpp"

int main()
{
    hptypes::DenseMatrix* mata;
    hptypes::DenseMatrix* matb;
    hptypes::DenseMatrix* matc;
    basicmatrixtest(mata, matb, matc);

    hptypes::DenseMatrix* matd;
    hptypes::DenseVector* vecd;
    vecmultest(matd, vecd);

    return 0;
}
