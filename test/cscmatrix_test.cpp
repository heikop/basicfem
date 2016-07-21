#include "../src/types/include/cscmatrix.hpp"
#include "../src/types/include/densevector.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/basicmatrixtests.hpp"
#include "include/vecmultest.hpp"

int main()
{
    hptypes::CscMatrix* mata;
    hptypes::CscMatrix* matb;
    hptypes::CscMatrix* matc;
    basicmatrixtest(mata, matb, matc);

    hptypes::CscMatrix* matd;
    hptypes::DenseVector* vecd;
    vecmultest(matd, vecd);

    return 0;
}
