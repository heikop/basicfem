#include "../src/types/include/lilmatrix.hpp"
#include "../src/types/include/densevector.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/basicmatrixtests.hpp"
#include "include/vecmultest.hpp"

int main()
{
    hptypes::LilMatrix* mata;
    hptypes::LilMatrix* matb;
    hptypes::LilMatrix* matc;
    basicmatrixtest(mata, matb, matc);

    hptypes::LilMatrix* matd;
    hptypes::DenseVector* vecd;
    vecmultest(matd, vecd);

    return 0;
}
