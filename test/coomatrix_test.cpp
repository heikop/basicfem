#include "../src/types/include/coomatrix.hpp"
#include "../src/types/include/densevector.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/basicmatrixtests.hpp"
#include "include/vecmultest.hpp"

int main()
{
    hptypes::CooMatrix* mata;
    hptypes::CooMatrix* matb;
    hptypes::CooMatrix* matc;
    basicmatrixtest(mata, matb, matc);

    hptypes::CooMatrix* matd;
    hptypes::DenseVector* vecd;
    vecmultest(matd, vecd);

    return 0;
}
