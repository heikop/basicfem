#include "../src/types/include/csrmatrix.hpp"
#include "../src/types/include/densevector.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/basicmatrixtests.hpp"
#include "include/vecmultest.hpp"

int main()
{
    hptypes::CsrMatrix* mata;
    hptypes::CsrMatrix* matb;
    hptypes::CsrMatrix* matc;
    basicmatrixtest(mata, matb, matc);

    std::cout << "basic tests passed" << std::endl;
    hptypes::CsrMatrix* matd;
    hptypes::DenseVector* vecd;
    vecmultest(matd, vecd);

    return 0;
}
