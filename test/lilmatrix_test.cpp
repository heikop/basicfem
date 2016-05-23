#include "../src/types/include/lilmatrix.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/sparsematrixtests.hpp"

int main()
{
    hptypes::LilMatrix* mata;
    hptypes::LilMatrix* matb;
    hptypes::LilMatrix* matc;
    sparsetest_basics(mata, matb, matc);

    return 0;
}
