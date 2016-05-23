#include "../src/types/include/cscmatrix.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/sparsematrixtests.hpp"

int main()
{
    hptypes::CscMatrix* mata;
    hptypes::CscMatrix* matb;
    hptypes::CscMatrix* matc;
    sparsetest_basics(mata, matb, matc);

    return 0;
}
