#include "../src/types/include/csrmatrix.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/sparsematrixtests.hpp"

int main()
{
    hptypes::CsrMatrix* mata;
    hptypes::CsrMatrix* matb;
    hptypes::CsrMatrix* matc;
    sparsetest_basics(mata, matb, matc);

    return 0;
}
