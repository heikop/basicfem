#include "../src/types/include/coomatrix.hpp"

#include "../src/helplib/include/helplib.hpp"
#include "include/sparsematrixtests.hpp"

int main()
{
    hptypes::CooMatrix* mata;
    hptypes::CooMatrix* matb;
    hptypes::CooMatrix* matc;
    sparsetest_basics(mata, matb, matc);

    return 0;
}
