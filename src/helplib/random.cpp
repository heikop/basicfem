#include "include/random.hpp"

namespace hphelp
{

void activate_random() { srand(time(NULL)); }

size_t randomsize(size_t min, size_t max)
{
    assert(min > 0 && max > 0);
    //srand(time(NULL));
    if (min <= max)
        return rand() % (max - min + 1) + min;
    return rand() % (min - max + 1) + max;
}

std::vector<size_t>& randompositions(size_t numdata, size_t maxpos)
{
    std::vector<size_t>* randomposvector = new std::vector<size_t>(numdata);
    for (auto x : (*randomposvector))
        x = rand() % maxpos;
    return *randomposvector;
}

std::vector<double>& randomvalues(size_t numdata, int maxnumdigits, int exponent)
{
    std::vector<double>* randomdatavector = new std::vector<double>(numdata);
    int maxnum(1);
    while (--maxnumdigits >= 0) maxnum *= 10;
    for (size_t i{0}; i < numdata; ++i)
        (*randomdatavector)[i] = static_cast<double>(rand() % maxnum) * static_cast<double>(std::pow(10, exponent)); // TODO make this better
    return *randomdatavector;
}


}//namespace hphelp
