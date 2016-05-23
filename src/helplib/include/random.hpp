#ifndef __RANDOM_HPP_
#define __RANDOM_HPP_

#include <vector>
#include <stdio.h> // NULL
#include <stdlib.h> // srand, rand, abs
#include <time.h> // time
#include <cmath> // std::pow, std::abs
#include <cassert>

namespace hphelp
{

// must be called in apps that use randomhelp functions
void activate_random();

size_t randomsize(size_t min, size_t max);

std::vector<size_t>& randompositions(size_t numdata, size_t maxpos);

std::vector<double>& randomvalues(size_t numdata, int maxnumdigits, int exponent);


}//namespace hphelp


#endif//__RANDOM_HPP_
