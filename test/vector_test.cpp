#include "../src/types/include/vector.hpp"

#include <stdexcept>
#include <vector>
#include <iostream>

#include "../src/helplib/include/helplib.hpp"

int main()
{
    hphelp::activate_random();
    for (size_t num_checks{0}; num_checks < 50; ++num_checks)
    {
        // 0-initialization
        size_t vecsize{hphelp::randomsize(10e2, 10e3)};
        hptypes::Vector veca(vecsize);
        if (veca.get_size_global() != vecsize)
            throw std::length_error("initialized with wrong size");
        //TODO check sum over all local lengths
        for (size_t i{0}; i < veca.get_size_local(); ++i)
            if (veca.get_local(i) != 0.0)
                throw std::invalid_argument("not zero-initialized/ get_local not correct");

        // set values
        std::vector<double> vecvalues(hphelp::randomvalues(vecsize, 10, -5));
        for (size_t i{0}; i < vecsize; ++i)
        {
            veca.set_local(i, vecvalues[i]);
            if (veca.get_local(i) != vecvalues[i])
                throw std::invalid_argument("set_local/ get_local not correct");
        }
        //TODO also with global set and get

        // copy-constructor
        hptypes::Vector vecb(veca);
        if (vecb.get_size_global() != veca.get_size_global())
            throw std::length_error("copy-constructor: wrong global size");
        if (vecb.get_size_local() != veca.get_size_local())
            throw std::length_error("copy-constructor: wrong local size");
        for (size_t i{0}; i < vecb.get_size_local(); ++i)
            if (vecb.get_local(i) != veca.get_local(i))
                throw std::invalid_argument("copy-constructor: value wrong assigned");

        // assignment-operator
        hptypes::Vector vecc(0);
        vecc = veca;
        if (vecc.get_size_global() != veca.get_size_global())
            throw std::length_error("assignment-operator: wrong global size");
        if (vecc.get_size_local() != veca.get_size_local())
            throw std::length_error("assignment-operator: wrong local size");
        for (size_t i{0}; i < vecc.get_size_local(); ++i)
            if (vecc.get_local(i) != veca.get_local(i))
                throw std::invalid_argument("assignment-operator: value wrong assigned");

        // binary + operator QUESTION: this test both move assignment and binary+ operators
        //vecc = veca + vecb;
        //for (size_t i{0}; i < vecc.get_size_local(); ++i)
        //    if (vecc.get_local(i) != veca.get_local(i) + vecb.get_local(i))
        //        throw std::invalid_argument("binary + operator: value wrong computed");
        //TODO unary-, +, +=, -, -= operators

        // norms
        if (veca.l2norm() != veca.lpnorm(2))
            throw std::invalid_argument("l2 and lp norm not equal");
        double veca_l1norm{0};
        //for (size_t i{0}; i < veca.get_size_local(); ++i)
        //    veca_l1norm += std::abs(veca.get_local(i));
        for (double x : vecvalues)
            veca_l1norm += std::abs(x);
        if (veca.lpnorm(1) != veca_l1norm)
            throw std::invalid_argument("l1norm is not computed right");
        double veca_maxnorm{0};
        for (double x : vecvalues)
            if (std::abs(x) > veca_maxnorm)
                veca_maxnorm = std::abs(x);
        if (veca.maxnorm() != veca_maxnorm)
            throw std::invalid_argument("maxnorm is not computed right");
    }

    return 0;
}
