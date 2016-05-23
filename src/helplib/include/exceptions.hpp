#ifndef __EXCEPTIONS_HPP_
#define __EXCEPTIONS_HPP_

#include <exception>
#include <string>

// many things TODO and even more to learn

namespace hphelp
{

class not_implemented : public std::exception
{
public:
    not_implemented() {};
    explicit not_implemented(const std::string& arg);
    explicit not_implemented(const char* arg);
    ~not_implemented() throw() {}
    const char* what() const throw() { return "NOT IMPLEMENTED"; }
};
//class not_implemented : public std::exception
//{
//public:
//    not_implemented() {}
//    ~not_implemented() throw () {}
//    const char* what() const throw () { return "NOT IMPLEMENTED"; }
//};

}//namespace hphelp

#endif//__EXCEPTIONS_HPP_
