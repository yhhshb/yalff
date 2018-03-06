#ifndef SGLIBEXCEPTIONS_HPP
#define SGLIBEXCEPTIONS_HPP

#include <exception>

namespace SGLib{

class invalid_fastx : public std::exception
{
    public:
        virtual const char* what() const noexcept;
};

class fastx_end : public std::exception
{
    public:
        virtual const char* what() const noexcept;
};

}//SGLib

#endif // SGLIBEXCEPTIONS_HPP
