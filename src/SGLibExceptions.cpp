#include "SGLibExceptions.hpp"

namespace SGLib{

const char* invalid_fastx::what() const noexcept
{
    return "The given input stream does not follow the supposed convention!";
}

const char* fastx_end::what() const noexcept
{
    return "The fast* file is finished, no more records in it";
}

}//SGLib
