#include "Fastx.hpp"
#include <memory>
#include <limits>

namespace SGLib{
namespace Fastx{

const FXA_t UnambiguousDNA = "ACGTacgt";

FastxRecord::FastxRecord(std::string id, std::string sequence) noexcept : Tracing::BaseRecord(id)
{
    genome = std::move(sequence);
}

const std::string& FastxRecord::getSequenceRef() const
{
    return genome;
}

FastxRecord::~FastxRecord()
{
    //dtor
}

//--------------------------------------------------------------------------------------------------------------

std::size_t countRecords(std::istream& is, delim_t record_delim)
{
    std::size_t result = 0;

    std::string line;
    while(std::getline(is, line))
    {
        if(line.length() != 0 && line[0] == record_delim) ++result;
        //is.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    return result;
}

}//Fastx
}//SGLib
