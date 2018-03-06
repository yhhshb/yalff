#ifndef FASTXCONSTANTS_HPP
#define FASTXCONSTANTS_HPP

#include <iostream>
#include "BaseRecord.hpp"

namespace SGLib{
namespace Fastx{

typedef const char delim_t;
typedef std::string FXA_t;

extern const FXA_t UnambiguousDNA;

class FastxRecord : public Tracing::BaseRecord
{
    public:
        FastxRecord(const FastxRecord&) = default;

        FastxRecord(FastxRecord&&) = default;

        FastxRecord(std::string id, std::string sequence) noexcept;

        const std::string& getSequenceRef() const;

        virtual ~FastxRecord();

    protected:
        std::string genome;

    private:
        friend std::ostream& operator<<(std::ostream& os,const FastxRecord& obj)
        {
            os << obj.getId() << "\n";
            os << obj.genome;
            return os;
        }
};

std::size_t countRecords(std::istream& is, delim_t record_delim);

}//Fastx
}//SGLib

#endif // FASTXCONSTANTS_HPP
