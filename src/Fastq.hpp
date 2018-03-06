#ifndef FASTQLOADER_HPP
#define FASTQLOADER_HPP

#include "SGLibBase.hpp"
#include "SGLibExceptions.hpp"
#include "Fastx.hpp"

namespace SGLib{
namespace Fastq{

extern Fastx::delim_t FASTQ_START_DELIM;
extern Fastx::delim_t FASTQ_MID_DELIM;

class FastqRecord : public Fastx::FastxRecord
{
    public:

        FastqRecord(const FastqRecord&) = default;

        FastqRecord(FastqRecord&&) = default;

        FastqRecord() noexcept;

        FastqRecord(std::string header, std::string sequence, std::string quality);

        const std::string& getQualityRef() const;

        char getQualityScore(std::string::size_type pos) const;

        void setQualityScore(std::string::size_type pos, char qs);

        virtual ~FastqRecord();

    protected:
        std::string qreads;

    private:
        friend std::ostream& operator<<(std::ostream& os,const FastqRecord& obj)
        {
            os << "@" << obj.id << " " << obj.description << "\n";
            os << obj.genome << "\n";
            os << "+" << "\n";
            os << obj.qreads;
            return os;
        }
};

FastqRecord readNextFastqRecord(std::istream& is);
void filterFastq(std::istream& is, std::ostream& os, std::function<bool(const std::string&)> genomic_predicate, std::function<bool(const std::string&)> quality_predicate);
void filterPairedEndFastq(std::istream& is1, std::istream& is2, std::ostream& os1, std::ostream& os2, std::function<bool(const std::string&)> genomic_predicate, std::function<bool(const std::string&)> quality_predicate);

}//Fastq
}//SGLib

#endif // FASTQLOADER_HPP
