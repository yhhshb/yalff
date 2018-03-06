#include "Fastq.hpp"
#include <sstream>

namespace SGLib{
namespace Fastq{

Fastx::delim_t FASTQ_START_DELIM = '@';
Fastx::delim_t FASTQ_MID_DELIM = '+';

FastqRecord::FastqRecord() noexcept : Fastx::FastxRecord("", "")
{
    qreads = "";
}

FastqRecord::FastqRecord(std::string header, std::string sequence, std::string quality) : Fastx::FastxRecord(header, std::move(sequence))
{
    if(genome.length() != quality.length())
    {
        throw std::runtime_error("The length of the quality string does not match the length of the sequence!");
    }
    qreads = std::move(quality);
}

const std::string& FastqRecord::getQualityRef() const
{
    return qreads;
}

char FastqRecord::getQualityScore(std::string::size_type pos) const
{
    return qreads.at(pos);
}

void FastqRecord::setQualityScore(std::string::size_type pos, char qs)
{
    qreads.at(pos) = qs;
}

FastqRecord::~FastqRecord()
{
    //dtor
}

//--------------------------------------------------------------------------------------------------------------



struct __ProtoFastqRecord{
    std::string iddesc;
    std::string sequence;
    std::string quality_reads;
};

__ProtoFastqRecord __readNextSimpleRecord(std::istream& is); //suppress warning
__ProtoFastqRecord __readNextSimpleRecord(std::istream& is)
{
    using namespace std;
    __ProtoFastqRecord toRet;
    std::string delim;

    is >> std::ws;

    if(std::getline(is, toRet.iddesc));
    else throw fastx_end();

    if(toRet.iddesc.length() != 0 && toRet.iddesc[0] == FASTQ_START_DELIM) toRet.iddesc = toRet.iddesc.erase(0, 1);
    else throw invalid_fastx();

    if(std::getline(is, toRet.sequence));
    else throw invalid_fastx();

    if(std::getline(is, delim));
    else throw invalid_fastx();

    if(delim.length() == 0 || delim[0] != FASTQ_MID_DELIM) throw invalid_fastx();

    if(std::getline(is, toRet.quality_reads));
    else throw invalid_fastx();

    //if(toRet.sequence.length() != toRet.quality_reads.length()) throw Incompatible_Quality_Reads_Exception();

    return toRet;
}

/**
 * Simple parser for fastq files containing multiple records.
 * All the records MUST be in the 4 line format.
 *
 * Read the next available record from the input stream.
 *
 * @param is The input stream.
 * @return The FastqRecord object.
 */
FastqRecord readNextFastqRecord(std::istream& is)
{
    __ProtoFastqRecord prec = __readNextSimpleRecord(is);

    std::istringstream iss(prec.iddesc);
    std::string id, description;
    iss >> id;
    iss >> description;
    FastqRecord toRet = FastqRecord(id, std::move(prec.sequence), std::move(prec.quality_reads));
    toRet.setDescription(description);
    return toRet;
}

/**
 * Filter a fastq stream.
 * Given two predicates, one for the sequence and the other for the quality the record goes in output iff both are true.
 *
 * @param is The input stream.
 * @param os The output stream.
 * @param genomic_predicate The predicate for the sequence.
 * @param quality_predicate The predicate for the quality sequence.
 */
void filterFastq(std::istream& is, std::ostream& os, std::function<bool(const std::string&)> genomic_predicate, std::function<bool(const std::string&)> quality_predicate)
{
    bool first = true;
    while(true)
    {
        std::unique_ptr<SGLib::Fastq::FastqRecord> fq_record;
        try{
            fq_record = std::make_unique<SGLib::Fastq::FastqRecord>(SGLib::Fastq::readNextFastqRecord(is));
        } catch(SGLib::fastx_end& file_ended) {
            static_cast<void>(file_ended);
            break;
        }

        if(genomic_predicate(fq_record->getSequenceRef()) && quality_predicate(fq_record->getQualityRef()))
        {
            if(first)
            {
                first = false;
                os << *fq_record;
            }
            else os << "\n" << *fq_record;
        }
    }
}

/**
 * Filter two fastq streams.
 * Given two predicates, one for the sequence and the other for the quality reads the record goes in output iff both are true.
 * This version takes two inputs and two outputs for paired end fastq files.
 * The two streams MUST be aligned and the records MUST be in the same relative order.
 *
 * @param is1 The first input stream.
 * @param is2 The second input stream.
 * @param os1 The first output stream.
 * @param os2 The second output stream.
 * @param genomic_predicate The predicate for the sequence.
 * @param quality_predicate The predicate for the quality sequence.
 */
void filterPairedEndFastq(std::istream& is1, std::istream& is2, std::ostream& os1, std::ostream& os2, std::function<bool(const std::string&)> genomic_predicate, std::function<bool(const std::string&)> quality_predicate)
{
    bool first = true;
    while(true)
    {
        std::unique_ptr<SGLib::Fastq::FastqRecord> fq_record1;
        std::unique_ptr<SGLib::Fastq::FastqRecord> fq_record2;
        try{
            fq_record1 = std::make_unique<SGLib::Fastq::FastqRecord>(SGLib::Fastq::readNextFastqRecord(is1));
            fq_record2 = std::make_unique<SGLib::Fastq::FastqRecord>(SGLib::Fastq::readNextFastqRecord(is2));
        } catch(SGLib::fastx_end& file_ended) {
            static_cast<void>(file_ended);
            break;
        }

        if(genomic_predicate(fq_record1->getSequenceRef()) && quality_predicate(fq_record1->getQualityRef()) &&
           genomic_predicate(fq_record2->getSequenceRef()) && quality_predicate(fq_record2->getQualityRef()))
        {
            if(first)
            {
                first = false;
                os1 << *fq_record1;
                os2 << *fq_record2;
            }
            else
            {
                os1 << "\n" << *fq_record1;
                os2 << "\n" << *fq_record2;
            }
        }
    }
}

}//Fastq
}//SGLib
