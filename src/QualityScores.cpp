#include "QualityScores.hpp"
#include <sstream>

namespace SGLib{
namespace QS{

unsigned short int getSangerQRToValue(char qs)
{
    return static_cast<unsigned short int>(qs - 33);
}

std::vector<unsigned short int> getSangerQRToValues(const std::string& quality_read)
{
    std::vector<unsigned short int> toRet(quality_read.size());
    for(std::size_t i = 0; i < quality_read.size(); ++i) toRet[i] = static_cast<unsigned short int>(quality_read[i]) - 33;
    return toRet;
}

//runtime data dependent error
std::vector<unsigned short int> getIlluminaQRToValues(const std::string& quality_read)
{
    std::vector<unsigned short int> toRet(quality_read.size());
    for(std::size_t i = 0; i < quality_read.size(); ++i)
    {
        short int value = static_cast<short int>(quality_read[i]) - 64;
        if(value < 0) throw std::invalid_argument("This function does not handle the Illumina 1.0 alphabet version");
        toRet[i] = static_cast<unsigned short int>(value);
    }
    return toRet;
}

std::vector<short int> getIlluminaQR1_0ToValues(const std::string& quality_read)
{
    std::vector<short int> toRet(quality_read.size());
    for(std::size_t i = 0; i < quality_read.size(); ++i)
    {
        toRet[i] = static_cast<short int>(quality_read[i]) - 64;
    }
    return toRet;
}

std::vector<double> getSangerQRToProbabilities(const std::string& quality_read)
{
    auto vals = getSangerQRToValues(quality_read);
    std::vector<double> toRet(vals.size());
    for(std::size_t i = 0; i < vals.size(); ++i) toRet[i] = std::pow(10, - static_cast<float>(vals[i]) / 10);
    return toRet;
}

std::vector<double> getIlluminaQRToProbabilities(const std::string& quality_read, Quality_t illumina_alphabet)
{
    if(illumina_alphabet == IlluminaQR_1_8) return getSangerQRToProbabilities(quality_read);
    else
    {
        auto vals = getIlluminaQRToValues(quality_read);
        std::vector<double> toRet(vals.size());
        if(illumina_alphabet == IlluminaQR_1_3 || illumina_alphabet == IlluminaQR_1_5)
        {
            for(std::size_t i = 0; i < vals.size(); ++i) toRet[i] = std::pow(10, - static_cast<float>(vals[i]) / 10);
        }
        else if(illumina_alphabet == IlluminaQR_1_0)
        {
            for(std::size_t i = 0; i < vals.size(); ++i)
            {
                double c = std::pow(10, - static_cast<float>(vals[i]) / 10);
                toRet[i] = c / (1 + c);
            }
        }
        else
        {
            throw std::invalid_argument("The alphabet is not one of the Illumina family");
        }
        return toRet;
    }
}

char getSangerQRFromValue(unsigned short int qval)
{
    return static_cast<char>(qval + 33);
}

std::string getSangerQRFromValues(const std::vector<unsigned short int>& qvals)
{
    std::stringstream sstr;
    for(unsigned short int qval : qvals) sstr << static_cast<char>(qval + 33);
    return sstr.str();
}

std::string getIlluminaQRFromValues(const std::vector<unsigned short int>& qvals)
{
    std::stringstream sstr;
    for(unsigned short int qval : qvals) sstr << static_cast<char>(qval + 64);
    return sstr.str();
}

std::string getIlluminaQR1_0FromValues(const std::vector<short int>& qvals)
{
    std::stringstream sstr;
    for(short int qval : qvals) sstr << static_cast<char>(qval + 64);
    return sstr.str();
}

std::string getSangerQRFromProbabilities(const std::vector<float>& pvals)
{
    std::stringstream sstr;
    for(float pval : pvals) sstr << static_cast<char>(-10 * std::log10(pval));
    return sstr.str();
}

std::string getIlluminaQRFromProbabilities(const std::vector<float>& pvals, Quality_t illumina_alphabet)
{
    if(illumina_alphabet == IlluminaQR_1_0)
    {
        std::stringstream sstr;
        for(float pval : pvals) sstr << static_cast<char>(-10 * std::log10(pval/(1-pval)));
        return sstr.str();
    }
    else return getSangerQRFromProbabilities(pvals);
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(const std::string& read, const std::vector<unsigned short>& qvals, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max, std::size_t k_min)
{
    return getSubsequencesBetweenQRVS(read, 0, read.length(), qvals, 0, qvals.size(), qr_threshold_min, qr_threshold_max, k_min);
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(const std::string& read, std::size_t rstart, std::size_t rstop, const std::vector<unsigned short>& qvals, std::size_t qstart, std::size_t qstop, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max, std::size_t k_min)
{
    if(rstart > read.length()) rstart = read.length();
    if(rstop > read.length()) rstop = read.length();

    if(qstart > qvals.size()) rstart = qvals.size();
    if(qstop > qvals.size()) rstop = qvals.size();

    if(rstart > rstop)
    {
        std::size_t buf = rstart;
        rstart = rstop;
        rstop = buf;
    }
    if(qstart > qstop)
    {
        std::size_t buf = qstart;
        qstart = qstop;
        qstop = buf;
    }
    return getSubsequencesBetweenQRVS(read.cbegin() + static_cast<long>(rstart), read.cbegin() + static_cast<long>(rstop), qvals.cbegin() + static_cast<long>(qstart), qvals.cbegin() + static_cast<long>(qstop), qr_threshold_min, qr_threshold_max, k_min);
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(std::string::const_iterator rstart, std::string::const_iterator rstop, std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max, std::size_t k_min)
{
    static_cast<void>(qstop);
    auto rpos = rstart;
    auto qpos = qstart;
    std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> toRet;
    while(static_cast<std::size_t>(rstop - rstart) >= k_min && rpos != rstop)
    {
        if(*qpos >= qr_threshold_min && *qpos < qr_threshold_max)
        {
            ++rpos;
            ++qpos;
        }
        else
        {
            if(rpos - rstart >= static_cast<long>(k_min)) toRet.push_back(std::pair<std::string::const_iterator, std::string::const_iterator>(rstart, rpos));
            rstart = ++rpos;
            qstart = ++qpos;
        }
    }
    if(rpos - rstart >= static_cast<long>(k_min)) toRet.push_back(std::pair<std::string::const_iterator, std::string::const_iterator>(rstart, rpos)); //check the last subsequence aligned to the end.
    return toRet;
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(const std::string& read, const std::vector<unsigned short>& qvals, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max, std::size_t k_min )
{
    return getSubsequencesWithQRVS(read, 0, read.length(), qvals, 0, qvals.size(), qr_threshold_min, qr_threshold_max, k_min);
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(const std::string& read, std::size_t rstart, std::size_t rstop, const std::vector<unsigned short>& qvals, std::size_t qstart, std::size_t qstop, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max,std::size_t k_min)
{
    if(rstart > read.length()) rstart = read.length();
    if(rstop > read.length()) rstop = read.length();

    if(qstart > qvals.size()) rstart = qvals.size();
    if(qstop > qvals.size()) rstop = qvals.size();

    if(rstart > rstop)
    {
        std::size_t buf = rstart;
        rstart = rstop;
        rstop = buf;
    }
    if(qstart > qstop)
    {
        std::size_t buf = qstart;
        qstart = qstop;
        qstop = buf;
    }
    return getSubsequencesWithQRVS(read.cbegin() + static_cast<long>(rstart), read.cbegin() + static_cast<long>(rstop), qvals.cbegin() + static_cast<long>(qstart), qvals.cbegin() + static_cast<long>(qstop), qr_threshold_min, qr_threshold_max, k_min);
}

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(std::string::const_iterator rstart, std::string::const_iterator rstop, std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop, unsigned short int qr_threshold_min, unsigned short int qr_threshold_max, std::size_t k_min)
{
    static_cast<void>(qstop);
    auto rpos = rstart;
    auto qpos = qstart;
    bool first = false;
    bool second = false;
    std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> toRet;
    std::pair<std::string::const_iterator, std::string::const_iterator> toPush;
    while(rpos != rstop)
    {
        //iterate until a qval such that qr_threshold_min <= qval < qr_threshold_max is found
        if(*qpos < qr_threshold_min || *qpos >= qr_threshold_max)
        {
            ++rpos;
            ++qpos;
        }
        else
        {
            //if the starting position of the pair is not yet assigned assign it such that it points to the suitable position minus k_min,
            //that is the starting value will have the position aligned to the right. If not possible the starting position is the start iterator.
            if(!first)
            {
                first = true;
                if(static_cast<std::size_t>(rpos - rstart) >= k_min)
                    toPush.first = rpos - static_cast<long>(k_min) + 1;
                else toPush.first = rstart;
            }
            //iterate through the quality values starting from the previously found position for a length of k_min,
            //if no other suitable positions are found assign the qpos + k_min position as end of pair otherwise
            //do nothing (wait for the next loop to find a second for the pair).
            if(!second)
            {
                size_t counter = 0;
                while((rpos != rstop) && (counter < k_min || (*qpos >= qr_threshold_min && *qpos < qr_threshold_max)))
                {
                    ++rpos;
                    ++qpos;
                    ++counter;
                }
                if(*(qpos - 1) < qr_threshold_min || *(qpos - 1) >= qr_threshold_max || rpos == rstop)
                {
                    toPush.second = rpos;
                    second = true;
                }
                else
                {
                    --qpos;
                    --rpos;
                }
            }
        }

        if(first && second)
        {
            first = second = false;
            toRet.push_back(toPush);
        }
    }
    if(first && !second)
    {
        toPush.second = rstop;
        toRet.push_back(toPush);
    }
    return toRet;
}

//--------------------------------------------------------------------------------------------------------------

namespace Compression{

const std::array<unsigned short int, 94> IL8B = {{
    0, 1, 6, 6, 6, 6, 6, 6, 6, 6, 15, 15, 15, 15, 15, 15,
    15, 15, 15, 15, 22, 22, 22, 22, 22, 27, 27, 27, 27, 27, 33, 33,
    33, 33, 33, 37, 37, 37, 37, 37, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40
}};

char illumina8Bin(char qs)
{
    return getSangerQRFromValue(IL8B[getSangerQRToValue(qs)]);
}

std::vector<unsigned short int> illumina8Bin(const std::vector<unsigned short>& qvals)
{
    return illumina8Bin(qvals.cbegin(), qvals.cend());
}

std::vector<unsigned short int> illumina8Bin(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop)
{
    std::vector<unsigned short int> toRet(static_cast<std::size_t>(qstop - qstart));
    for(std::size_t i = 0; qstart != qstop; ++qstart) toRet[i++] = IL8B[*qstart];
    return toRet;
}

std::vector<unsigned short int> rblock(const std::vector<unsigned short>& qvals, double theta)
{
    return rblock(qvals.begin(), qvals.end(), theta);
}

std::vector<unsigned short int> rblock(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop,  double theta)
{
    std::vector<unsigned short int> toRet(static_cast<std::size_t>(qstop - qstart));
    if(qstart != qstop)
    {
        unsigned short int qr_min = *qstart;
        unsigned short int qr_max = *qstart;
        std::size_t pos = 0;
        std::size_t old_pos = 0;
        while(qstart != qstop)
        {
            if(*qstart <= qr_max && *qstart >= qr_min)
            {
                ++qstart;
                ++pos;
                continue;
            }
            if(*qstart > qr_max)
            {
                unsigned short representative = static_cast<unsigned short>(round(sqrt(qr_min * (*qstart))));
                if(static_cast<double>(representative) / qr_min < theta && static_cast<double>(*qstart) / representative < theta)
                {
                    qr_max = *qstart++;
                    ++pos;
                    continue;
                }
            }
            else if(*qstart < qr_min)
            {
                unsigned short representative = static_cast<unsigned short>(round(sqrt(qr_max * (*qstart))));
                if(static_cast<double>(representative) / *qstart < theta && static_cast<double>(qr_max) / representative < theta)
                {
                    qr_min = *qstart++;
                    ++pos;
                    continue;
                }
            }
            unsigned short representative = static_cast<unsigned short>(round(sqrt(qr_min * qr_max)));
            for(; old_pos < pos; ++old_pos) toRet[old_pos] = representative;
            qr_min = *qstart++;
            qr_max = qr_min;
            ++pos;
        }
        unsigned short representative = static_cast<unsigned short>(round(sqrt(qr_min * qr_max)));
        if(pos != toRet.size()) throw std::length_error("pos is not equal to the size of the returning vector");
        for(; old_pos < pos; ++old_pos) toRet[old_pos] = representative;
    }
    return toRet;
}

std::vector<unsigned short int> pblock(const std::vector<unsigned short>& qvals, unsigned short two_p)
{
    return pblock(qvals.cbegin(), qvals.cend(), two_p);
}

std::vector<unsigned short int> pblock(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop,  unsigned short two_p)
{
    unsigned short qr_min = *qstart;
    unsigned short qr_max = qr_min;
    std::size_t pos = 0;
    std::size_t old_pos = 0;

    std::vector<unsigned short int> toRet(static_cast<std::size_t>(qstop - qstart));
    while(qstart != qstop)
    {
        if(*qstart <= qr_max && *qstart >= qr_min)
        {
            ++qstart;
            ++pos;
            continue;
        }
        if(*qstart > qr_max && *qstart - qr_min <= two_p)
        {
            qr_max = *qstart++;
            ++pos;
            continue;
        }
        if(*qstart < qr_min && qr_max - *qstart <= two_p)
        {
            qr_min = *qstart++;
            ++pos;
            continue;
        }
        unsigned short int representative = (qr_max + qr_min) / 2;
        for(; old_pos < pos; ++old_pos) toRet[old_pos] = representative;
        qr_min = *qstart++;
        qr_max = qr_min;
        ++pos;
    }
    unsigned short int representative = (qr_max + qr_min) / 2;
    for(; old_pos < pos; ++old_pos) toRet[old_pos] = representative;
    return toRet;
}

}//Compression

}//QS
}//SGLib
