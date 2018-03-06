#ifndef QUALITYSCORES_HPP
#define QUALITYSCORES_HPP

#include <cmath>
#include <algorithm>
#include "SGLibBase.hpp"

namespace SGLib{
namespace QS{

unsigned short int getSangerQRToValue(char qs);
std::vector<unsigned short int> getSangerQRToValues(const std::string& quality_read);
std::vector<unsigned short int> getIlluminaQRToValues(const std::string& quality_read);
std::vector<short int> getIlluminaQR1_0ToValues(const std::string& quality_read);
std::vector<double> getSangerQRToProbabilities(const std::string& quality_read);
std::vector<double> getIlluminaQRToProbabilities(const std::string& quality_read, Quality_t illumina_alphabet);

char getSangerQRFromValue(unsigned short int qval);
std::string getSangerQRFromValues(const std::vector<unsigned short int>& qvals);
std::string getIlluminaQRFromValues(const std::vector<unsigned short int>& qvals);
std::string getIlluminaQR1_0FromValues(const std::vector<short int>& qvals);
std::string getSangerQRFromProbabilities(const std::vector<float>& qvals);
std::string getIlluminaQRFromProbabilities(const std::vector<float>& quality_read, Quality_t illumina_alphabet);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(const std::string& read,
                                                                                                            const std::vector<unsigned short>& qvals,
                                                                                                            unsigned short int qr_threshold_min,
                                                                                                            unsigned short int qr_threshold_max,
                                                                                                            std::size_t k_min = 0);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(const std::string& read,
                                                                                                            std::size_t rstart,
                                                                                                            std::size_t rstop,
                                                                                                            const std::vector<unsigned short>& qvals,
                                                                                                            std::size_t qstart,
                                                                                                            std::size_t qstop,
                                                                                                            unsigned short int qr_threshold_min,
                                                                                                            unsigned short int qr_threshold_max,
                                                                                                            std::size_t k_min = 0);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesBetweenQRVS(std::string::const_iterator rstart,
                                                                                                            std::string::const_iterator rstop,
                                                                                                            std::vector<unsigned short>::const_iterator qstart,
                                                                                                            std::vector<unsigned short>::const_iterator qstop,
                                                                                                            unsigned short int qr_threshold_min,
                                                                                                            unsigned short int qr_threshold_max,
                                                                                                            std::size_t k_min = 0);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(const std::string& read,
                                                                                                         const std::vector<unsigned short>& qvals,
                                                                                                         unsigned short int qr_threshold_min,
                                                                                                         unsigned short int qr_threshold_max,
                                                                                                         std::size_t k_min);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(const std::string& read,
                                                                                                         std::size_t rstart,
                                                                                                         std::size_t rstop,
                                                                                                         const std::vector<unsigned short>& qvals,
                                                                                                         std::size_t qstart,
                                                                                                         std::size_t qstop,
                                                                                                         unsigned short int qr_threshold_min,
                                                                                                         unsigned short int qr_threshold_max,
                                                                                                         std::size_t k_min);

std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>> getSubsequencesWithQRVS(std::string::const_iterator rstart,
                                                                                                         std::string::const_iterator rstop,
                                                                                                         std::vector<unsigned short>::const_iterator qstart,
                                                                                                         std::vector<unsigned short>::const_iterator qstop,
                                                                                                         unsigned short int qr_threshold_min,
                                                                                                         unsigned short int qr_threshold_max,
                                                                                                         std::size_t k_min = 0);

template<typename CounterType, typename T>
std::vector<CounterType> countValues(const std::vector<T>& qvals);

template<typename CounterType, typename T>
std::vector<CounterType> countValues(typename std::vector<T>::const_iterator qstart, typename std::vector<T>::const_iterator qstop);

template<typename FloatType>
FloatType getProbability(const std::vector<double>& probabilities);

template<typename FloatType>
FloatType getProbability(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop);

template<typename FloatType>
FloatType getLogProbability(const std::vector<double>& probabilities, std::size_t base);

template<typename FloatType>
FloatType getLogProbability(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t base);

template<typename FloatType>
FloatType getProbability(FloatType prev_prob, FloatType out_prob, FloatType in_prob);

template<typename FloatType>
std::vector<FloatType> getProbabilities(const std::vector<double>& probabilities, std::size_t k, FloatType epsilon);

template<typename FloatType>
std::vector<FloatType> getProbabilities(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t k, FloatType epsilon);

template<typename FloatType>
std::vector<FloatType> getLogProbabilities(const std::vector<double>& probabilities, std::size_t k, std::size_t base);

template<typename FloatType>
std::vector<FloatType> getLogProbabilities(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t k, std::size_t base);

template<typename CounterType, typename FloatType>
CounterType countBadProbabilities(const std::vector<FloatType>& probabilities, FloatType pthr);

template<typename CounterType, typename FloatType>
CounterType countBadProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, FloatType pthr);

template<typename CounterType, typename FloatType>
std::vector<CounterType> countBadProbabilities(const std::vector<FloatType>& probabilities, std::size_t k, FloatType pthr);

template<typename CounterType, typename FloatType>
std::vector<CounterType> countBadProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, std::size_t k, FloatType pthr);

template<typename CounterType, typename FloatType>
CounterType countGoodProbabilities(const std::vector<FloatType>& probabilities, FloatType pthr);

template<typename CounterType, typename FloatType>
CounterType countGoodProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, FloatType pthr);

template<typename CounterType, typename FloatType>
std::vector<CounterType> countGoodProbabilities(const std::vector<FloatType>& probabilities, std::size_t k, FloatType pthr);

template<typename CounterType, typename FloatType>
std::vector<CounterType> countGoodProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, std::size_t k, FloatType pthr);

//--------------------------------------------------------------------------------------------------------------

namespace Compression{

extern const std::array<unsigned short int, 94> IL8B;

char illumina8Bin(char qs);

std::vector<unsigned short int> illumina8Bin(const std::vector<unsigned short>& qvals);

std::vector<unsigned short int> illumina8Bin(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop);

std::vector<unsigned short int> rblock(const std::vector<unsigned short>& qvals, double theta);

std::vector<unsigned short int> rblock(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop,  double theta);

std::vector<unsigned short int> pblock(const std::vector<unsigned short>& qvals, unsigned short two_p);

std::vector<unsigned short int> pblock(std::vector<unsigned short>::const_iterator qstart, std::vector<unsigned short>::const_iterator qstop,  unsigned short two_p);

}//Compression

//--------------------------------------------------------------------------------------------------------------

template<typename CounterType, typename T>
std::vector<CounterType> countValues(const std::vector<T>& qvals)
{
    return countValues<CounterType, T>(qvals.cbegin(), qvals.cend());
}

template<typename CounterType, typename T>
std::vector<CounterType> countValues(typename std::vector<T>::const_iterator qstart, typename std::vector<T>::const_iterator qstop)
{
    std::vector<CounterType> toRet;
    if(std::is_unsigned<T>::value) toRet.resize(94, 0);
    else toRet.resize(68, 0);
    for(;qstart != qstop; ++qstart)
    {
        if(*qstart >= 0) ++toRet[*qstart];
        else
        {
            unsigned short val = static_cast<unsigned short>(*qstart + 5);
            ++toRet[val];
        }
    }
    return toRet;
}

template<typename FloatType>
FloatType getProbability(const std::vector<double>& probabilities, FloatType epsilon)
{
    return getProbability<FloatType>(probabilities.cbegin(), probabilities.cend(), epsilon);
}

template<typename FloatType>
FloatType getProbability(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop)
{
    FloatType toRet = 1;
    for(; pstart != pstop; ++pstart) toRet *= 1.0 - *pstart;
    return toRet;
}

template<typename FloatType>
FloatType getLogProbability(const std::vector<double>& probabilities, std::size_t base)
{
    return getLogProbability<FloatType>(probabilities.cbegin(), probabilities.cend(), base);
}

template<typename FloatType>
FloatType getLogProbability(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t base)
{
    FloatType toRet = 1.0;
    for(; pstart != pstop; ++pstart) toRet += log(1.0 - *pstart) / log(base);
    return toRet;
}

template<typename FloatType>
FloatType getProbability(FloatType prev_prob, FloatType out_prob, FloatType in_prob)
{
    return prev_prob / (1.0 - out_prob) * (1.0 - in_prob);
}

template<typename FloatType>
std::vector<FloatType> getProbabilities(const std::vector<double>& probabilities, std::size_t k, FloatType epsilon)
{
    return getProbabilities<FloatType>(probabilities.cbegin(), probabilities.cend(), k, epsilon);
}

template<typename FloatType>
std::vector<FloatType> getProbabilities(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t k, FloatType epsilon)
{
    const long k1 = static_cast<long>(k);
    const long read_len = pstop - pstart;
    std::vector<FloatType> toRet;
    if(read_len < k1) return toRet;
    toRet.resize(static_cast<std::size_t>(read_len - k1 + 1));
    auto next_itr = pstart + k1;
    toRet[0] = getProbability<FloatType>(pstart, next_itr);
    for(std::size_t i = 1; i < static_cast<std::size_t>(read_len - k1 + 1); ++i)
    {
        if(toRet[i - 1] > epsilon)
        {
            toRet[i] = getProbability(toRet[i - 1], *pstart++, *next_itr++);
            //toRet[i] = toRet[i - 1] / (1.0 - *pstart++) * (1.0 - *next_itr++);
        }
        else
        {
            toRet[i] = getProbability<FloatType>(++pstart, ++next_itr);
        }
    }
    return toRet;
}

template<typename FloatType>
std::vector<FloatType> getLogProbabilities(const std::vector<double>& probabilities, std::size_t k, std::size_t base)
{
    return getLogProbabilities<FloatType>(probabilities.cbegin(), probabilities.cend(), k, base);
}

template<typename FloatType>
std::vector<FloatType> getLogProbabilities(std::vector<double>::const_iterator pstart, std::vector<double>::const_iterator pstop, std::size_t k, std::size_t base)
{
    const long k1 = static_cast<long>(k);
    const long read_len = pstop - pstart;
    std::vector<FloatType> toRet;
    if(read_len < k1) return toRet;
    toRet.resize(static_cast<std::size_t>(read_len - k1 + 1));
    auto next_itr = pstart + k1;
    toRet[0] = getLogProbability<FloatType>(pstart, next_itr, base);
    for(std::size_t i = 1; i < static_cast<std::size_t>(read_len - k1 + 1); ++i)
    {
        toRet[i] = toRet[i - 1] - log(1.0 - *pstart++) / log(base) + log(1.0 - *next_itr++) / log(base);
    }
    return toRet;
}

template<typename CounterType, typename FloatType> //Here the use of CounterType is for space purposes.
CounterType countBadProbabilities(const std::vector<FloatType>& probabilities, FloatType pthr)
{
    return countBadProbabilities<CounterType>(probabilities.cbegin(), probabilities.cend(), pthr);
}

template<typename CounterType, typename FloatType>
CounterType countBadProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, FloatType pthr)
{
    CounterType toRet = 0;
    while(pstart != pstop) if(*pstart++ > pthr) ++toRet;
    return toRet;
}

template<typename CounterType, typename FloatType>
std::vector<CounterType> countBadProbabilities(const std::vector<FloatType>& probabilities, std::size_t k, FloatType pthr)
{
    return countBadProbabilities<CounterType>(probabilities.cbegin(), probabilities.cend(), k, pthr);
}

template<typename CounterType, typename FloatType>
std::vector<CounterType> countBadProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, std::size_t k, FloatType pthr)
{
    if(pstop - pstart < 0)
    {
        auto buf = pstop;
        pstop = pstart;
        pstart = buf;
    }
    const CounterType read_len = static_cast<CounterType>(pstop - pstart);
    std::vector<CounterType> toRet;
    if(read_len < k) return toRet;
    toRet.resize(read_len - k + 1, 0);
    auto next_itr = pstart + static_cast<long>(k);
    toRet[0] = countBadProbabilities<CounterType>(pstart, next_itr, pthr);
    for(CounterType i = 1; i < static_cast<CounterType>(read_len - k + 1); ++i)
    {
        if(*pstart++ > pthr) toRet[i] = toRet[i-1] - 1;
        else toRet[i] = toRet[i-1];
        if(*next_itr++ > pthr) ++toRet[i];
    }
    return toRet;
}

template<typename CounterType, typename FloatType>
CounterType countGoodProbabilities(const std::vector<FloatType>& probabilities, FloatType pthr)
{
    return countGoodProbabilities<CounterType>(probabilities.cbegin(), probabilities.cend(), pthr);
}

template<typename CounterType, typename FloatType>
CounterType countGoodProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, FloatType pthr)
{
    CounterType toRet = 0;
    while(pstart != pstop) if(*pstart++ < pthr) ++toRet;
    return toRet;
}

template<typename CounterType, typename FloatType>
std::vector<CounterType> countGoodProbabilities(const std::vector<FloatType>& probabilities, std::size_t k, FloatType pthr)
{
    return countGoodProbabilities<CounterType>(probabilities.cbegin(), probabilities.cend(), k, pthr);
}

template<typename CounterType, typename FloatType>
std::vector<CounterType> countGoodProbabilities(typename std::vector<FloatType>::const_iterator pstart, typename std::vector<FloatType>::const_iterator pstop, std::size_t k, FloatType pthr)
{
    if(pstop - pstart < 0)
    {
        auto buf = pstop;
        pstop = pstart;
        pstart = buf;
    }
    const CounterType read_len = static_cast<CounterType>(pstop - pstart);
    std::vector<CounterType> toRet;
    if(read_len < k) return toRet;
    toRet.resize(read_len - k + 1, 0);
    auto next_itr = pstart + static_cast<long>(k);
    toRet[0] = countGoodProbabilities<CounterType>(pstart, next_itr, pthr);
    for(CounterType i = 1; i < static_cast<CounterType>(read_len - k + 1); ++i)
    {
        if(*pstart++ < pthr) toRet[i] = toRet[i-1] - 1;
        else toRet[i] = toRet[i-1];
        if(*next_itr++ < pthr) ++toRet[i];
    }
    return toRet;
}

}//QS
}//SGLib

#endif // QUALITYSCORES_HPP
