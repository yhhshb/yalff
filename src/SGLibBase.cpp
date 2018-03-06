#include "SGLibBase.hpp"
#include <algorithm>

namespace SGLib{

/**
 * Find the next available subsequence defined over the given alphabet.
 *
 * The value returned MUST be checked if the search is unsuccessful (because there are no substrings of length k strictly defined over the alphabet).
 *
 * @param seq The entire sequence with alphabet a strictly superset of the given one.
 * @param start The start position of the search.
 * @param alphabet The alphabet over which the subsequences must be defined to be considered correct.
 * @param k The length of the desired substring.
 * @return The starting position of the desired substring. seq.length() if no suitable subsequence found.
 */
std::size_t fetchNextSubsequence(const std::string& seq, std::size_t start, std::size_t stop, const IUPAC_t& alphabet, std::size_t k)
{
    if(start > stop)
    {
        auto buffer = start;
        start = stop;
        stop = buffer;
    }
    if(start > seq.length() || stop > seq.length()) return seq.length();
    return static_cast<std::size_t>(fetchNextSubsequence(seq.begin() + static_cast<long>(start), seq.begin() + static_cast<long>(stop), alphabet, k) - seq.begin());
}

/**
 * Find the next available subsequence defined over the given alphabet.
 *
 * The value returned MUST be checked if the search is unsuccessful (because there are no substrings of length k strictly defined over the alphabet).
 *
 * @param start The start iterator.
 * @param stop The stop iterator.
 * @param alphabet The alphabet over which the subsequences must be defined to be considered correct.
 * @param k The length of the desired substring.
 * @return A constant iterator pointing to the starting position of the substring. stop if no suitable subsequence found.
 */
std::string::const_iterator fetchNextSubsequence(std::string::const_iterator start, std::string::const_iterator stop, const IUPAC_t& alphabet, std::size_t k)
{
    auto in = [](char c, const IUPAC_t& alpha)
    {
        return std::find(alpha.begin(), alpha.end(), c) != alpha.end();
    };
    auto pos = start;
    std::size_t counter = 0;
    while(counter != k && static_cast<std::size_t>(stop - start) >= k)
    {
        if(in(*pos++, alphabet)) ++counter;
        else
        {
            start = pos;
            counter = 0;
        }
    }
    if(counter == k) return start;
    else return stop;
}

std::pair<std::size_t, std::size_t> getPositionsFromIterators(const std::string& read, const std::pair<std::string::const_iterator, std::string::const_iterator>& itrs)
{
    return getPositionsFromIterators(read.cbegin(), itrs);
}

std::pair<std::size_t, std::size_t> getPositionsFromIterators(std::string::const_iterator begin_, const std::pair<std::string::const_iterator, std::string::const_iterator>& itrs)
{
    return std::pair<std::size_t, std::size_t>(itrs.first - begin_, itrs.second - begin_);
}

std::vector<std::pair<std::size_t, std::size_t>> getPositionsFromIterators(const std::string& read, const std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>>& itrs)
{
    return getPositionsFromIterators(read.cbegin(), itrs);
}

std::vector<std::pair<std::size_t, std::size_t>> getPositionsFromIterators(std::string::const_iterator begin_, const std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>>& itrs)
{
    std::vector<std::pair<std::size_t, std::size_t>> toRet(itrs.size());
    for(std::size_t i = 0; i < itrs.size(); ++i)
    {
        toRet[i] = std::pair<std::size_t, std::size_t>(itrs[i].first - begin_, itrs[i].second - begin_);
    }
    return toRet;
}

namespace Predicates{

std::function<bool(const std::string&, const IUPAC_t&)> isDefinedOverAlphabet = [](const std::string& str , const IUPAC_t& alphabet)
{
    for(char c : str) if(std::find(alphabet.cbegin(), alphabet.cend(), c) == alphabet.cend()) return false;
    return true;
};

std::function<bool(const std::string&)> alwaysTrue = [](const std::string&)
{
    return true;
};

std::function<bool(const std::string&)> isUnambigousDNA = std::bind(isDefinedOverAlphabet, std::placeholders::_1, std::ref(IUPACUnambiguousDNA));

}//Predicates

}//SGLib
