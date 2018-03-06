#ifndef SGLIBBASE_HPP
#define SGLIBBASE_HPP

#include <vector>
#include "SGLibConstants.hpp"

namespace SGLib{

std::size_t fetchNextSubsequence(const std::string& seq, std::size_t start, std::size_t stop, const IUPAC_t& alphabet, std::size_t k);
std::string::const_iterator fetchNextSubsequence(std::string::const_iterator start, std::string::const_iterator stop, const IUPAC_t& alphabet, std::size_t k);

std::pair<std::size_t, std::size_t> getPositionsFromIterators(const std::string& read, const std::pair<std::string::const_iterator, std::string::const_iterator>& itrs);
std::pair<std::size_t, std::size_t> getPositionsFromIterators(std::string::const_iterator begin_, const std::pair<std::string::const_iterator, std::string::const_iterator>& itrs);
std::vector<std::pair<std::size_t, std::size_t>> getPositionsFromIterators(const std::string& read, const std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>>& itrs);
std::vector<std::pair<std::size_t, std::size_t>> getPositionsFromIterators(std::string::const_iterator begin_, const std::vector<std::pair<std::string::const_iterator, std::string::const_iterator>>& itrs);

template<class UIntType>
inline UIntType rotl(UIntType v, std::size_t shift){
    static_assert(std::is_unsigned<UIntType>::value, "Rotation only makes sense for unsigned types");
    std::size_t m = sizeof(v) * 8; //std::numeric_limits<UIntType>::digits;
    //std::cout << "divisore: " << m << std::endl;
    UIntType s = shift % m;
    //std::cout << "shift: " << s << std::endl;
    return (v<<s) | (v>>(m-s));
}

template<class UIntType>
inline UIntType rotr(UIntType v, std::size_t shift){
    static_assert(std::is_unsigned<UIntType>::value, "Rotation only makes sense for unsigned types");
    std::size_t m = sizeof(v) * 8;
    UIntType s = shift % m;
    return (v>>s) | (v<<(m-s));
}

namespace Predicates{

extern std::function<bool(const std::string&, const IUPAC_t&)> isDefinedOverAlphabet;
extern std::function<bool(const std::string&)> alwaysTrue;
extern std::function<bool(const std::string&)> isUnambigousDNA;

}//Predicates

}//SGLib

#endif // SGLIBBASE_HPP
