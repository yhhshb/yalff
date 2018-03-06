#ifndef SGLIBCONSTANTS_HPP
#define SGLIBCONSTANTS_HPP

#include <memory>
#include <string>

namespace SGLib{

typedef std::string IUPAC_t; //alphabet type
typedef std::string Quality_t; //quality reads type

extern const IUPAC_t ExtendedIUPACProtein; //Extended uppercase IUPAC protein single letter alphabet including X etc.
extern const IUPAC_t IUPACProtein; //Uppercase IUPAC protein single letter alphabet of the 20 standard amino acids.
extern const IUPAC_t IUPACAmbiguousDNA; //Uppercase IUPAC ambiguous DNA.
extern const IUPAC_t IUPACUnambiguousDNA; //Uppercase IUPAC unambiguous DNA (letters ACGT only).
extern const IUPAC_t ExtendedIUPACDNA; //Extended IUPAC DNA alphabet.
extern const IUPAC_t IUPACAmbiguousRNA; //Uppercase IUPAC ambiguous RNA.
extern const IUPAC_t IUPACUnambiguousRNA; //Uppercase IUPAC unambiguous RNA (letters ACGU only).

extern const std::shared_ptr<const IUPAC_t> ExtendedIUPACProteinPtr;
extern const std::shared_ptr<const IUPAC_t> IUPACProteinPtr;
extern const std::shared_ptr<const IUPAC_t> IUPACAmbiguousDNAPtr;
extern const std::shared_ptr<const IUPAC_t> IUPACUnambiguousDNAPtr;
extern const std::shared_ptr<const IUPAC_t> ExtendedIUPACDNAPtr;
extern const std::shared_ptr<const IUPAC_t> IUPACAmbiguousRNAPtr;
extern const std::shared_ptr<const IUPAC_t> IUPACUnambiguousRNAPtr;

extern const Quality_t SangerQR;
extern const Quality_t IlluminaQR_1_0;
extern const Quality_t IlluminaQR_1_3;
extern const Quality_t IlluminaQR_1_5;
extern const Quality_t IlluminaQR_1_8;

extern const std::shared_ptr<const Quality_t> SangerQRPtr;
extern const std::shared_ptr<const Quality_t> IlluminaQR_1_0Ptr;
extern const std::shared_ptr<const Quality_t> IlluminaQR_1_3Ptr;
extern const std::shared_ptr<const Quality_t> IlluminaQR_1_5Ptr;
extern const std::shared_ptr<const Quality_t> IlluminaQR_1_8Ptr;

extern const std::array<uint8_t, 256> base_to_int;

}//SGLib

#endif //SGLIBCONSTANTS_HPP
