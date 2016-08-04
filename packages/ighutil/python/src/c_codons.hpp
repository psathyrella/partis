#ifndef CODONS_H
#define CODONS_H
#include <string>

namespace codons
{

unsigned short pack_codon(const char n1, const char n2, const char n3);

char translate_codon(const char n1, const char n2, const char n3);

/// simple 1st frame forward translation of a given DNA string
std::string translate_dna(const std::string& dnastr);

/// Convert a string to upper case
void toupper(std::string& s);

bool codon_table_init();

} // namespace codonalign

#endif
