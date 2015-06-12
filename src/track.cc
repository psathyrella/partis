#include "track.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Track::Track(string name, vector<string> symbols, string ambiguous_char) :
  name_(name),
  ambiguous_char_(ambiguous_char)
{
  AddSymbols(symbols);
}

// ----------------------------------------------------------------------------------------
void Track::AddSymbol(string symbol) {
  assert(alphabet_.size() < max_alphabet_size_);  // cannot (at the moment) have more than 255 symbols in an alphabet
  alphabet_.push_back(symbol);
  symbol_indices_[symbol] = alphabet_.size() - 1;
}

// ----------------------------------------------------------------------------------------
void Track::AddSymbols(vector<string> &symbols) {
  for(size_t i = 0; i < symbols.size(); ++i)
    AddSymbol(symbols[i]);
}

// ----------------------------------------------------------------------------------------
uint8_t Track::symbol_index(const string &symbol) {
  if(ambiguous_char_ != "" && symbol == ambiguous_char_)
    return ambiguous_index_;
  if(symbol_indices_.count(symbol) == 0)
    throw runtime_error("ERROR symbol '" + symbol + "' not found among " + Stringify());
  return symbol_indices_[symbol];
}

// ----------------------------------------------------------------------------------------
string Track::Stringify() {
  string return_str;
  for(auto & symbol : alphabet_)
    return_str += " " + symbol;
  if(ambiguous_char_ != "")
    return_str += "   ambiguous: " + ambiguous_char_;
  return return_str;
}

}
