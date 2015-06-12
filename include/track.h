#ifndef HAM_TRACK_H
#define HAM_TRACK_H

#include <map>
#include <cassert>
#include <stdexcept>

#include "text.h"
#include "mathutils.h"

using namespace std;
namespace ham {

// ----------------------------------------------------------------------------------------
class Track {
public:
  Track() {}
  Track(string name, vector<string> symbols, string ambiguous_char = "");
  void set_name(string nm) { name_ = nm; }
  void AddSymbol(string symbol);
  void AddSymbols(vector<string> &symbols);
  void SetAmbiguous(string amb) { ambiguous_char_ = amb; }
  string ambiguous_char() { return ambiguous_char_; }
  uint8_t ambiguous_index() { return ambiguous_index_; }

  string name() { return name_; }
  size_t alphabet_size() { return alphabet_.size(); }
  string symbol(size_t iter) { return alphabet_.at(iter); }  // return <iter>th element of <alphabet_> (which is probably a letter, but could be several letters or something else). NOTE throws std:out_of_range exception if <iter> is invalide
  uint8_t symbol_index(const string &symbol);
  string Stringify();
private:
  string name_;
  vector<string> alphabet_;  // vector of this track's allowed symbols (eg {A,C,G,T})
  map<string, uint8_t> symbol_indices_;
  string ambiguous_char_;
  const static uint8_t max_alphabet_size_ = 255;
  const static uint8_t ambiguous_index_ = max_alphabet_size_ - 1;
};

}
#endif
