#ifndef StochHMM_sequences_h
#define StochHMM_sequences_h

#include "sequence.h"

namespace StochHMM {

class sequences{
public:
  sequences():sequence_length_(0){}
  ~sequences();
      
  short getValue(int , size_t); // Get digitized value for sequence track(i) in jth position
  sequence* getSeq(size_t); // Return ith sequence
  sequence& operator[](size_t index) {return *seqs[index]; }
  std::string& getHeader(size_t iter=0);
  std::string* getUndigitized(size_t iter) { return seqs[iter]->getUndigitized(); }
  size_t n_seqs() { return seqs.size(); }
  size_t GetSequenceLength() { return sequence_length_;}
  sequences getSubSequences(size_t pos, size_t len);
  
  std::string stringifyWOHeader();  //! Get string of sequences
  std::string stringify();  //! Get string of sequences
  void print() { std::cout << stringify() << std::endl; }
  std::string undigitize(); //! Get sequence based on alphabet
  void addSeq(sequence* sq);
private:
  std::vector<sequence*> seqs; 
  size_t sequence_length_;  // length of the sequences (required to be the same for all)
};
}
#endif
