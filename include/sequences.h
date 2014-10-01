#ifndef STOCHHHMM_SEQUENCES_H
#define STOCHHHMM_SEQUENCES_H

#include "sequence.h"

using namespace std;
namespace stochhmm {

class Sequences{
public:
  Sequences():sequence_length_(0){}
  ~Sequences();
      
  short getValue(int , size_t); // Get digitized value for sequence track(i) in jth position
  Sequence* getSeq(size_t); // Return ith sequence
  Sequence& operator[](size_t index) {return *seqs[index]; }
  string& getHeader(size_t iter=0);
  string* getUndigitized(size_t iter) { return seqs[iter]->getUndigitized(); }
  size_t n_seqs() { return seqs.size(); }
  size_t GetSequenceLength() { return sequence_length_;}
  Sequences getSubSequences(size_t pos, size_t len);
  
  string stringifyWOHeader();  //! Get string of sequences
  string stringify();  //! Get string of sequences
  void print() { cout << stringify() << endl; }
  string undigitize(); //! Get sequence based on alphabet
  void addSeq(Sequence* sq);
private:
  vector<Sequence*> seqs; 
  size_t sequence_length_;  // length of the sequences (required to be the same for all)
};
}
#endif
