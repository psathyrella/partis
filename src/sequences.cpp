#include "sequences.h"

namespace stochhmm {
Sequences::~Sequences() {
  for(size_t i=0;i<seqs.size();i++){
    delete seqs[i];
    seqs[i] = NULL;
  }
}
  
//!Get the digitized value from the sequence at trck at position
//! \param trck Sequence track to use
//! \param position Position in sequence to get the value from
//! \return short digitized value of the sequence based on track type 
short Sequences::getValue(int trck, size_t position){
  return seqs[trck]->getValue(position);
}
  
//!Get pointer to ith sequence from sequences
//! \param iter Iterator to use for extracting sequence;
//! \return sequence* pointer to sequence
//! \return NULL if no sequence exists at iter
Sequence* Sequences::getSeq(size_t iter) {
  if(iter<seqs.size()){
    return seqs[iter];
  }
  return NULL;
}

// ----------------------------------------------------------------------------------------
string& Sequences::getHeader(size_t iter) {
  assert(iter < seqs.size());
  return seqs[iter]->header_;
}

// ----------------------------------------------------------------------------------------
string Sequences::stringifyWOHeader() {
  string tmp_str;
  for(size_t i=0; i<seqs.size(); i++)
    tmp_str += seqs[i]->stringifyWOHeader();
  return tmp_str;
}

// ----------------------------------------------------------------------------------------
string Sequences::stringify() {
  string tmp_str;
  for(size_t i=0; i<seqs.size(); i++) {
    track* trk = seqs[i]->getTrack();
    tmp_str += ">" + trk->getName() + "\n";
    tmp_str += seqs[i]->stringifyWOHeader();
  }
  return tmp_str;
}

// ----------------------------------------------------------------------------------------
string Sequences::undigitize() {
  string output;
  for(size_t i=0; i<seqs.size(); i++) {
    track* trk = seqs[i]->getTrack();
    output += ">" + trk->getName();
    output += seqs[i]->undigitize();
  }
  return output;
}
  
// ----------------------------------------------------------------------------------------
//! \return return the collection of subsequences from <pos> with size <len>
Sequences Sequences::getSubSequences(size_t pos, size_t len) {
  Sequences new_seqs;  // init with *zero* seqs because we push back below
  for (size_t is=0; is<seqs.size(); is++) {
    Sequence *tmp_seq = new(nothrow) Sequence(seqs[is]->getSubSequence(pos,len));
    assert(tmp_seq);
    new_seqs.addSeq(tmp_seq);
  }
  return new_seqs;
}    

// ----------------------------------------------------------------------------------------
void Sequences::addSeq(Sequence* sq) {
  if (sequence_length_==0) {  // make sure <sq> has the same length as any previously added sequences
    sequence_length_ = sq->size();
  } else {
    assert(sq->size() == sequence_length_);
  }
  seqs.push_back(sq);  // NOTE we now own this sequence, i.e. we will delete it when we die
}
}
