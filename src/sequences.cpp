#include "sequences.h"

namespace StochHMM {
//!Destroy sequences
sequences::~sequences() {
  for(size_t i=0;i<seqs.size();i++){
    delete seqs[i];
    seqs[i] = NULL;
  }
}
  
//!Get the digitized value from the sequence at trck at position
//! \param trck Sequence track to use
//! \param position Position in sequence to get the value from
//! \return short digitized value of the sequence based on track type 
short sequences::seqValue(int trck, size_t position){
  return seqs[trck]->seqValue(position);
}
  
//!Get pointer to ith sequence from sequences
//! \param iter Iterator to use for extracting sequence;
//! \return sequence* pointer to sequence
//! \return NULL if no sequence exists at iter
sequence* sequences::getSeq(size_t iter) {
  if(iter<seqs.size()){
    return seqs[iter];
  }
  return NULL;
}

// ----------------------------------------------------------------------------------------
std::string& sequences::getHeader(size_t iter) {
  assert(iter < seqs.size());
  return seqs[iter]->header;
}

// ----------------------------------------------------------------------------------------
std::string sequences::stringifyWOHeader() {
  std::string tmp_str;
  for(size_t i=0; i<size(); i++)
    tmp_str += seqs[i]->stringifyWOHeader();
  return tmp_str;
}

// ----------------------------------------------------------------------------------------
std::string sequences::stringify() {
  std::string tmp_str;
  for(size_t i=0; i<size(); i++) {
    track* trk = seqs[i]->getTrack();
    tmp_str += ">" + trk->getName() + "\n";
    tmp_str += seqs[i]->stringifyWOHeader();
  }
  return tmp_str;
}

// ----------------------------------------------------------------------------------------
std::string sequences::undigitize() {
  std::string output;
  for(size_t i=0; i<size(); i++) {
    track* trk = seqs[i]->getTrack();
    output += ">" + trk->getName();
    output += seqs[i]->undigitize();
  }
  return output;
}
  
// ----------------------------------------------------------------------------------------
//! \return return the collection of subsequences from <pos> with size <len>
sequences sequences::getSubSequences(size_t pos, size_t len) {
  sequences new_seqs;  // init with *zero* seqs because we push back below
  for (size_t is=0; is<size(); is++) {
    sequence *tmp_seq = new(std::nothrow) sequence(seqs[is]->getSubSequence(pos,len));
    assert(tmp_seq);
    new_seqs.addSeq(tmp_seq);
  }
  return new_seqs;
}    

// ----------------------------------------------------------------------------------------
void sequences::addSeq(sequence* sq) {
  if (length==0) {  // make sure <sq> has the same length as any previously added sequences
    length = sq->size();
  } else {
    assert(sq->size() == length);
  }
  seqs.push_back(sq);  // NOTE we now own this sequence, i.e. we will delete it when we die
}
}
