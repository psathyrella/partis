#ifndef HAM_SEQUENCE_H
#define HAM_SEQUENCE_H

#include "track.h"

using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
class Sequence {
  friend class Sequences;
public:
  Sequence(string name, string seq_str, Track* trk);
  Sequence(const Sequence&);
  ~Sequence();

  inline string name() const { return name_; }
  inline uint8_t operator[](size_t index) { return seq_->at(index); }  // digitized value at position <index>
  inline uint8_t value(size_t pos) const { return (*seq_)[pos]; }  // get digitized value at <pos>
  inline string symbol(size_t pos) const { return track_->symbol((*seq_)[pos]); }  // get undigitized value at <pos>
  inline size_t size() const { return seq_->size(); }  // length of sequence
  inline Track* track() const { return track_; }
  inline string undigitized() { return undigitized_; }
  Sequence GetSubSequence(size_t pos, size_t len);

  void Print(string separator=" ");  // if separator is specified, print it between each element in the sequence
private:
  string name_;
  string header_;
  string undigitized_;  // undigitized sequence
  Track* track_; // track describing alphabet and type. NOTE we don't own this pointer, i.e. we don't delete it when we die
  vector<uint8_t>* seq_; // digitized Sequence
};

// ----------------------------------------------------------------------------------------
// NOTE Sequences owns its constituent Sequence pointers, i.e. it deletes them when it dies
class Sequences {
public:
  Sequences() : sequence_length_(0) {}
  void AddSeq(Sequence* sq);
  ~Sequences();

  short value(size_t iseq, size_t ipos) { return seqs_.at(iseq)->value(ipos); }  // return digitized value of <iseq>th sequence at position <ipos>
  Sequence& operator[](size_t index) {return *(seqs_.at(index)); }
  size_t n_seqs() { return seqs_.size(); }
  size_t GetSequenceLength() { return sequence_length_;}
  Sequences GetSubSequences(size_t pos, size_t len);

  void Print();
private:
  vector<Sequence*> seqs_;
  size_t sequence_length_; // length of the sequences (required to be the same for all)
};

}
#endif
