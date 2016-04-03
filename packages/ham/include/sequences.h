#ifndef HAM_SEQUENCE_H
#define HAM_SEQUENCE_H

#include <set>
#include "track.h"

using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
class Sequence {
  friend class Sequences;
public:
  Sequence();  // NOTE don't use this! It's only so I can use stl maps without crashing
  Sequence(Track* trk, string name, string &undigitized);
  Sequence(Track* trk, string name, string &undigitized, size_t pos, size_t len);  // create the subsequence from <pos> of length <len>
  Sequence(Sequence &rhs, size_t pos, size_t len);
  Sequence(const Sequence &rhs);
  ~Sequence();

  inline string name() const { return name_; }
  inline void set_name(string name)  { name_ = name; }
  inline uint8_t operator[](size_t index) { return seqq_.at(index); }  // digitized value at position <index>
  inline uint8_t value(size_t pos) const { return seqq_[pos]; }  // get digitized value at <pos>
  inline string symbol(size_t pos) const { return track_->symbol(seqq_[pos]); }  // get undigitized value at <pos>
  inline size_t size() const { return seqq_.size(); }
  inline Track* track() const { return track_; }
  inline vector<uint8_t> *seqq() { return &seqq_; }
  inline string undigitized() const { return undigitized_; }
  Sequence GetSubSequence(size_t pos, size_t len);

  void Print(string separator = " "); // if separator is specified, print it between each element in the sequence
private:
  void Digitize();  // convert the string in <undigitized_> to a vector<uint8_t> in <seqq_>
  void CheckPosLen(string name, string undigitized, size_t pos, size_t len);
  string name_;
  string header_;
  string undigitized_;  // undigitized sequence
  Track* track_; // track describing alphabet and type. NOTE we don't own this pointer, i.e. we don't delete it when we die
  vector<uint8_t> seqq_; // digitized Sequence
};

// ----------------------------------------------------------------------------------------
class Sequences {
public:
  Sequences() : sequence_length_(0) {}
  // Sequences(const Sequences &rhs);
  Sequences(Sequences &rhs, size_t pos, size_t len);  // copy <seqs> from <pos> to <pos> + <len>
  // Sequences(vector<Sequence> &seqs);
  void AddSeq(Sequence sq);

  short value(size_t iseq, size_t ipos) { return seqs_.at(iseq).value(ipos); }  // return digitized value of <iseq>th sequence at position <ipos>
  Sequence &operator[](size_t index) { return seqs_.at(index); }
  // Sequence GetAtConst(size_t index) { return seqs_.at(index); }  // 
  Sequence *get_ptr(size_t index) { return &seqs_.at(index); }
  size_t n_seqs() const { return seqs_.size(); }
  size_t GetSequenceLength() { return sequence_length_;}
  Sequences Union(Sequences &otherseqs);  // return union set of self and <otherseqs>
  // Sequences GetSubSequences(size_t pos, size_t len);

  void Print();
  string name_str(string delimiter = " ") {
    string name_str;
    for(size_t iseq = 0; iseq < n_seqs(); ++iseq) {
      if(iseq > 0) name_str += delimiter;
      name_str += seqs_[iseq].name();
    }
    return name_str;
  }
  string seq_str(string delimiter = " ") {
    string seq_str;
    for(size_t iseq = 0; iseq < n_seqs(); ++iseq) {
      if(iseq > 0) seq_str += delimiter;
      seq_str += seqs_[iseq].undigitized_;
    }
    return seq_str;
  }
private:
  vector<Sequence> seqs_;
  size_t sequence_length_; // length of the sequences (required to be the same for all)
};

}
#endif
