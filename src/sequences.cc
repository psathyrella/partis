#include "sequences.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Sequence::Sequence(string name, string undigitized, Track* trk):  //, vector<uint8_t> *seq, int pos, int len):
  name_(name),
  undigitized_(undigitized),
  track_(trk)
{
  ClearWhitespace("\n", &undigitized_);
  // EDIT DAMMIT THIS ISN'T WORKING. SCREW OPTIMIZING
  // ClearWhitespace("\n", &undigitized);
  // if (pos > 0) {
  //   undigitized_ = undigitized.substr(pos, len);
  //   seq_ = new vector<uint8_t>(len);
  //   // for(int i=pos; i<pos+len; ++i)
  //   //   (*seq_)[i-pos] = (*seq)[i];
  //   for(int i=0; i<len; ++i)
  //     (*seq_)[i] = (*seq)[pos+i];
  //   // cout << undigitized << endl;
  //   // cout << undigitized_ << endl;
  //   // assert(0);
  // } else {
  //   undigitized_ = undigitized;
  seq_ = new vector<uint8_t>(undigitized_.size());
  for(size_t i = 0; i < undigitized_.size(); ++i) {
    uint8_t symbl = track_->symbol_index(undigitized_.substr(i, 1));
    (*seq_)[i] = symbl;
  }
  // }
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(const Sequence& rhs) {
  name_ = rhs.name_;
  header_  = rhs.header_;
  track_  = rhs.track_;
  undigitized_ = rhs.undigitized_;
  seq_ = rhs.seq_;
  if(rhs.seq_)
    seq_ = new vector<uint8_t>(*rhs.seq_);
}

// ----------------------------------------------------------------------------------------
Sequence::~Sequence() {
  delete seq_;
}

// ----------------------------------------------------------------------------------------
//! \return a copy of the sequence from <pos> of size <len>
Sequence Sequence::GetSubSequence(size_t pos, size_t len) {
  assert(pos < undigitized_.size());
  assert(len < undigitized_.size());
  assert(pos + len <= undigitized_.size()); // arg. if pos+len overflows this still passes
  string subseq_str = undigitized_.substr(pos, len);
  return Sequence(name_, subseq_str, track_);

  // return Sequence(name_, undigitized_, track_, seq_, pos, len);

  // vector<uint8_t> digisubseq(len);
  // for (size_t i=pos; i<len; ++i)
  //   digisubseq[i] = seq_->at(i);
  // return Sequence(name_, undigitized_.substr(pos, len), track_, &digisubseq);
}

// ----------------------------------------------------------------------------------------
void Sequence::Print(string separator) {
  if(!header_.empty())
    cout << header_ << endl;
  for(auto & symbol : undigitized_)
    cout << symbol << separator;
  cout << endl;
}

// ****************************************************************************************
// ----------------------------------------------------------------------------------------
void Sequences::Print() {
  for(size_t i = 0; i < seqs_.size(); ++i) {
    Track* trk = seqs_[i].track();
    cout << ">" << trk->name() << endl;
    seqs_[i].Print();
  }
}

// ----------------------------------------------------------------------------------------
//! \return return the collection of subsequences from <pos> with size <len>
Sequences Sequences::GetSubSequences(size_t pos, size_t len) {
  Sequences new_seqs;  // init with *zero* seqs because we push back below
  for(size_t is = 0; is < seqs_.size(); is++)
    new_seqs.AddSeq(Sequence(seqs_[is].GetSubSequence(pos, len)));
  return new_seqs;
}

// ----------------------------------------------------------------------------------------
void Sequences::AddSeq(Sequence sq) {
  if(sequence_length_ == 0) {   // make sure <sq> has the same length as any previously added sequences
    sequence_length_ = sq.size();
  } else {
    assert(sq.size() == sequence_length_);  // all sequences must have the same length
  }
  seqs_.push_back(sq);  // NOTE we now own this sequence, i.e. we will delete it when we die
}

}
