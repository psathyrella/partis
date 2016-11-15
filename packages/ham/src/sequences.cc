#include "sequences.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Sequence::Sequence() : track_(nullptr)
{
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Track* trk, string name, string &undigitized):
  name_(name),
  track_(trk),
  seqq_(undigitized.size())
{
  undigitized_ = string(undigitized);
  ClearWhitespace("\n", &undigitized_);
  Digitize();
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Track* trk, string name, string &undigitized, size_t pos, size_t len):
  name_(name),
  track_(trk),
  seqq_(undigitized.size())
{
  CheckPosLen(name, undigitized, pos, len);
  undigitized_ = undigitized.substr(pos, len);  // <len> better be greater than zero
  ClearWhitespace("\n", &undigitized_);
  Digitize();
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Sequence &rhs, size_t pos, size_t len) {
  name_ = rhs.name_;
  CheckPosLen(name_, rhs.undigitized(), pos, len);
  header_  = rhs.header_;
  track_  = rhs.track_;
  undigitized_ = rhs.undigitized().substr(pos, len);  // <len> better be greater than zero
  seqq_ = vector<uint8_t>(rhs.seqq_.begin() + pos, rhs.seqq_.begin() + pos + len);
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(const Sequence &rhs) {
  name_ = rhs.name_;
  header_  = rhs.header_;
  undigitized_ = rhs.undigitized_;
  track_  = rhs.track_;
  seqq_ = rhs.seqq_;
}

// ----------------------------------------------------------------------------------------
void Sequence::CheckPosLen(string name, string undigitized, size_t pos, size_t len) {
  // if(pos < 0 || len < 0)  // can't actually be less than zero, but I have great hopes of eventually making them signed
  //   throw runtime_error("invalid pos " + to_string(pos) + " len " + to_string(len) + " for " + name);
  if(pos >= undigitized.size() || pos + len > undigitized.size())
    throw runtime_error("len " + to_string(len) + " too large for " + undigitized + " in " + name);
}

// ----------------------------------------------------------------------------------------
void Sequence::Digitize() {
  assert(seqq_.size() == undigitized_.size());
  for(size_t i = 0; i < undigitized_.size(); ++i) {
    string symbol = undigitized_.substr(i, 1);
    seqq_[i] = track_->symbol_index(symbol);
  }
}

// ----------------------------------------------------------------------------------------
Sequence::~Sequence() {
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
// // ----------------------------------------------------------------------------------------
// Sequences::Sequences(const Sequences &rhs) {
//   for(size_t is=0; is<rhs.n_seqs(); ++is)
//     AddSeq(rhs.GetAtConst(is));
// }

// ----------------------------------------------------------------------------------------
Sequences Sequences::Union(Sequences &otherseqs) {
  Sequences union_seqs;
  set<string> already_added;  // set of sequences names to make sure we don't add the same sequence twice (this is really just a neurotic double check of whatever code is sending us sequences)

  for(size_t is=0; is<seqs_.size(); ++is) {
    if(already_added.count(seqs_[is].name()))
      throw runtime_error("tried to add sequence with name " + seqs_[is].name() + " twice in Sequences::Union");
    else
      already_added.insert(seqs_[is].name());
    union_seqs.AddSeq(seqs_[is]);
  }

  for(size_t is=0; is<otherseqs.n_seqs(); ++is) {
    if(already_added.count(otherseqs[is].name()))
      throw runtime_error("tried to add sequence with name " + otherseqs[is].name() + " twice in Sequences::Union");
    else
      already_added.insert(otherseqs[is].name());
    union_seqs.AddSeq(otherseqs[is]);
  }

  return union_seqs;
}
// ----------------------------------------------------------------------------------------
Sequences::Sequences(Sequences &seqs, size_t pos, size_t len) : sequence_length_(0) {
  for(auto & seq : seqs.seqs_)
    AddSeq(Sequence(seq, pos, len));
}

// // ----------------------------------------------------------------------------------------
// Sequences::Sequences(vector<Sequence> &seqs) {
//   for(auto & seq : seqs)
//     AddSeq(seq);
// }

// ----------------------------------------------------------------------------------------
void Sequences::Print() {
  for(size_t i = 0; i < seqs_.size(); ++i) {
    Track* trk = seqs_[i].track();
    cout << ">" << trk->name() << endl;
    seqs_[i].Print();
  }
}

// ----------------------------------------------------------------------------------------
void Sequences::AddSeq(Sequence sq) {
  if(n_seqs() == 0) {  // if this is the first sequence, set <sequence_length_>
    sequence_length_ = sq.size();
  } else {
    if(sq.size() != sequence_length_)  // all sequences must have the same length
      throw runtime_error("Sequences::AddSeq() sequences must all have the same length, but got " + to_string(sq.size()) + " and " + to_string(sequence_length_));
  }
  seqs_.push_back(sq);  // NOTE we now own this sequence, i.e. we will delete it when we die
}

}
