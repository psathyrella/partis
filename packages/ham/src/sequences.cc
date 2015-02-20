#include "sequences.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Track* trk, string name, string &undigitized):
  name_(name),
  track_(trk) {
  undigitized_ = undigitized;
  ClearWhitespace("\n", &undigitized_);
  Digitize();
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Track* trk, string name, string &undigitized, size_t pos, size_t len):
  name_(name),
  track_(trk) {
  undigitized_ = undigitized.substr(pos, len);  // <len> better be greater than zero
  ClearWhitespace("\n", &undigitized_);
  Digitize();
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(Sequence &rhs, size_t pos, size_t len) {
  name_ = rhs.name_;
  header_  = rhs.header_;
  track_  = rhs.track_;
  undigitized_ = rhs.undigitized().substr(pos, len);  // <len> better be greater than zero
  seq_ = new vector<uint8_t>(rhs.seq_->begin() + pos, rhs.seq_->begin() + pos + len);
}

// ----------------------------------------------------------------------------------------
Sequence::Sequence(const Sequence &rhs) {
  name_ = rhs.name_;
  header_  = rhs.header_;
  track_  = rhs.track_;
  undigitized_ = rhs.undigitized();
  seq_ = new vector<uint8_t>(*rhs.seq_);
}

// ----------------------------------------------------------------------------------------
void Sequence::Digitize() {
  seq_ = new vector<uint8_t>(undigitized_.size());
  for(size_t i = 0; i < undigitized_.size(); ++i) {
    string symbol = undigitized_.substr(i, 1);
    (*seq_)[i] = track_->symbol_index(symbol);
  }
}

// ----------------------------------------------------------------------------------------
Sequence::~Sequence() {
  delete seq_;
}

// ----------------------------------------------------------------------------------------
Sequences::Sequences(Sequences &seqs, size_t pos, size_t len) : sequence_length_(0) {
  for(auto & seq : seqs.seqs_)
    AddSeq(Sequence(seq, pos, len));
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
void Sequences::AddSeq(Sequence sq) {
  if(n_seqs() == 0) {  // if this is the first sequence, set <sequence_length_>
    sequence_length_ = sq.size();
  } else {
    if(sq.size() != sequence_length_)  // all sequences must have the same length
      throw runtime_error("ERROR sequences must all have the same length, but got " + to_string(sq.size()) + " and " + to_string(sequence_length_));
  }
  seqs_.push_back(sq);  // NOTE we now own this sequence, i.e. we will delete it when we die
}

}
