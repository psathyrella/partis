#include "lexicaltable.h"
namespace ham {

// ----------------------------------------------------------------------------------------
LexicalTable::LexicalTable() {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::Init() {  // TODO doesn't do anything
}

// ----------------------------------------------------------------------------------------
void LexicalTable::ReplaceLogProbs(vector<vector<double> > new_log_probs) {
  assert(log_probs_.size() == new_log_probs.size());
  original_log_probs_ = log_probs_;
  log_probs_ = new_log_probs;
  double total(0.0); // make sure things add to 1.0
  size_t icol(0);
  for(size_t irow=0; irow<log_probs_[icol].size(); ++irow)  { // or icol, it's kind arbitrary
    double prob = exp(log_probs_[icol][irow]);
    total += prob;
  }
  if(fabs(total - 1.0) >= EPS)  // make sure transition probs sum to 1.0
    throw runtime_error("ERROR bad normalization after replacement " + to_string(total) + " in LexicalTable::ReplaceLogProbs()\n");
}

// ----------------------------------------------------------------------------------------
void LexicalTable::UnReplaceLogProbs() {
  assert(original_log_probs_.size() == log_probs_.size());
  log_probs_ = original_log_probs_;
}

// ----------------------------------------------------------------------------------------
double LexicalTable::LogProb(Sequence *seq, size_t pos) {
  assert(pos < (*seq).size());
  if(log_probs_.size() != 1) {
    cout << "lp s " << log_probs_.size() << endl;
    assert(0);
  }
  assert((*seq)[pos] < log_probs_[0].size());
  return log_probs_[0][(*seq)[pos]];
}

// ----------------------------------------------------------------------------------------
LexicalTable::~LexicalTable() {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::AddColumn(vector<double> logprobs) {
  assert(log_probs_.size() == 0);  // a.t.m. we only want to support one column, i.e. no joint emission
  log_probs_.push_back(logprobs);
  assert(log_probs_.size() == 1); // TODO remove this
}

}
