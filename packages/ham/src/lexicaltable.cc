#include "lexicaltable.h"
namespace ham {

// ----------------------------------------------------------------------------------------
LexicalTable::LexicalTable() {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::Init() {
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
    throw runtime_error("ERROR bad normalization after replacement " + to_string(total) + "\n");
}

// ----------------------------------------------------------------------------------------
void LexicalTable::UnReplaceLogProbs() {
  log_probs_ = original_log_probs_;
}

// ----------------------------------------------------------------------------------------
LexicalTable::~LexicalTable() {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::AddColumn(vector<double> logprobs) {
  assert(log_probs_.size() == 0);  // a.t.m. we only want to support one column, i.e. no joint emission
  log_probs_.push_back(logprobs);
}

}
