#include "lexicaltable.h"
namespace ham {

// ----------------------------------------------------------------------------------------
LexicalTable::LexicalTable() : track_(nullptr) {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::Init(Track *track) {
  track_ = track;
}

// ----------------------------------------------------------------------------------------
void LexicalTable::ReplaceLogProbs(vector<double> new_log_probs) {
  assert(log_probs_.size() == new_log_probs.size());
  original_log_probs_ = log_probs_;
  log_probs_ = new_log_probs;
  double total(0.0); // make sure things add to 1.0
  for(size_t ip=0; ip<log_probs_.size(); ++ip)  {
    double prob = exp(log_probs_[ip]);
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
double LexicalTable::LogProb(Sequence *seq, size_t pos) {  // todo profile and improve checking
  assert(pos < (*seq).size());
  assert((*seq)[pos] < log_probs_.size());
  return log_probs_[(*seq)[pos]];
}

// ----------------------------------------------------------------------------------------
LexicalTable::~LexicalTable() {
}

}
