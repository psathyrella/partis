#include "lexicaltable.h"
namespace ham {

// ----------------------------------------------------------------------------------------
LexicalTable::LexicalTable() : log_probs_(NULL) {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::Init() {
  log_probs_ = new vector<vector<double> >;
}

// ----------------------------------------------------------------------------------------
LexicalTable::~LexicalTable() {
  delete log_probs_;
}

// ----------------------------------------------------------------------------------------
void LexicalTable::AddColumn(vector<double> logprobs) {
  log_probs_->push_back(logprobs);
}

}
