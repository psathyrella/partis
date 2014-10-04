#include "lexicaltable.h"
namespace ham {
  
// ----------------------------------------------------------------------------------------
LexicalTable::LexicalTable() : prob(NULL), counts(NULL), logProb(NULL) {
}

// ----------------------------------------------------------------------------------------
void LexicalTable::init() {
  prob = new vector<vector<double> >;
  counts = new vector<vector<double> >;
  logProb = new vector<vector<double> >;
}
  
// ----------------------------------------------------------------------------------------
LexicalTable::~LexicalTable() {
  delete logProb;
  delete prob;
  delete counts;
}

// ----------------------------------------------------------------------------------------
void LexicalTable::AddColumn(vector<double> logprobs) {
  logProb->push_back(logprobs);
  prob->push_back(get_exp_vector(logprobs));
}
  
}
