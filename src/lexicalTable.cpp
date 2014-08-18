#include "lexicalTable.h"
namespace StochHMM {
  
// ----------------------------------------------------------------------------------------
lexicalTable::lexicalTable() : logProb(NULL), prob(NULL), counts(NULL) {
}

// ----------------------------------------------------------------------------------------
void lexicalTable::init() {
  logProb = new vector<vector<double> >;
  prob = new vector<vector<double> >;
  counts = new vector<vector<double> >;
}
  
// ----------------------------------------------------------------------------------------
lexicalTable::~lexicalTable() {
  delete logProb;
  delete prob;
  delete counts;
}

// ----------------------------------------------------------------------------------------
void lexicalTable::AddColumn(vector<double> logprobs) {
  logProb->push_back(logprobs);
  prob->push_back(get_exp_vector(logprobs));
}
  
}
