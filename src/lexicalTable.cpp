#include "lexicalTable.h"
namespace StochHMM {
  
// ----------------------------------------------------------------------------------------
lexicalTable::lexicalTable() : prob(NULL), counts(NULL), logProb(NULL) {
}

// ----------------------------------------------------------------------------------------
void lexicalTable::init() {
  prob = new vector<vector<double> >;
  counts = new vector<vector<double> >;
  logProb = new vector<vector<double> >;
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
