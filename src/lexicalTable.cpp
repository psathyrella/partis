#include "lexicalTable.h"
namespace StochHMM {
  
// ----------------------------------------------------------------------------------------
lexicalTable::lexicalTable() {
  logProb = new(nothrow) vector<vector<double> >;
  prob = new(nothrow) vector<vector<double> >;
  counts = new(nothrow) vector<vector<double> >;
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
