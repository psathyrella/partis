#ifndef JOB_HOLDER_H
#define JOB_HOLDER_H

#include <map>
#include <string>
#include "hmm.h"
#include "sequences.h"

using namespace std;
using namespace StochHMM;

class JobHolder {
public:
  JobHolder(string input_dir, vector<string> regions, string seqfname);
  map<string,sequences> GetSubSeqs(size_t k_v, size_t k_d);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  model *hmm(string region) { assert(hmms_.find(region) != hmms_.end()); return &hmms_[region]; };
private:
  vector<string> regions_;
  map<string,model> hmms_;
  map<string,sequences> seqs_;
  StateFuncs default_functions_;  // automatically initialize all the Univariate and Multivariate PDFs. Only included for backward compatibility
};
#endif
