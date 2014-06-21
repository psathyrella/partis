#ifndef JOB_HOLDER_H
#define JOB_HOLDER_H

#include <map>
#include <string>
#include "StochHMMlib.h"
#include "StochHMM_usage.h"
#include "germlines.h"

using namespace std;
using namespace StochHMM;

// ----------------------------------------------------------------------------------------
class JobHolder {
public:
  JobHolder(string hmmtype, string hmm_dir, string seqfname, size_t n_max_versions=0);
  void Run(size_t k_v, size_t k_d, string algorithm);

private:
  map<string,sequences> GetSubSeqs(size_t k_v, size_t k_d);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  size_t GetInsertLength(vector<string> labels);
  size_t GetErosionLength(string side, vector<string> path_labels, string gene_name);

  string hmm_dir_;  // location of .hmm files
  size_t n_max_versions_;    // only look at the first n gene versions (speeds things up for testing)
  GermLines gl_;
  vector<string> regions_;
  track track_;
  sequences seqs_;
  vector<string>::iterator i_current_region_;  // region and position of the *next* job we will pass out with GetNextJob()
  vector<string>::iterator i_current_gene_;
};
#endif
