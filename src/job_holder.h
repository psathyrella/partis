#ifndef JOB_HOLDER_H
#define JOB_HOLDER_H

#include <map>
#include <string>
#include "hmm.h"
#include "sequences.h"
#include "germlines.h"

using namespace std;
using namespace StochHMM;

// ----------------------------------------------------------------------------------------
class JobHolder {
public:
  JobHolder(string hmmtype, string hmm_dir, string seqfname);
  void InitJobs(size_t k_v, size_t k_d);
  void GetNextHMM();

  sequences *current_seqs_;
  model current_hmm_;
  bool finished_;
private:
  map<string,sequences> GetSubSeqs(size_t k_v, size_t k_d);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  string sanitize_name(string gene_name);  // should put this somewhere else

  string hmm_dir_;  // location of .hmm files
  GermLines gl_;
  vector<string> regions_;
  track track_;
  sequences seqs_;
  map<string,sequences> subseqs_;
  vector<string>::iterator i_current_region_;  // region and position of the *next* job we will pass out with GetNextJob()
  vector<string>::iterator i_current_gene_;
};
#endif
