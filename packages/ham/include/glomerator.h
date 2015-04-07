#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>

#include "smctc.hh"
#include "args.h"
#include "jobholder.h"
#include "text.h"

using namespace std;
namespace ham {

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track);
  void Cluster();
private:
  void ReadCachedLogProbs(Track *track);
  vector<string> GetClusterList(map<string, Sequences> &partinfo);
  void GetSoloLogProb(string key);
  double LogProbOfPartition(vector<string> &clusters);
  void PrintPartition(vector<string> &clusters, string extrastr);
  void WriteCachedLogProbs();
  void WritePartitions();
  int NaiveHammingDistance(string key_a, string key_b);
  int HammingDistance(string seq_a, string seq_b);  // dammit I don't like having two functions, but *@!(#ing Sequence class is not stl safe.
  int HammingDistance(Sequence seq_a, Sequence seq_b);
  int MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b);
  void GetNaiveSeq(string key);
  void GetLogProb(JobHolder &jh, string name, Sequences &seqs, KBounds &kbounds);
  void Merge();

  // input info
  HMMHolder hmms_;
  GermLines gl_;
  Args *args_;
  ofstream ofs_;

  // current cluster info (modified at each step)
  map<string, Sequences> info_;
  map<string, vector<string> > only_genes_;
  map<string, KBounds> kbinfo_;
  double max_log_prob_of_partition_;
  vector<string> best_partition_;
  bool finished_;

  // cumulative information
  map<string, double> log_probs_;  // includes cached info from previous runs
  map<string, string> naive_seqs_;  // includes cached info from previous runs
  map<string, string> errors_;
  vector<pair<double, vector<string> > > list_of_partitions_;
};

}
#endif
