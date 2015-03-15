#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>

#include "args.h"
#include "jobholder.h"
#include "text.h"

using namespace std;
namespace ham {

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args);
  void Cluster();
private:
  void ReadCachedLogProbs(string fname);
  vector<string> GetClusterList(map<string, Sequences> &partinfo);
  void GetSoloLogProb(string key);
  double LogProbOfPartition(vector<string> &clusters);
  void PrintPartition(vector<string> &clusters, string extrastr);
  void WriteCachedLogProbs();
  void WritePartitions();
  int MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b);
  void GetResult(JobHolder &jh, string name, Sequences &seqs, KBounds &kbounds, string &errors);
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
  map<string, double> log_probs_;
  vector<pair<double, vector<string> > > list_of_partitions_;
};

}
#endif
