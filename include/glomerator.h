#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>

#include "smctc.hh"
#include "args.h"
#include "dphandler.h"
#include "clusterpath.h"
#include "text.h"

using namespace std;
namespace ham {

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track);
  void Cluster();
  double LogProbOfPartition(Partition &clusters);
  smc::particle<ClusterPath> SMCInit(smc::rng *rgen);
  void SMCMove(long time, smc::particle<ClusterPath> &particle, smc::rng *rgen);
private:
  bool AllFinished();
  void ReadCachedLogProbs(Track *track);
  Partition GetClusterList(map<string, Sequences> &partinfo);
  void GetSoloLogProb(string key);
  void PrintPartition(Partition &clusters, string extrastr);
  void WriteCachedLogProbs();
  void WritePartitions();
  int NaiveHammingDistance(string key_a, string key_b);
  int HammingDistance(string seq_a, string seq_b);  // dammit I don't like having two functions, but *@!(#ing Sequence class is not stl safe.
  int HammingDistance(Sequence seq_a, Sequence seq_b);
  // int MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b);
  void GetNaiveSeq(string key, ClusterPath &path);
  void GetLogProb(DPHandler &dph, string name, Sequences &seqs, KBounds &kbounds);
  string JoinNames(string name1, string name2);
  void Merge(ClusterPath &path);

  // input info
  HMMHolder hmms_;
  GermLines gl_;
  Args *args_;
  ofstream ofs_;

  Partition initial_partition_;
  double initial_logprob_;

  map<string, Sequences> info_;  // NOTE it would be more memory-efficient to just keep track of vectors of keys here, and have Glomerator keep all the actual info
  map<string, vector<string> > only_genes_;
  map<string, KBounds> kbinfo_;

  map<string, double> log_probs_;  // includes cached info from previous runs
  map<string, string> naive_seqs_;  // includes cached info from previous runs
  map<string, string> errors_;
  smc::sampler<ClusterPath> sampler_;
  smc::moveset<ClusterPath> moveset_;
  vector<ClusterPath> clusterpaths_;

};

}
#endif
