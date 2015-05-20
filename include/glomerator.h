#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

#include "smctc.hh"
#include "args.h"
#include "dphandler.h"
#include "clusterpath.h"
#include "text.h"

using namespace std;
namespace ham {

// ----------------------------------------------------------------------------------------
class Query {
public:
  string name_;
  Sequences seqs_;
  KBounds kbounds_;
  vector<string> only_genes_;
  double mean_mute_freq_;
  pair<string, string> parents_;  // queries that were joined to make this
};

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track);
  ~Glomerator();
  void Cluster();
  double LogProbOfPartition(Partition &clusters);
  void Merge(ClusterPath *path, smc::rng *rgen=nullptr);
  Partition GetAnInitialPartition(int &tmpi, double &logweight);  // get the next initial partition in the list of initial partitions, and remove it from the list
  void WritePartitions(vector<ClusterPath> &paths);
private:
  void ReadCachedLogProbs(Track *track);
  Partition GetPartitionFromMap(map<string, Sequences> &partinfo);
  void GetSoloLogProb(string key);
  void PrintPartition(Partition &clusters, string extrastr);
  void WriteCachedLogProbs();
  int NaiveHammingDistance(string key_a, string key_b);
  int HammingDistance(string seq_a, string seq_b);  // dammit I don't like having two functions, but *@!(#ing Sequence class is not stl safe.
  int HammingDistance(Sequence seq_a, Sequence seq_b);
  // int MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b);
  void GetNaiveSeq(string key);
  void GetLogProb(DPHandler &dph, string name, Sequences &seqs, KBounds &kbounds, vector<string> &only_genes, double mean_mute_freq);
  Query GetMergedQuery(string name_a, string name_b);
  string JoinNames(string name1, string name2);
  Query *ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen);

  // input info
  HMMHolder hmms_;
  GermLines gl_;
  Args *args_;
  ofstream ofs_;

  vector<Partition> initial_partitions_;
  vector<double> initial_logprobs_;
  vector<double> initial_logweights_;

  int i_initial_partition_;  // index of the next inital paritition to grab (for smc stuff)

  map<string, Sequences> info_;  // NOTE it would be more memory-efficient to just keep track of vectors of keys here, and have Glomerator keep all the actual info
  map<string, vector<string> > only_genes_;
  map<string, KBounds> kbinfo_;
  map<string, float> mute_freqs_;  // overall mute freq for single sequences, mean overall mute freq for n-sets of sequences

  map<string, double> log_probs_;  // includes cached info from previous runs
  map<string, string> naive_seqs_;  // includes cached info from previous runs
  map<string, string> errors_;

  int n_cached_, n_calculated_;
};

}
#endif
