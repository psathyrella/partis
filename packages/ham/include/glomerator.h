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

string SeqStr(vector<Sequence> &seqs, string delimiter = " ");
string SeqNameStr(vector<Sequence> &seqs, string delimiter = " ");

// ----------------------------------------------------------------------------------------
class Query {
public:
  string name_;
  vector<Sequence> seqs_;
  KBounds kbounds_;
  vector<string> only_genes_;
  double mean_mute_freq_;
  vector<int> cyst_positions_;
  pair<string, string> parents_;  // queries that were joined to make this
};

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track);
  ~Glomerator();
  void Cluster();
  double LogProbOfPartition(Partition &clusters);
  void Merge(ClusterPath *path, smc::rng *rgen=nullptr);

  // Return the next (i.e. the <i_initial_partition_>th) initial partition in the list of initial partitions, and increment <i_initial_partition_>.
  // Also sets arguments <initial_path_index> and <logweight> to correspond to the returned partition.
  Partition GetAnInitialPartition(int &initial_path_index, double &logweight);

  void WritePartitions(vector<ClusterPath> &paths);
private:
  void ReadCachedLogProbs();
  void GetSoloLogProb(string key);
  void PrintPartition(Partition &clusters, string extrastr);
  void WriteCachedLogProbs();
  double NaiveHammingFraction(string key_a, string key_b);
  int HammingDistance(string seq_a, string seq_b);  // dammit I don't like having two functions, but *@!(#ing Sequence class is not stl safe.
  int HammingDistance(Sequence seq_a, Sequence seq_b);
  void GetNaiveSeq(string key);
  void GetLogProb(string name, vector<Sequence> &seqs, KBounds &kbounds, vector<string> &only_genes, double mean_mute_freq);
  vector<Sequence> MergeSeqVectors(string name_a, string name_b);
  bool SameLength(vector<Sequence> &seqs);
  Query GetMergedQuery(string name_a, string name_b);
  string JoinNames(string name1, string name2);
  // string JoinNameStrings(vector<Sequence> &strlist, string delimiter=":");
  string JoinSeqStrings(vector<Sequence> &strlist, string delimiter=":");
  Query *ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen);

  // input info
  Track *track_;
  Args *args_;
  DPHandler vtb_dph_, fwd_dph_;
  ofstream ofs_;

  vector<Partition> initial_partitions_;
  vector<double> initial_logprobs_;
  vector<double> initial_logweights_;

  int i_initial_partition_;  // index of the next inital paritition to grab (for smc stuff)

  map<string, vector<Sequence> > seq_info_;  // NOTE it would be more memory-efficient to just keep track of vectors of keys here, and have Glomerator keep all the actual info
  map<string, vector<string> > only_genes_;
  map<string, KBounds> kbinfo_;
  map<string, float> mute_freqs_;  // overall mute freq for single sequences, mean overall mute freq for n-sets of sequences

  // NOTE they keys for these two maps are colon-separated lists of *query* *sequences*, whereas all the other maps are of query names. This is because we need logprobs and naive seqs for each truncation length
  // NOTE also that I don't keep track of the order, which I kinda should do since I might be calculating some things twice
  map<string, double> log_probs_;  // includes cached info from previous runs
  map<string, Sequence> naive_seqs_;  // includes cached info from previous runs
  map<string, string> errors_;

  int n_fwd_cached_, n_fwd_calculated_, n_vtb_cached_, n_vtb_calculated_;
};

}
#endif
