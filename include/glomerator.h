#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <ctime>
#include <algorithm>

#include "smctc.hh"
#include "args.h"
#include "dphandler.h"
#include "clusterpath.h"
#include "text.h"

using namespace std;
namespace ham {

typedef pair<vector<string>, vector<string> > ClusterPair;

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
  double LogProbOfPartition(Partition &clusters, bool debug=false);
  void Merge(ClusterPath *path, smc::rng *rgen=nullptr);

  void CacheNaiveSeqs();
  void PrintClusterSizes(set<vector<string> > &clusters);
  ClusterPair GetClustersToMerge(set<vector<string> > &clusters, int max_per_cluster, bool merge_whatever_you_got);
  ClusterPair GetSmallBigClusters(set<vector<string> > &clusters);
  void NaiveSeqGlomerate(int n_clusters);
  
  // Return the next (i.e. the <i_initial_partition_>th) initial partition in the list of initial partitions, and increment <i_initial_partition_>.
  // Also sets arguments <initial_path_index> and <logweight> to correspond to the returned partition.
  Partition GetAnInitialPartition(int &initial_path_index, double &logweight);

  void WritePartitions(vector<ClusterPath> &paths);
  void WriteAnnotations(vector<ClusterPath> &paths);
private:
  void ReadCachedLogProbs();
  void GetSoloLogProb(string key);
  void PrintPartition(Partition &clusters, string extrastr);
  void WriteCacheLine(ofstream &ofs, string query);
  void WriteCachedLogProbs();
  string ParentalString(pair<string, string> *parents);
  int CountMembers(string namestr);
  string ClusterSizeString(ClusterPath *path);
  void WriteStatus(ClusterPath *path);  // write some progress info to file
  double NaiveHammingFraction(string key_a, string key_b);
  double HammingFraction(Sequence seq_a, Sequence seq_b);
  void GetNaiveSeq(string key, pair<string, string> *parents=nullptr);
  void GetLogProb(string name, vector<Sequence> &seqs, KBounds &kbounds, vector<string> &only_genes, double mean_mute_freq);
  vector<Sequence> MergeSeqVectors(string name_a, string name_b);
  bool SameLength(vector<Sequence> &seqs, bool debug=false);
  Query GetMergedQuery(string name_a, string name_b);
  string JoinNames(string name1, string name2);
  string JoinNameStrings(vector<Sequence> &strlist, string delimiter=":");
  string JoinSeqStrings(vector<Sequence> &strlist, string delimiter=":");
  Query ChooseMerge(ClusterPath *path, smc::rng *rgen, double *chosen_lratio);
  Query *ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen);

  Track *track_;
  Args *args_;
  DPHandler vtb_dph_, fwd_dph_;
  ofstream ofs_;

  vector<Partition> initial_partitions_;
  vector<double> initial_logprobs_;
  vector<double> initial_logweights_;

  map<string, double> hamming_fractions_;  // cached hamming fractions NOTE key is query names, *not* sequences

  int i_initial_partition_;  // index of the next inital paritition to grab (for smc stuff)

  map<string, vector<Sequence> > seq_info_;  // NOTE it would be more memory-efficient to just keep track of vectors of keys here, and have Glomerator keep all the actual info
  map<string, vector<string> > only_genes_;
  map<string, KBounds> kbinfo_;
  map<string, float> mute_freqs_;  // overall mute freq for single sequences, mean overall mute freq for n-sets of sequences

  // NOTE they keys for these two maps are colon-separated lists of *query* *sequences*, whereas all the other maps are of query names. This is because we need logprobs and naive seqs for each truncation length
  // NOTE also that I don't keep track of the order, which I kinda should do since I might be calculating some things twice.
  // These all include cached info from previous runs
  map<string, double> log_probs_;  
  map<string, double> naive_hfracs_;  // NOTE since this uses the joint key, it assumes there's only *one* way to get to a given cluster
  map<string, Sequence> naive_seqs_;
  map<string, RecoEvent> events_;  // annotations corresponding to the naive seqs. NOTE keeping it separate, at least for now, since I only want the full annotations for the final partition, but I need naive seqs for loads and loads of groups of sequences
  map<string, string> errors_;
  map<string, double> naive_hamming_fractions_;

  set<string> initial_log_probs_, initial_naive_hfracs_, initial_naive_seqs_;  // keep track of the ones we read from the initial cache file so we can write only the new ones to the output cache file

  int n_fwd_calculated_, n_vtb_calculated_, n_hfrac_calculated_, n_hamming_merged_;

  time_t last_status_write_time_;  // last time that we wrote our progress to a file
  FILE *progress_file_;
};

}
#endif
