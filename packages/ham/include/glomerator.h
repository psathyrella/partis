#ifndef HAM_GLOMERATOR_H
#define HAM_GLOMERATOR_H

#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <functional>
#include <pthread.h>

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
  Query() {}
  Query(string name, vector<Sequence*> seqs, bool seed_missing, vector<string> only_genes, KBounds kbounds, float mute_freq, size_t cdr3_length, string p1="", string p2="") :
    name_(name),
    seqs_(seqs),
    seed_missing_(seed_missing),
    only_genes_(only_genes),
    kbounds_(kbounds),
    mute_freq_(mute_freq),
    cdr3_length_(cdr3_length)
  {
    // if(cdr3_length > 300)
    //   throw runtime_error("cdr3 length too big " + to_string(cdr3_length_) + " for " + name + "\n");
    if(p1 != "" and p2 != "")
      parents_ = pair<string, string>(p1, p2);
    for(auto *pseq : seqs)
      if(pseq == nullptr)
	throw runtime_error("null sequence pointer passed to Query constructor for " + name);
  }

  string name_;
  vector<Sequence*> seqs_;
  bool seed_missing_;
  vector<string> only_genes_;
  KBounds kbounds_;
  float mute_freq_;
  size_t cdr3_length_;
  pair<string, string> parents_;  // queries that were joined to make this
};

// ----------------------------------------------------------------------------------------
class Glomerator {
public:
  Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track);
  ~Glomerator();
  void Cluster();
  double LogProbOfPartition(Partition &clusters, bool debug=false);
  void Merge(ClusterPath *path);

  void CacheNaiveSeqs();
  // NOTE don't remove these (yet, at least)
  // ClusterPair GetClustersToMergeForNaiveSeqGlomerate(set<vector<string> > &clusters, int max_per_cluster, bool merge_whatever_you_got);
  // void PrintClusterSizes(set<vector<string> > &clusters);
  // ClusterPair GetSmallBigClusters(set<vector<string> > &clusters);
  // void NaiveSeqGlomerate(int n_clusters);
  
  // Return the next (i.e. the <i_initial_partition_>th) initial partition in the list of initial partitions, and increment <i_initial_partition_>.
  // Also sets arguments <initial_path_index> and <logweight> to correspond to the returned partition.
  // Partition GetAnInitialPartition(int &initial_path_index, double &logweight);

  void WritePartitions(ClusterPath &cp);
  void WriteAnnotations(ClusterPath &cp);
private:
  void ReadCacheFile();
  void WriteCacheLine(ofstream &ofs, string query);
  void WriteCacheFile();

  void PrintPartition(Partition &clusters, string extrastr);
  string CacheSizeString();
  string FinalString(bool newline=false);
  string GetStatusStr(time_t current_time);
  void WriteStatus();  // write some progress info to file

  string ParentalString(pair<string, string> *parents);
  int CountMembers(string namestr, bool exclude_extra_seeds=false);
  unsigned LargestClusterSize(Partition &partition);
  string ClusterSizeString(Partition *partition);
  string JoinNames(string name1, string name2, string delimiter=":");
  string JoinNameStrings(vector<Sequence*> &strlist, string delimiter=":");
  string JoinSeqStrings(vector<Sequence*> &strlist, string delimiter=":");
  string PrintStr(string queries);
  bool SeedMissing(string queries);

  double CalculateHfrac(string &seq_a, string &seq_b);
  double NaiveHfrac(string key_a, string key_b);

  string ChooseSubsetOfNames(string queries, int n_max);
  string GetNaiveSeqNameToCalculate(string actual_queries);  // convert between the actual queries/key we're interested in and the one we're going to calculate
  string GetLogProbNameToCalculate(string queries, int n_max);
  pair<string, string> GetLogProbPairOfNamesToCalculate(string actual_queries, pair<string, string> actual_parents);  // convert between the actual queries/key we're interested in and the one we're going to calculate
  bool FirstParentMuchBigger(string queries, string queries_other, int nmax);
  string FindNaiveSeqNameReplace(pair<string, string> *parents);
  string &GetNaiveSeq(string key, pair<string, string> *parents=nullptr);
  // double NormFactor(string name);
  double GetLogProb(string queries);
  double GetLogProbRatio(string key_a, string key_b);
  string CalculateNaiveSeq(string key, RecoEvent *event=nullptr);
  double CalculateLogProb(string queries);

  bool check_cache(string queries) {
    if(cachefo_.find(queries) != cachefo_.end())
      return true;
    else if(tmp_cachefo_.find(queries) != tmp_cachefo_.end())
      return true;
    else
      throw false;
  }

  Query &cachefo(string queries);

  bool SameLength(vector<Sequence*> &seqs, bool debug=false);
  void AddFailedQuery(string queries, string error_str);
  void UpdateLogProbTranslationsForAsymetrics(Query &qmerge);
  vector<Sequence*> GetSeqs(string query);
  void MoveSubsetsFromTmpCache(string query);
  void CopyToPermanentCache(string translated_query, string superquery);
  Query &GetMergedQuery(string name_a, string name_b);

  bool LikelihoodRatioTooSmall(double lratio, int candidate_cluster_size);
  Partition GetSeededClusters(Partition &partition);
  pair<double, Query> FindHfracMerge(ClusterPath *path);
  pair<double, Query> FindLRatioMerge(ClusterPath *path);
  pair<double, Query> *ChooseRandomMerge(vector<pair<double, Query> > &potential_merges);

  Track *track_;
  Args *args_;
  GermLines &gl_;
  HMMHolder &hmms_;
  ofstream ofs_;

  Partition initial_partition_;

  map<string, string> naive_seq_name_translations_;
  map<string, pair<string, string> > logprob_name_translations_;
  map<string, string> logprob_asymetric_translations_;
  map<string, string> name_subsets_;

  map<string, Sequence> single_seqs_;  // only place that we keep the actual sequences (rather than pointers/references)
  map<string, Query> single_seq_cachefo_;  // keep some (approximate) single-sequence info to help us build missing cache entries
  map<string, Query> cachefo_;  // cache info for clusters we've actually merged
  map<string, Query> tmp_cachefo_;  // cache info for clusters we're only considering merging

  // These all include cached info from previous runs
  map<string, double> log_probs_;  
  map<string, double> naive_hfracs_;  // NOTE since this uses the joint key, it assumes there's only *one* way to get to a given cluster (this is similar to, but not quite the same as, the situation for log probs and naive seqs)
  map<string, double> lratios_;
  map<string, string> naive_seqs_;
  map<string, string> errors_;

  set<string> failed_queries_;

  set<string> initial_log_probs_, initial_naive_hfracs_, initial_naive_seqs_;  // keep track of the ones we read from the initial cache file so we can write only the new ones to the output cache file

  int n_fwd_calculated_, n_vtb_calculated_, n_hfrac_calculated_, n_hfrac_merges_, n_lratio_merges_;

  double asym_factor_;

  bool force_merge_;  // this gets set to true if args_->n_final_clusters() is set, and we've got to keep going past the most likely partition in order to get down to the requested number of clusters

  Partition *current_partition_;  // (a.t.m. only used for writing to status file)
  time_t last_status_write_time_;  // last time that we wrote our progress to a file
  FILE *progress_file_;

  string empty_string_;  // ok this is horrible, but i need to return a reference to an empty string from GetNaiveSeq(), and it's been so long since i edited this code i can't figure out a better way to do it than this
};

}
#endif
