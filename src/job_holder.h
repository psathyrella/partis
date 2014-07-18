#ifndef JOB_HOLDER_H
#define JOB_HOLDER_H

#include <map>
#include <string>
#include <math.h>
#include <set>
#include "StochHMMlib.h"
#include "StochHMM_usage.h"
#include "germlines.h"

using namespace std;
using namespace StochHMM;

typedef pair<size_t,size_t> KSet;  // pair of k_v,k_d values specifying how to chop up the query sequence into v+insert, d+insert, j
typedef pair<string,string> StrPair;

// ----------------------------------------------------------------------------------------
class HMMHolder {
public:
  HMMHolder(string hmm_dir, size_t n_seqs_per_track):hmm_dir_(hmm_dir),n_seqs_per_track_(n_seqs_per_track) {}
  ~HMMHolder();
  model *Get(string gene);
private:
  string hmm_dir_;
  size_t n_seqs_per_track_;
  GermLines gl_;  // TODO kind of hackey to have a separate one of these in HMMHolder. Then again, I don't think it's really that expensive.
  map<string,model*> hmms_;  // map of gene name to hmm pointer
};

// ----------------------------------------------------------------------------------------
class JobHolder {
public:
  JobHolder(size_t n_seqs_per_track, string algorithm, sequences *seqs, HMMHolder *hmms, string only_genes="");
  ~JobHolder();
  void Run(size_t k_v_start, size_t n_k_v, size_t k_d_start, size_t n_k_d);
  void SetDebug(int debug) { debug_ = debug; };
  void FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score);
  void PushBackRecoEvent(KSet kset, map<string,string> &best_genes, double score);
  RecoEvent FillRecoEvent(KSet kset, map<string,string> &best_genes, double score);
  void StreamOutput();  // print csv event info to stderr
  StrPair GetQueryStrs(KSet kset, string region);
  void RunKSet(KSet kset);

private:
  void PrintPath(StrPair query_strs, string gene, string extra_str="");
  sequences GetSubSeqs(KSet kset, string region);
  map<string,sequences> GetSubSeqs(KSet kset);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  size_t GetInsertLength(vector<string> labels);
  size_t GetErosionLength(string side, vector<string> path_labels, string gene_name);
  double AddWithMinusInfinities(double first, double second);

  string hmm_dir_;  // location of .hmm files
  size_t n_seqs_per_track_;  // 1 for plain hmms, 2 for pair hmms
  map<string, set<string> > only_genes_;
  string algorithm_;
  GermLines gl_;
  vector<RecoEvent> events_;
  map<KSet,double> best_scores_;  // map from kset to best score for that kset (summed over regions)
  map<KSet,double> total_scores_;  // map from kset to total score for that kset (summed over regions)
  map<KSet,map<string,string> > best_genes_;  // map from a kset to its corresponding triplet of best genes
  sequences *seqs_;
  HMMHolder *hmms_;
  map<string, map<StrPair,trellis*> > trellisi_;  // collection of the trellises we've calculated, so we can reuse them. eg: trellisi_["IGHV1-18*01"]["ACGGGTCG"] for single hmms, or trellisi_["IGHV1-18*01"][("ACGGGTCG","ATGGTTAG")] for pair hmms
  map<string, map<StrPair,traceback_path*> > paths_;  // collection of the paths. 
  vector<string>::iterator i_current_region_;  // region and position of the *next* job we will pass out with GetNextJob()
  vector<string>::iterator i_current_gene_;
  int debug_;
  size_t n_best_events_; // print and return this many events
};
#endif
