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
  HMMHolder(string hmm_dir, size_t n_seqs_per_track, GermLines &gl):hmm_dir_(hmm_dir),n_seqs_per_track_(n_seqs_per_track),gl_(gl) {}
  ~HMMHolder();
  model *Get(string gene);
private:
  string hmm_dir_;
  size_t n_seqs_per_track_;
  GermLines &gl_;  // TODO kind of hackey to have a separate one of these in HMMHolder. Then again, I don't think it's really that expensive.
  map<string,model*> hmms_;  // map of gene name to hmm pointer
};

// ----------------------------------------------------------------------------------------
class Result {
public:
  Result() : boundary_error_(""),total_score_(-INFINITY) {}
  string boundary_error_;
  double total_score_;
  vector<RecoEvent> events_;
};

// ----------------------------------------------------------------------------------------
class JobHolder {
public:
  JobHolder(GermLines &gl, HMMHolder &hmms, string algorithm, string only_gene_str="");
  ~JobHolder();
  void Clear();
  Result Run(sequences &seqs, size_t k_v_min, size_t k_v_max, size_t k_d_min, size_t k_d_max);
  Result Run(sequence &seq, size_t k_v_min, size_t k_v_max, size_t k_d_min, size_t k_d_max);
  void RunKSet(sequences &seqs, KSet kset, map<KSet,double> *best_scores, map<KSet,double> *total_scores, map<KSet,map<string,string> > *best_genes);
  void SetDebug(int debug) { debug_ = debug; };
  void SetNBestEvents(size_t n_best) { n_best_events_ = n_best; }
  void FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score);
  void PushBackRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score, vector<RecoEvent> *events);
  RecoEvent FillRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score);
  void StreamOutput(double test);  // print csv event info to stderr
  StrPair GetQueryStrs(sequences &seqs, KSet kset, string region);
  void WriteBestGeneProbs(ofstream &ofs, string query_name);
  string errors() { return errors_; }

private:
  void PrintPath(StrPair query_strs, string gene, double score, string extra_str="");
  sequences GetSubSeqs(sequences &seqs, KSet kset, string region);
  map<string,sequences> GetSubSeqs(sequences &seqs, KSet kset);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  size_t GetInsertLength(vector<string> labels);
  size_t GetErosionLength(string side, vector<string> path_labels, string gene_name);
  double AddWithMinusInfinities(double first, double second);

  string hmm_dir_;  // location of .hmm files
  string errors_;
  GermLines &gl_;
  HMMHolder &hmms_;
  string algorithm_;
  int debug_;
  size_t n_best_events_; // print and return this many events
  map<string, set<string> > only_genes_;

  map<string, map<StrPair,trellis*> > trellisi_;  // collection of the trellises we've calculated, so we can reuse them. eg: trellisi_["IGHV1-18*01"]["ACGGGTCG"] for single hmms, or trellisi_["IGHV1-18*01"][("ACGGGTCG","ATGGTTAG")] for pair hmms
  map<string, map<StrPair,traceback_path*> > paths_;  // collection of the paths. 
  map<string, map<StrPair,double> > all_scores_;
  map<string,double> best_per_gene_scores_;
};
#endif
