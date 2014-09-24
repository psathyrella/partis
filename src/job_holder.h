#ifndef JOB_HOLDER_H
#define JOB_HOLDER_H

#include <map>
#include <string>
#include <sstream>
#include <math.h>
#include <set>
#include <iomanip>
#include "StochHMMlib.h"
#include "StochHMM_usage.h"
#include "germlines.h"

using namespace std;
using namespace StochHMM;

typedef pair<string,string> StrPair;

// ----------------------------------------------------------------------------------------
class KSet {  // pair of k_v,k_d values specifying how to chop up the query sequence into v+insert, d+insert, j
public:
  KSet(size_t k_v, size_t k_d) : v(k_v), d(k_d) {}
  bool equals(KSet rhs) { return v == rhs.v && d == rhs.d; }
  bool operator< (const KSet rhs) const { return  v < rhs.v && d < rhs.d; }

  size_t v;
  size_t d;
};

// ----------------------------------------------------------------------------------------
class KBounds {
public:
  KBounds(KSet kmin, KSet kmax) : vmin(kmin.v), dmin(kmin.d), vmax(kmax.v), dmax(kmax.d) {}
  bool equals(KBounds rhs) { return vmin==rhs.vmin && vmax==rhs.vmax && dmin==rhs.dmin && dmax==rhs.dmax; }
  string stringify() {
    stringstream ss;
    ss << vmin << "-" << vmax << ", " << dmin << "-" << dmax;
    return ss.str();
  }
  size_t vmin,dmin,vmax,dmax;
};
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
  Result(KBounds kbounds) : total_score_(-INFINITY), no_path_(false), better_kbounds_(kbounds), boundary_error_(false), could_not_expand_(false) {}
  void check_boundaries(KSet best, KBounds kbounds);  // and if you find errors, put expanded bounds in better_[kmin,kmax]_
  bool boundary_error() { return boundary_error_; } // is the best kset on boundary of k space?
  bool could_not_expand() { return could_not_expand_; }
  KBounds better_kbounds() { return better_kbounds_; }
  double total_score() { return total_score_; }
  double total_score_;  // TODO move this to private
  bool no_path_;
  vector<RecoEvent> events_;  // TODO move this to private

private:
  KBounds better_kbounds_;
  bool boundary_error_;
  bool could_not_expand_;
};

// ----------------------------------------------------------------------------------------
class JobHolder {
public:
  JobHolder(GermLines &gl, HMMHolder &hmms, string algorithm, string only_gene_str="");
  ~JobHolder();
  void Clear();
  Result Run(sequences &seqs, KBounds kbounds);  // run all over the kspace specified by bounds in kmin and kmax
  Result Run(sequence &seq, KBounds kbounds);
  void RunKSet(sequences &seqs, KSet kset, map<KSet,double> *best_scores, map<KSet,double> *total_scores, map<KSet,map<string,string> > *best_genes);
  void SetDebug(int debug) { debug_ = debug; };
  void SetNBestEvents(size_t n_best) { n_best_events_ = n_best; }
  void FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score);
  void PushBackRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score, vector<RecoEvent> *events);
  RecoEvent FillRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score);
  void StreamOutput(double test);  // print csv event info to stderr
  StrPair GetQueryStrs(sequences &seqs, KSet kset, string region);
  void WriteBestGeneProbs(ofstream &ofs, string query_name);

private:
  void PrintPath(StrPair query_strs, string gene, double score, string extra_str="");
  sequences GetSubSeqs(sequences &seqs, KSet kset, string region);
  map<string,sequences> GetSubSeqs(sequences &seqs, KSet kset);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  size_t GetInsertLength(vector<string> labels);
  size_t GetErosionLength(string side, vector<string> path_labels, string gene_name);
  double AddWithMinusInfinities(double first, double second);

  string hmm_dir_;  // location of .hmm files
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
