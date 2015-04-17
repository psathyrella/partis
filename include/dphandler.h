#ifndef HAM_DPHANDLER_H
#define HAM_DPHANDLER_H

#include <map>
#include <string>
#include <sstream>
#include <math.h>
#include <set>
#include <iomanip>
#include <stdexcept>

#include "trellis.h"
#include "mathutils.h"
#include "bcrutils.h"
#include "args.h"

using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
class DPHandler {
public:
  DPHandler(Args *args, GermLines &gl, HMMHolder &hmms, vector<string> only_genes = {});
  ~DPHandler();
  // void PrintHMMS() { hmms_.Print(); }
  void Clear();
  Result Run(string algorithm, Sequences seqs, KBounds kbounds, double overall_mute_freq = -INFINITY);  // run all over the kspace specified by bounds in kmin and kmax
  Result Run(string algorithm, Sequence seq, KBounds kbounds, double overall_mute_freq = -INFINITY);
  void StreamOutput(double test);  // print csv event info to stderr
  // void WriteBestGeneProbs(ofstream &ofs, string query_name);

private:
  void RunKSet(string algorithm, Sequences &seqs, KSet kset, map<KSet, double> *best_scores, map<KSet, double> *total_scores, map<KSet, map<string, string> > *best_genes);
  void FillTrellis(string algorithm, Sequences query_seqs, vector<string> query_strs, string gene, double *score, string &origin);
  void PushBackRecoEvent(Sequences &seqs, KSet kset, map<string, string> &best_genes, double score, vector<RecoEvent> *events);
  RecoEvent FillRecoEvent(Sequences &seqs, KSet kset, map<string, string> &best_genes, double score);
  vector<string> GetQueryStrs(Sequences &seqs, KSet kset, string region);

  void PrintPath(vector<string> query_strs, string gene, double score, string extra_str = "");
  Sequences GetSubSeqs(Sequences &seqs, KSet kset, string region);
  map<string, Sequences> GetSubSeqs(Sequences &seqs, KSet kset);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  void SetInsertions(string region, string query_str, vector<string> path_names, RecoEvent *event);
  size_t GetInsertStart(string side, size_t path_length, size_t insert_length);
  size_t GetInsertLength(string side, vector<string> names);
  size_t GetErosionLength(string side, vector<string> names, string gene_name);

  Args *args_;
  string hmm_dir_;  // location of .hmm files
  GermLines &gl_;
  HMMHolder &hmms_;
  map<string, set<string> > only_genes_;

  map<string, map<vector<string>, trellis*> > trellisi_; // collection of the trellises we've calculated, so we can reuse them. eg: trellisi_["IGHV1-18*01"]["ACGGGTCG"] for single hmms, or trellisi_["IGHV1-18*01"][("ACGGGTCG","ATGGTTAG")] for pair hmms
  map<string, map<vector<string>, TracebackPath*> > paths_; // collection of the paths.
  map<string, map<vector<string>, double> > all_scores_;
  map<string, double> best_per_gene_scores_;
};
}
#endif
