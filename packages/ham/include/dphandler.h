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
  DPHandler(string algorithm, Args *args, GermLines &gl, HMMHolder &hmms);
  ~DPHandler();
  void Clear();
  Result Run(vector<Sequence*> pseqvector, KBounds kbounds, vector<string> only_gene_list = {}, double overall_mute_freq = -INFINITY, bool clear_cache = true);  // run all over the kspace specified by bounds in kmin and kmax
  Result Run(vector<Sequence> seqvector, KBounds kbounds, vector<string> only_gene_list = {}, double overall_mute_freq = -INFINITY, bool clear_cache = true);
  Result Run(Sequence seq, KBounds kbounds, vector<string> only_gene_list = {}, double overall_mute_freq = -INFINITY, bool clear_cache = true);
  void HandleFishyAnnotations(Result &multi_seq_result, vector<Sequence*> pqry_seqs, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq);
  void HandleFishyAnnotations(Result &multi_seq_result, vector<Sequence> qry_seqs, KBounds kbounds, vector<string> only_gene_list, double overall_mute_freq);
  // void StreamOutput(double test);  // print csv event info to stderr
  // void WriteBestGeneProbs(ofstream &ofs, string query_name);
  void PrintCachedTrellisSize();

private:
  void RunKSet(Sequences &seqs, KSet kset, map<string, set<string> > &only_genes, map<KSet, double> *best_scores, map<KSet, double> *total_scores, map<KSet, map<string, string> > *best_genes);
  KSet FindPartialCacheMatch(string region, string gene, KSet kset);
  void InitCache(string gene);
  void FillTrellis(KSet kset, Sequences query_seqs, vector<string> query_strs, string gene, string &origin);
  RecoEvent FillRecoEvent(Sequences &seqs, KSet kset, map<string, string> &best_genes, double score);
  vector<string> GetQueryStrs(Sequences &seqs, KSet kset, string region);

  void PrintPath(KSet kset, vector<string> query_strs, string gene, double score, string extra_str = "");
  Sequences GetSubSeqs(Sequences &seqs, KSet kset, string region);
  map<string, Sequences> GetSubSeqs(Sequences &seqs, KSet kset);  // get the subsequences for the v, d, and j regions given a k_v and k_d
  void SetInsertions(string region, vector<string> path_names, RecoEvent *event);
  size_t GetInsertStart(string side, size_t path_length, size_t insert_length);
  string GetInsertion(string side, vector<string> names);
  size_t GetErosionLength(string side, vector<string> names, string gene_name);

  string algorithm_;
  Args *args_;
  GermLines &gl_;
  HMMHolder &hmms_;

  // NOTE BEWARE DRAGONS AND ALL THAT SHIT!
  // if you add something new here you *must* clear it in Clear(), because we reuse the dphandler for different sequences UPDATE kind of don't do that any more
  // NOTE also that the vector<string> key can take up a ton of memory for multi-hmms with large k UPDATE dammit, no, I don't think that's where the memory was going
  map<string, map<vector<string>, Trellis> > scratch_cachefo_;  // collection of the trellises that  we've calculated from scratch, so we can reuse them. eg: scratch_cachefo_["IGHV1-18*01"]["ACGGGTCG"] for single hmms, or scratch_cachefo_["IGHV1-18*01"][("ACGGGTCG","ATGGTTAG")] for pair hmms
  map<string, map<KSet, TracebackPath> > paths_;
  map<string, map<KSet, double> > scores_;
  map<string, double> per_gene_support_;  // log prob of the best (full) annotation for each gene
};
}
#endif
