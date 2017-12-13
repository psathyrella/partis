#ifndef BCRUTILS_H
#define BCRUTILS_H

#include <string>
#include <map>
#include <set>
#include <cassert>
#include <sstream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include "model.h"
#include "text.h"

using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
// A collection of utilities for the b cell clustering project.
// Yeah, I would love to put them in a separate dir, but I can't get the darn thing
// to compile when I do that, so I'm punting on it for the moment.

// class to allow easy sorting of per-gene support vectors
class SupportPair {
public:
  SupportPair(string gene, double logprob) : pr_(gene, logprob) {}
  bool operator < (const SupportPair &rhs) const { return logprob() < rhs.logprob(); }  // return true if rhs is more likely than self
  string gene() const { return pr_.first; }
  double logprob() const { return pr_.second; }
  pair<string, double> pr_;
};

// ----------------------------------------------------------------------------------------
class Insertions {
public:
  Insertions() : v_ {"fv"}, d_ {"vd"}, j_ {"dj", "jf"} {
    insertions_["v"] = v_;
    insertions_["d"] = d_;
    insertions_["j"] = j_;
  }
  vector<string> operator[](string region) { return insertions_[region]; }

private:
  vector<string> v_, d_, j_;
  map<string, vector<string> > insertions_;
};

// ----------------------------------------------------------------------------------------
class TermColors {
public:
  TermColors();
  string Color(string col, string seq);  // a verb! it's a verb!
  string RedifyIfMuted(char germline_nuc, char nuc);
  string GetRegion(string gene);
  string ColorChars(char ch_to_color, string color, string seq);
  string ColorMutants(string color, string seq, string ref_1 = "", vector<string> other_refs = {}, string ambiguous_char = "");
  string ColorGene(string gene);

  map<string, string> codes_;
};

// ----------------------------------------------------------------------------------------
class GermLines {
public:
  GermLines(string gldir, string locus);
  string SanitizeName(string gene_name);
  string GetRegion(string gene);

  string locus_;
  vector<string> regions_;
  string dummy_d_gene;  // e.g. for light chain
  map<string, vector<string> > names_;
  map<string, string> seqs_;
  map<string, int> cyst_positions_, tryp_positions_;
};

// ----------------------------------------------------------------------------------------
class RecoEvent {  // keeps track of recombination event. Initially, just to allow printing. Translation of print_reco_event in utils.py
public:
  RecoEvent();
  map<string, string> genes_;
  map<string, size_t> deletions_;
  map<string, string> insertions_;
  string naive_seq_;
  float score_;
  int cyst_position_, tryp_position_, cdr3_length_;
  map<string, vector<SupportPair> >  per_gene_support_;  // for each region, a sorted list of (gene, logprob) pairs

  bool operator < (const RecoEvent& rhs) const { return (score_ < rhs.score_); }
  void SetGenes(string vgene, string dgene, string jgene) { genes_["v"] = vgene; genes_["d"] = dgene; genes_["j"] = jgene; }
  void SetGene(string region, string gene) { genes_[region] = gene; }
  void SetDeletion(string name, size_t len) { deletions_[name] = len; }
  void SetInsertion(string name, string insertion) { insertions_[name] = insertion; }
  void SetNaiveSeq(GermLines &gl);  // NOTE this probably duplicates some code in Print() below, but I don't want to mess with that code at the moment (doesn't really get used any more)
  void SetScore(double score) { score_ = score; }
  void Clear() { genes_.clear(); deletions_.clear(); insertions_.clear(); }
  void Print(GermLines &germlines, size_t cyst_position = 0, size_t final_tryp_position = 0, bool one_line = false, string extra_indent = "");
};

// ----------------------------------------------------------------------------------------
class KSet {  // pair of k_v,k_d values specifying how to chop up the query sequence into v+insert, d+insert, j []
public:
  KSet(size_t k_v, size_t k_d) : v(k_v), d(k_d) {}
  bool isnull() { return v == 0 && d == 0; }
  bool equals(KSet rhs) { return v == rhs.v && d == rhs.d; }
  bool operator< (const KSet rhs) const  // kind of weird, but it's just for std::map
  {
    if(v == rhs.v)
      return d < rhs.d;
    return v < rhs.v;
  }

  size_t v;
  size_t d;
};

// ----------------------------------------------------------------------------------------
class KBounds {
public:
  KBounds() {}  // need this so we can use stl vectors
  KBounds(KSet kmin, KSet kmax) : vmin(kmin.v), dmin(kmin.d), vmax(kmax.v), dmax(kmax.d) {}
  bool equals(KBounds rhs) { return vmin == rhs.vmin && vmax == rhs.vmax && dmin == rhs.dmin && dmax == rhs.dmax; }
  KBounds LogicalOr(KBounds rhs);  // return the "logical OR" of <kb1> and <kb2>, i.e. the area encompassed by either of 'em
  string stringify() {
    stringstream ss;
    ss << vmin << "-" << vmax << ", " << dmin << "-" << dmax;
    return ss.str();
  }
  size_t vmin, dmin, vmax, dmax;
};
// ----------------------------------------------------------------------------------------
class HMMHolder {
public:
  HMMHolder(string hmm_dir, GermLines &gl, Track *track): hmm_dir_(hmm_dir), gl_(gl), track_(track) {}
  ~HMMHolder();
  Model *Get(string gene);
  Track *track() { return track_; }
  // Rescale, within each hmm, the emission probabilities to reflect <overall_mute_freq> instead of the mute freq which was recorded in the hmm file.
  // If <overall_mute_freq> is -INFINITY, we re-rescale them to what they were originally
  void RescaleOverallMuteFreqs(map<string, set<string> > &only_genes, double overall_mute_freq);  // WOE BETIDE THEE WHO FORGETETH TO RE-RESET THESE
  void UnRescaleOverallMuteFreqs(map<string, set<string> > &only_genes);
  void CacheAll();  // read all available hmms into memory
  string NameString(map<string, set<string> > *only_genes=nullptr, int max_to_print=-1);  // if more than <max_to_print> for any region, only print the number of genes for each region
private:
  string hmm_dir_;
  GermLines &gl_;
  map<string, Model*> hmms_; // map of gene name to hmm pointer
  Track *track_;  // each of the models has a track... but they should all be the same, so just toss one here for easy access
};

// ----------------------------------------------------------------------------------------
class Result {
public:
  Result(KBounds kbounds, string locus) : total_score_(-INFINITY), no_path_(false), locus_(locus), better_kbounds_(kbounds), boundary_error_(false), could_not_expand_(false), finalized_(false) {}
  void PushBackRecoEvent(RecoEvent event) { events_.push_back(event); }
  void Finalize(GermLines &gl, map<string, double> &unsorted_per_gene_support, KSet best_kset, KBounds kbounds);
  RecoEvent &best_event() { assert(finalized_); return best_event_; }
  bool boundary_error() { return boundary_error_; } // is the best kset on boundary of k space?  // TODO boundary error stuff is deprectated (since sw does a much smarter job of choosing kbounds), so it can be removed
  bool could_not_expand() { return could_not_expand_; }
  KBounds better_kbounds() { return better_kbounds_; }
  double total_score() { return total_score_; }
  double total_score_;
  bool no_path_;

private:
  void check_boundaries(KSet best, KBounds kbounds);  // and if you find errors, put expanded bounds in better_[kmin,kmax]_

  string locus_;
  KBounds better_kbounds_;
  bool boundary_error_;
  bool could_not_expand_;
  bool finalized_;

  vector<RecoEvent> events_;  // one reco event for each kset, sorted by score in Finalize()
  RecoEvent best_event_;  // most likely event, among those in events_ (this event has its per_gene_support_ set). Set by Finalize().
};

void StreamHeader(ofstream &ofs, string algorithm);
void StreamErrorput(ofstream &ofs, string algorithm, vector<Sequence> &seqs, string errors);
void StreamErrorput(ofstream &ofs, string algorithm, vector<Sequence*> &pseqs, string errors);
string PerGeneSupportString(vector<SupportPair> &support);
void StreamViterbiOutput(ofstream &ofs, RecoEvent &event, vector<Sequence> &seqs, string errors);
void StreamViterbiOutput(ofstream &ofs, RecoEvent &event, vector<Sequence*> &pseqs, string errors);
void StreamForwardOutput(ofstream &ofs, vector<Sequence> &seqs, double total_score, string errors);
void StreamForwardOutput(ofstream &ofs, vector<Sequence*> &pseqs, double total_score, string errors);

string SeqStr(vector<Sequence*> &pseqs, string delimiter = " ");
string SeqStr(vector<Sequence> &seqs, string delimiter = " ");
string SeqNameStr(vector<Sequence*> &pseqs, string delimiter = " ");
string SeqNameStr(vector<Sequence> &seqs, string delimiter = " ");

bool HasDGene(string locus);

// The star-tree assumption causes a systematic bias towards too-long insertions/deletions (since each mutation in each sequence is viewed as the result of an independent mutation event).
// Since the accuracy of the inferred naive sequence does not suffer from significant inaccuracy as a result of this, though (it's largely just taking the consensus sequence in these situations), we can get a better annotation by rerunning with just the naive sequence as input.
// conclusion after validation: fixes the too-long deletion/insertion thing, but is not, globally more accurate (although, it's only *less* accurate within insertins, and see note in HandleFishyAnnotations(), i.e. maybe with some more work that could be fixed)
bool FishyMultiSeqAnnotation(size_t n_seqs, RecoEvent &event);

vector<Sequence> GetSeqVector(vector<Sequence*> pseqvector);

void runps();
int GetMemVal(string name, string path);  // kB
int GetRss();
int GetMemTot();

}

#endif
