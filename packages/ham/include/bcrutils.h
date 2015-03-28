#ifndef BCRUTILS_H
#define BCRUTILS_H

#include <string>
#include <map>
#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
namespace ham {
// ----------------------------------------------------------------------------------------
// A collection of utilities for the b cell clustering project.
// Yeah, I would love to put them in a separate dir, but I can't get the darn thing
// to compile when I do that, so I'm punting on it for the moment.

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
  string ColorMutants(string color, string seq, string ref_1 = "", vector<string> other_refs = {});
  string ColorGene(string gene);

  map<string, string> codes_;
};

// ----------------------------------------------------------------------------------------
class GermLines {
public:
  GermLines(string input_dir);
  string SanitizeName(string gene_name);
  string GetRegion(string gene);

  vector<string> regions_;
  map<string, vector<string> > names_;
  map<string, string> seqs_;
};

// ----------------------------------------------------------------------------------------
class RecoEvent {  // class to keep track of recombination event. Initially, just to allow printing. Translation of print_reco_event in utils.py
public:
  RecoEvent();
  map<string, string> genes_;
  map<string, size_t> deletions_;
  map<string, string> insertions_;
  string seq_name_;
  string seq_;
  string naive_seq_;
  vector<string> auxiliary_seq_names_;
  vector<string> auxiliary_seqs_;
  float score_;

  bool operator < (const RecoEvent& rhs) const { return (score_ < rhs.score_); }  // return true if this event is more likely than the rhs event
  void SetGenes(string vgene, string dgene, string jgene) { genes_["v"] = vgene; genes_["d"] = dgene; genes_["j"] = jgene; }
  void SetGene(string region, string gene) { genes_[region] = gene; }
  void SetDeletion(string name, size_t len) { deletions_[name] = len; }
  void SetInsertion(string name, string insertion) { insertions_[name] = insertion; }
  void SetSeq(string seq_name, string seq) { seq_name_ = seq_name; seq_ = seq; }
  void SetNaiveSeq(GermLines &gl);  // NOTE this probably duplicates some code in Print() below, but I don't want to mess with that code at the moment (doesn't really get used any more)
  void AddAuxiliarySeqs(string name, string seq);  // NOTE this class should in general be treated as representing a *single* event with a *single* sequence. It's just that we allow the possiblity here of attaching auxiliary sequences, but e.g. the insertions should *not* be assumed to correspond to these other sequences
  void SetScore(double score) { score_ = score; }
  void Clear() { genes_.clear(); deletions_.clear(); insertions_.clear(); }
  void Print(GermLines &germlines, size_t cyst_position = 0, size_t final_tryp_position = 0, bool one_line = false, string extra_indent = "");
};
}
#endif
