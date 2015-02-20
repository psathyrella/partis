#ifndef GERMLINES_H
#define GERMLINES_H

#include <string>
#include <map>
#include <cassert>
#include <vector>

using namespace std;

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
  TermColors() {
    codes_["bold"] = "\033[1m";
    codes_["reverse"] = "\033[7m";
    codes_["purple"] = "\033[95m";
    codes_["blue"] = "\033[94m";
    codes_["green"] = "\033[92m";
    codes_["yellow"] = "\033[93m";
    codes_["red"] = "\033[91m";
    codes_["end"] = "\033[0m";
  }

  // ----------------------------------------------------------------------------------------
  string Color(string col, string seq) {  // a verb! it's a verb!
    assert(codes_.find(col) != codes_.end());
    return codes_[col] + seq + codes_["end"];
  }
  // ----------------------------------------------------------------------------------------
  string RedifyIfMuted(char germline_nuc, char nuc) {
    if(nuc == germline_nuc)
      return string(1, nuc);
    else
      return Color("red", string(1, nuc));
  }

  // ----------------------------------------------------------------------------------------
  // get region from gene name, e.g. IGHD1-1*01 --> d
  string GetRegion(string gene) {
    assert(gene.find("IGH") != string::npos);
    char region_char(tolower(gene[3]));
    string region(1, region_char);
    assert(region == "v" || region == "d" || region == "j");
    return region;
  }
  // ----------------------------------------------------------------------------------------
  string ColorMutants(string color, string seq, string ref_1 = "", vector<string> other_refs = {}) {
    // Return <seq> with mutant bases w.r.t. <ref_1> escaped to appear red (and 'i', inserts, yellow) in bash terminal.
    // If <refs> are specified, use bold text and reverse video to show if <seq> is muted w/respect to more than one ref
    if(ref_1 != "")
      other_refs.push_back(ref_1);  // only doing it this way so we can call it without specifying <other_refs>
    for(auto & ref : other_refs) {
      if(ref.size() != seq.size()) {
        throw runtime_error("ERROR seqs not same length in color_mutants: " + ref + "\n                                              " + seq);
      }
    }

    string return_str;
    for(size_t inuc = 0; inuc < seq.size(); ++inuc) {
      if(seq[inuc] == 'i') {
        return_str += Color("yellow", seq.substr(inuc, 1));
      } else {
        int ndiff(0);  // number of reference sequences that differ at base inuc
        for(auto & ref : other_refs) {
          if(seq[inuc] != ref[inuc])
            ndiff += 1;
        }
        if(ndiff == 0)
          return_str += seq[inuc];
        else if(ndiff == 1)
          return_str += Color(color, seq.substr(inuc, 1));
        else
          return_str += Color("reverse", Color(color, seq.substr(inuc, 1)));
      }
    }
    return return_str;
  }
  // ----------------------------------------------------------------------------------------
  string ColorGene(string gene) {
    string return_str = gene.substr(0, 3) + Color("purple", gene.substr(3, 1));
    string n_version = gene.substr(4, gene.find("-") - 4);
    string n_subversion = gene.substr(gene.find("-") + 1, gene.find("*") - gene.find("-") - 1);
    if(GetRegion(gene) == "j") {
      n_version = gene.substr(4, gene.find("*") - 4);
      n_subversion = "";
      return_str += Color("red", n_version);
    } else {
      return_str += Color("red", n_version) + "-" + Color("red", n_subversion);
    }
    string allele = gene.substr(gene.find("*") + 1, gene.find("_") - gene.find("*") - 1);
    return_str += "*" + Color("yellow", allele);
    if(gene.find("_") != string::npos)   // _F or _P in j gene names
      return_str += gene.substr(gene.find("_"));
    return return_str;
  }

  // ----------------------------------------------------------------------------------------
  map<string, string> codes_;
};

// ----------------------------------------------------------------------------------------
class GermLines {
public:
  // ----------------------------------------------------------------------------------------
  GermLines(string input_dir):
    regions_ {"v", "d", "j"} {
    for(auto & region : regions_) {
      names_[region] = vector<string>();
      ifstream ifs(input_dir + "/igh" + region + ".fasta");
      assert(ifs.is_open());
      string line, name, seq;
      while(getline(ifs, line)) {
        if(line[0] == '>') {   // read header lines
          assert(line[1] != ' ');   // make my life hard, will you?
          name = line.substr(1, line.find(" ") - 1);  // skip the '>', and run until the first blank. It *should* be the gene name. We'll find out later when we look for the file.
          names_[region].push_back(name);
        } else {
          // line.replace(line.find("\n"), 1, "");
          seq = (seq == "") ? line : seq + line;
          if(ifs.peek() == '>' || ifs.peek() == EOF) { // if we're done with this sequence, add it to the map
            seqs_[name] = seq;
            seq = "";
          }
        }
      }
      ifs.close();
    }
  }

  // ----------------------------------------------------------------------------------------
  // replace * with _star_ and /OR15 with _slash_
  string SanitizeName(string gene_name) {
    size_t istar(gene_name.find("*"));
    if(istar != string::npos)
      gene_name.replace(istar, 1, "_star_");

    size_t islash(gene_name.find("/"));
    if(islash != string::npos)
      gene_name.replace(islash, 1, "_slash_");

    return gene_name;
  }

  // ----------------------------------------------------------------------------------------
  // get region from gene name, e.g. IGHD1-1*01 --> d
  string GetRegion(string gene) {
    assert(gene.find("IGH") != string::npos);
    char region_char(tolower(gene[3]));
    string region(1, region_char);
    assert(region == "v" || region == "d" || region == "j");
    return region;
  }

  // ----------------------------------------------------------------------------------------
  vector<string> regions_;
  map<string, vector<string> > names_;
  map<string, string> seqs_;
};

// ----------------------------------------------------------------------------------------
class RecoEvent {  // class to keep track of recombination event. Initially, just to allow printing. Translation of print_reco_event in utils.py
public:
  RecoEvent() { score_ = 999; }
  map<string, string> genes_;
  map<string, size_t> deletions_;
  map<string, string> insertions_;
  string seq_name_;
  string seq_;
  vector<string> auxiliary_seq_names_;
  vector<string> auxiliary_seqs_;
  float score_;

  bool operator < (const RecoEvent& rhs) const { return (score_ < rhs.score_); }  // return true if this event is more likely than the rhs event
  void SetGenes(string vgene, string dgene, string jgene) { genes_["v"] = vgene; genes_["d"] = dgene; genes_["j"] = jgene; }
  void SetGene(string region, string gene) { genes_[region] = gene; }
  void SetDeletion(string name, size_t len) { deletions_[name] = len; }
  void SetInsertion(string name, string insertion) { insertions_[name] = insertion; }
  void SetSeq(string seq_name, string seq) { seq_name_ = seq_name; seq_ = seq; }
  void AddAuxiliarySeqs(string name, string seq) {  // NOTE this class should in general be treated as representing a *single* event with a *single* sequence. It's just that we allow the possiblity here of attaching auxiliary sequences, but e.g. the insertions should *not* be assumed to correspond to these other sequences
    assert(seq.size() == seq_.size());  // make sure we already have a first seq, and also that the new seq is the same size
    auxiliary_seq_names_.push_back(name);
    auxiliary_seqs_.push_back(seq);
  }

  void SetScore(double score) { score_ = score; }
  void Clear() { genes_.clear(); deletions_.clear(); insertions_.clear(); }
  void Print(GermLines &germlines, size_t cyst_position = 0, size_t final_tryp_position = 0, bool one_line = false, string extra_indent = "") {
    assert(0);  // this needs to be updated to match the python version
    string print_seq = seq_;  // this variable is a relic of past transgressions, should be redundantized when I finish what I'm doing
    // need to update this to allow printing of auxiliary sequences... or not. Maybe just do it in the python version. I don't really use this anymore, actually
    // if(print_second_seq)
    //   print_seq = second_seq_;

    if(score_ == -INFINITY) {
      cout << "    " << extra_indent << score_ << endl;
      return;
    }
    // make sure we remembered to fill all the necessary information
    vector<string> deletion_names {"v_3p", "d_5p", "d_3p", "j_5p"};
    vector<string> insertion_names {"vd", "dj"};
    for(auto & region : germlines.regions_) assert(genes_.find(region) != genes_.end());
    for(auto & name : deletion_names) assert(deletions_.find(name) != deletions_.end());
    for(auto & name : insertion_names) assert(insertions_.find(name) != insertions_.end());
    assert(print_seq.size() > 0);

    map<string, string> original_seqs;
    for(auto & region : germlines.regions_)
      original_seqs[region] = germlines.seqs_[genes_[region]];

    if(deletions_.find("v_5p") == deletions_.end()) {  // try to infer the left-hand v "deletion"
      deletions_["v_5p"] = original_seqs["v"].size() + original_seqs["d"].size() + original_seqs["j"].size()
                           - deletions_["v_3p"] - deletions_["d_5p"] - deletions_["d_3p"] - deletions_["j_5p"]
                           + insertions_["vd"].size() + insertions_["dj"].size() - print_seq.size();
    }
    original_seqs["v"] = original_seqs["v"].substr(deletions_["v_5p"]);

    size_t v_length = original_seqs["v"].size() - deletions_["v_3p"];
    size_t d_length = original_seqs["d"].size() - deletions_["d_5p"] - deletions_["d_3p"];
    size_t j_length = original_seqs["j"].size() - deletions_["j_5p"];

    map<string, string> eroded_seqs;
    eroded_seqs["v"] = original_seqs["v"].substr(0, original_seqs["v"].size() - deletions_["v_3p"]);
    eroded_seqs["d"] = original_seqs["d"].substr(deletions_["d_5p"], original_seqs["d"].size() - deletions_["d_3p"] - deletions_["d_5p"]);
    eroded_seqs["j"] = original_seqs["j"].substr(deletions_["j_5p"]);

    size_t germline_v_end = original_seqs["v"].size() - 1;
    size_t germline_d_start = original_seqs["v"].size() - deletions_["v_3p"] + insertions_["vd"].size() - deletions_["d_5p"];
    size_t germline_d_end = germline_d_start + original_seqs["d"].size();
    size_t germline_j_start = germline_d_end + 1 - deletions_["d_3p"] + insertions_["dj"].size() - deletions_["j_5p"];

    // see if we've "eroded" left side of v or right side of j, and if so add some more dots
    // if (deletions_.find("v_5p") != deletions_.end()) {  // truncating this part for the moment
    //   for (size_t i=0; i<deletions_["v_5p"]; ++i)
    // 	print_seq = "." + print_seq;
    // }
    if(deletions_.find("j_3p") != deletions_.end()) {
      for(size_t i = 0; i < deletions_["j_3p"]; ++i)
        print_seq += ".";
    }

    string final_seq;
    TermColors tc;
    for(size_t inuc = 0; inuc < print_seq.size(); ++inuc) {
      size_t ilocal = inuc;
      string new_nuc;
      if(ilocal < v_length) {
        new_nuc = tc.RedifyIfMuted(eroded_seqs["v"][ilocal], print_seq[inuc]);
      } else {
        ilocal -= v_length;
        if(ilocal < insertions_["vd"].size()) {
          new_nuc = tc.RedifyIfMuted(insertions_["vd"][ilocal], print_seq[inuc]);
        } else {
          ilocal -= insertions_["vd"].size();
          if(ilocal < d_length) {
            new_nuc = tc.RedifyIfMuted(eroded_seqs["d"][ilocal], print_seq[inuc]);
          } else {
            ilocal -= d_length;
            if(ilocal < insertions_["dj"].size()) {
              new_nuc = tc.RedifyIfMuted(insertions_["dj"][ilocal], print_seq[inuc]);
            } else {
              ilocal -= insertions_["dj"].size();
              new_nuc = tc.RedifyIfMuted(eroded_seqs["j"][ilocal], print_seq[inuc]);
            }
          }
        }
      }
      // reverse video the conserved codons, if we have them
      if(cyst_position > 0 && final_tryp_position > 0) {
        if(inuc == cyst_position || inuc == final_tryp_position) {
          new_nuc = "\033[7m" + new_nuc;
        } else if(inuc == cyst_position + 2 || inuc == final_tryp_position + 2) {
          new_nuc = new_nuc + "\033[m";
        }
      }
      // and finally tack it onto the final sequence
      final_seq += new_nuc;
    }

    // pad with dots
    for(size_t i = 0; i < deletions_["v_3p"]; ++i)
      eroded_seqs["v"] += ".";
    for(size_t i = 0; i < deletions_["d_5p"]; ++i)
      eroded_seqs["d"] = "." + eroded_seqs["d"];
    for(size_t i = 0; i < deletions_["d_3p"]; ++i)
      eroded_seqs["d"] += ".";
    for(size_t i = 0; i < deletions_["j_5p"]; ++i)
      eroded_seqs["j"] = "." + eroded_seqs["j"];

    string insertion_str;
    for(size_t iv = 0; iv < v_length; ++iv)
      insertion_str += " ";
    insertion_str += insertions_["vd"];
    for(size_t id = 0; id < d_length; ++id)
      insertion_str += " ";
    insertion_str += insertions_["dj"];
    for(size_t ij = 0; ij < j_length; ++ij)
      insertion_str += " ";
    string d_str;
    for(size_t ig = 0; ig < germline_d_start; ++ig)
      d_str += " ";
    d_str += eroded_seqs["d"];
    int ij_max = original_seqs["j"].size() - deletions_["j_5p"] + insertions_["dj"].size() - deletions_["d_3p"];
    // assert(ij_max>=0);
    if(ij_max < 0)
      ij_max = 0;
    for(int ij = 0; ij < ij_max; ++ij)
      d_str += " ";
    string vj_str(eroded_seqs["v"]);
    int ivj_max = germline_j_start - germline_v_end - 2;
    // assert(ivj_max >= 0);
    if(ivj_max < 0)
      ivj_max = 0;
    for(int ivj = 0; ivj < ivj_max; ++ivj)
      vj_str += " ";
    vj_str += eroded_seqs["j"];

    if(!one_line) {
      cout << extra_indent << "    " << insertion_str << "   inserts" << endl;
      cout << extra_indent <<  "    " << d_str << "   " << tc.ColorGene(genes_["d"]) << endl;
      cout << extra_indent << "    " << vj_str << "   " << tc.ColorGene(genes_["v"]) << "," << tc.ColorGene(genes_["j"]) << endl;
    }
    cout << extra_indent << "    " << final_seq << "   " << score_ << endl;
  }
};
#endif
