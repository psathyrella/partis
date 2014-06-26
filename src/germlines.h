#ifndef GERMLINES_H
#define GERMLINES_H

#include <string>
#include <map>
#include <cassert>
#include <vector>

using namespace std;

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
  string RedifyIfMuted(char germline_nuke, char nuke) {
    if (nuke==germline_nuke)
      return string(1,nuke);
    else
      return Color("red", string(1,nuke));
  }
  // ----------------------------------------------------------------------------------------
  string ColorMutants(string color, string ref_1, string seq, string ref_2="") {  // return <seq> with mutant bases w.r.t. <ref_1> escaped to appear red (and 'i', inserts, yellow) in bash terminal
    if (ref_1.size() == 0)
      return seq;
    if (ref_1.size() != seq.size()) {
      cout << "ERROR seqs not same length in color_mutants: " << ref_1 << endl
	   << "                                              " << seq << endl;
    }
    string return_str;
    for (size_t inuke=0; inuke<seq.size(); ++inuke) {
      if (seq[inuke] == 'i') {
	return_str += Color("yellow", seq.substr(inuke,1));
      } else if (seq[inuke] == ref_1[inuke]) {  // nuke same as ref 1 
	if (ref_2.size()==0 || seq[inuke]==ref_2[inuke])
	  return_str += seq[inuke];
	else
	  return_str += Color(color, seq.substr(inuke,1));
      } else {  // nuke different to ref 1
	if (ref_2.size()==0 || seq[inuke]==ref_2[inuke])
	  return_str += Color(color, seq.substr(inuke,1));
	else
	  return_str += Color("reverse", Color(color, seq.substr(inuke,1)));
      }
    }
    return return_str;
  }

  map<string,string> codes_;
};

// ----------------------------------------------------------------------------------------
class GermLines {
 public:
  // ----------------------------------------------------------------------------------------
  GermLines(string input_dir="/home/dralph/Dropbox/work/recombinator/data"):
    regions_{"v","d","j"}
    {
      for (auto &region : regions_) {
	names_[region] = vector<string>();
	ifstream ifs(input_dir + "/igh" + region + ".fasta");
	assert(ifs.is_open());
	string line,name,seq;
	while (getline(ifs,line)) {
	  if (line[0] == '>') {  // read header lines
	    assert(line[1] != ' ');   // make my life hard, will you?
	    name = line.substr(1, line.find(" ") - 1);  // skip the '>', and run until the first blank. It *should* be the gene name. We'll find out later when we look for the file.
	    names_[region].push_back(name);
	  } else {
	    // line.replace(line.find("\n"), 1, "");
	    seq = (seq=="") ? line : seq + line;
	    if (ifs.peek()=='>' || ifs.peek()==EOF) {  // if we're done with this sequence, add it to the map
	      seqs_[name] = seq;
	      seq = "";
	    }
	  }
	}
	ifs.close();
      }
    }

  // ----------------------------------------------------------------------------------------
  // replace * with _star_ and / with _slash_
  string SanitizeName(string gene_name) {
    size_t istar(gene_name.find("*"));
    if (istar != string::npos)
      gene_name.replace(istar, 1, "_star_");

    size_t islash(gene_name.find("/"));
    if (islash != string::npos)
      gene_name.replace(islash, 1, "_slash_");

    return gene_name;
  }

  // ----------------------------------------------------------------------------------------
  // get region from gene name, e.g. IGHD1-1*01 --> d
  string GetRegion(string gene) {
    assert(gene.find("IGH") != string::npos);
    char region_char(tolower(gene[3]));
    string region(1, region_char);
    assert(region=="v" || region=="d" || region=="j");
    return region;
  }

  // ----------------------------------------------------------------------------------------
  vector<string> regions_;
  map<string,vector<string> > names_;
  map<string,string> seqs_;
};

// ----------------------------------------------------------------------------------------
class RecoEvent {  // class to keep track of recombination event. Initially, just to allow printing. Translation of print_reco_event in utils.py
public:
  map<string,string> genes_;
  map<string,size_t> deletions_;
  map<string,string> insertions_;
  string seq_;

  void SetGenes(string vgene, string dgene, string jgene) { genes_["v"] = vgene; genes_["d"] = dgene; genes_["j"] = jgene; }
  void SetGene(string region, string gene) { genes_[region] = gene; }
  void SetDeletion(string name, size_t len) { deletions_[name] = len; }
  void SetInsertion(string name, string insertion) { insertions_[name] = insertion; }
  void SetSeq(string seq) { seq_ = seq; }
  void Clear() { genes_.clear(); deletions_.clear(); insertions_.clear(); }
  void Print(GermLines &germlines, size_t cyst_position=0, size_t final_tryp_position=0, bool one_line=false) {
    // make sure we remembered to fill all the necessary information
    vector<string> deletion_names{"v_3p", "d_5p", "d_3p", "j_5p"};
    vector<string> insertion_names{"vd", "dj"};
    for (auto &region: germlines.regions_) assert(genes_.find(region) != genes_.end());
    for (auto &name: deletion_names) assert(deletions_.find(name) != deletions_.end());
    for (auto &name: insertion_names) assert(insertions_.find(name) != insertions_.end());
    assert(seq_.size() > 0);

    map<string,string> original_seqs;
    for (auto &region: germlines.regions_)
      original_seqs[region] = germlines.seqs_[genes_[region]];
    size_t v_length = original_seqs["v"].size() - deletions_["v_3p"];
    size_t d_length = original_seqs["d"].size() - deletions_["d_5p"] - deletions_["d_3p"];
    size_t j_length = original_seqs["j"].size() - deletions_["j_5p"];
    map<string,string> eroded_seqs;
    eroded_seqs["v"] = original_seqs["v"].substr(0, original_seqs["v"].size() - deletions_["v_3p"]);
    eroded_seqs["d"] = original_seqs["d"].substr(deletions_["d_5p"], original_seqs["d"].size() - deletions_["d_3p"] - deletions_["d_5p"]);
    eroded_seqs["j"] = original_seqs["j"].substr(deletions_["j_5p"]);

    size_t germline_v_end = original_seqs["v"].size() - 1;
    size_t germline_d_start = original_seqs["v"].size() - deletions_["v_3p"] + insertions_["vd"].size() - deletions_["d_5p"];
    size_t germline_d_end = germline_d_start + original_seqs["d"].size();
    size_t germline_j_start = germline_d_end + 1 - deletions_["d_3p"] + insertions_["dj"].size() - deletions_["j_5p"];

    // see if we've "eroded" left side of v or right side of j, and if so add some more dots
    if (deletions_.find("v_5p") != deletions_.end()) {
      for (size_t i=0; i<deletions_["v_5p"]; ++i)
	seq_ = "." + seq_;
    }
    if (deletions_.find("j_3p") != deletions_.end()) {
      for (size_t i=0; i<deletions_["j_3p"]; ++i)
	seq_ += ".";
    }
    
    string final_seq;
    TermColors tc;
    for (size_t inuke=0; inuke<seq_.size(); ++inuke) {
      size_t ilocal = inuke;
      string new_nuke;
      if (ilocal < v_length) {
	new_nuke = tc.RedifyIfMuted(eroded_seqs["v"][ilocal], seq_[inuke]);
      } else {
	ilocal -= v_length;
	if (ilocal < insertions_["vd"].size()) {
	  new_nuke = tc.RedifyIfMuted(insertions_["vd"][ilocal], seq_[inuke]);
        } else {
	  ilocal -= insertions_["vd"].size();
	  if (ilocal < d_length) {
	    new_nuke = tc.RedifyIfMuted(eroded_seqs["d"][ilocal], seq_[inuke]);
	  } else {
	    ilocal -= d_length;
	    if (ilocal < insertions_["dj"].size()) {
	      new_nuke = tc.RedifyIfMuted(insertions_["dj"][ilocal], seq_[inuke]);
	    } else {
	      ilocal -= insertions_["dj"].size();
	      new_nuke = tc.RedifyIfMuted(eroded_seqs["j"][ilocal], seq_[inuke]);
	    }
	  }
	}
      }
      // reverse video the conserved codons, if we have them
      if (cyst_position > 0 && final_tryp_position > 0) {
	if (inuke == cyst_position || inuke == final_tryp_position) {
	  new_nuke = "\033[7m" + new_nuke;
	} else if (inuke == cyst_position + 2 || inuke == final_tryp_position + 2) {
	  new_nuke = new_nuke + "\033[m";
	}
      }
      // and finally tack it onto the final sequence
      final_seq += new_nuke;
    }

    // pad with dots
    for (size_t i=0; i<deletions_["v_3p"]; ++i)
      eroded_seqs["v"] += ".";
    for (size_t i=0; i<deletions_["d_5p"]; ++i)
      eroded_seqs["d"] = "." + eroded_seqs["d"];
    for (size_t i=0; i<deletions_["d_3p"]; ++i)
      eroded_seqs["d"] += ".";
    for (size_t i=0; i<deletions_["j_5p"]; ++i)
      eroded_seqs["j"] = "." + eroded_seqs["j"];

    string insertion_str;
    for (size_t iv=0; iv<v_length; ++iv)
      insertion_str += " ";
    insertion_str += insertions_["vd"];
    for (size_t id=0; id<d_length; ++id)
      insertion_str += " ";
    insertion_str += insertions_["dj"];
    for (size_t ij=0; ij<j_length; ++ij)
      insertion_str += " ";
    string d_str;
    for (size_t ig=0; ig<germline_d_start; ++ig)
      d_str += " ";
    d_str += eroded_seqs["d"];
    for (size_t ij=0; ij<(original_seqs["j"].size() - deletions_["j_5p"] + insertions_["dj"].size() - deletions_["d_3p"]); ++ij)
      d_str += " ";
    string vj_str(eroded_seqs["v"]);
    for (size_t ivj=0; ivj<(germline_j_start - germline_v_end - 2); ++ivj)
      vj_str += " ";
    vj_str += eroded_seqs["j"];

    if (!one_line) {
      cout << "    " << insertion_str << "   inserts" << endl;
      cout <<  "    " << d_str << "   ighd" << endl;
      cout << "    " << vj_str << "   ighv,ighj\n" << endl;
    }
    cout << "    " << final_seq << endl;
  }
};
#endif
