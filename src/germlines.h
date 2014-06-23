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
  TermColors():
    head_("\033[95m"),
    blue_("\033[94m"),
    green_("\033[92m"),
    yellow_("\033[93m"),
    red_("\033[91m"),
    end_("\033[0m") {
  }

  string red(string seq) { return red_ + seq + end_; }  // return seq with the proper bash escape strings tacked on left and right so it prints as red
  string yellow(string seq) { return yellow_ + seq + end_; }  // same for yellow
  string redify_mutants(string germline, string seq) {  // return <seq> with mutant bases w.r.t. <germline> escaped to appear red (and 'i', inserts, yellow) in bash terminal
    if (germline.size() != seq.size()) {
      cout << "ERROR seqs not same length in redify_mutants: " << germline << endl
	   << "                                              " << seq << endl;
    }
    string return_str;
    for (size_t inuke=0; inuke<seq.size(); ++inuke) {
      if (seq[inuke] == 'i') {
	return_str += yellow(seq.substr(inuke,1));
      } else if (seq[inuke] == germline[inuke]) {
	return_str += seq[inuke];
      } else {
	return_str += red(seq.substr(inuke,1));
      }
    }
    return return_str;
  }

  string head_;
  string blue_;
  string green_;
  string yellow_;
  string red_;
  string end_;
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
  vector<string> regions_;
  map<string,vector<string> > names_;
  map<string,string> seqs_;
};
#endif
