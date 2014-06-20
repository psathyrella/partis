#ifndef GERMLINES_H
#define GERMLINES_H

#include <string>
#include <map>
#include <cassert>
#include <vector>

using namespace std;

// ----------------------------------------------------------------------------------------
class GermLines {
 public:
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
	    if (ifs.peek() == '>') {  // if we're done with this sequence, add it to the map
	      seqs_[name] = seq;
	      seq = "";
	    }
	  }
	}
	ifs.close();
      }
    }

  vector<string> regions_;
  map<string,vector<string> > names_;
  map<string,string> seqs_;
};

#endif
