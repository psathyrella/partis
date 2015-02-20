#include <iostream>
#include <fstream>

#include "model.h"
#include "trellis.h"
#include "text.h"
#include "tclap/CmdLine.h"

using namespace ham;
using namespace TCLAP;
using namespace std;

// ----------------------------------------------------------------------------------------
string read_file(string fname) {
  ifstream ifs(fname);
  string seq;
  string line;
  while(getline(ifs, line)) {
    seq += line;
  }
  return seq;
}

// ----------------------------------------------------------------------------------------
int main(int argc, const char *argv[]) {
  // set up command line arguments
  ValueArg<string> hmmfname_arg("f", "hmmfname", "hmm (.yaml) model file", true, "", "string");
  ValueArg<string> seqs_arg("s", "seqs", "colon-separated list of sequences", true, "", "string");
  try {
    CmdLine cmd("ham -- the fantabulous HMM compiler", ' ', "");
    cmd.add(hmmfname_arg);
    cmd.add(seqs_arg);
    cmd.parse(argc, argv);
  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  // read hmm model file
  Model hmm;
  hmm.Parse(hmmfname_arg.getValue());
  string seqstr(read_file(seqs_arg.getValue()));
  printf("length %.0e\n", (float)seqstr.size());
  Sequence seq(hmm.track(0), "seq", seqstr);
  trellis trell(&hmm, seq);
  trell.Forward();
  cout << "\nforward log prob: " << trell.ending_forward_log_prob() << endl;

}
