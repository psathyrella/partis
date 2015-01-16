#include <iostream>
#include <fstream>

#include "model.h"
#include "trellis.h"
#include "tclap/CmdLine.h"

using namespace ham;
using namespace TCLAP;
using namespace std;

// ----------------------------------------------------------------------------------------
// split a colon-separated list in a string into a vector of strings, e.g. "a:b:c" --> {"a", "b", "c"}
// NOTE copied from bcrham.cc
vector<string> split_argument(string argstr) {
  vector<string> arglist;
  while(true) {
    size_t i_next_colon(argstr.find(":"));
    string arg = argstr.substr(0, i_next_colon); // get the next arg in the colon-separated list
    arglist.push_back(arg); // add it to arglist
    argstr = argstr.substr(i_next_colon + 1); // then excise it from argstr
    if(i_next_colon == string::npos)
      break;
  }
  return arglist;
}

// ----------------------------------------------------------------------------------------
int main(int argc, const char *argv[]) {

  // set up command line arguments
  ValueArg<string> hmmfname_arg("f", "hmmfname", "hmm (.yaml) model file", true, "", "string");
  ValueArg<string> seqs_arg("s", "seqs", "colon-separated list of sequences", true, "", "string");
  ValueArg<string> outfile_arg("o", "outfile", "output text file", false, "", "string");
  try {
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmfname_arg);
    cmd.add(seqs_arg);
    cmd.add(outfile_arg);
    cmd.parse(argc, argv);
  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  // read hmm model file
  Model hmm;
  hmm.Parse(hmmfname_arg.getValue());
  vector<string> seqstrs = split_argument(seqs_arg.getValue());
  
  // create sequences from command line
  Sequences seqs;
  for(auto &seqstr : seqstrs)
    seqs.AddSeq(Sequence(hmm.track(0), "seq", seqstr));

  // make the trellis, a wrapper for holding the DP tables and running the algorithms
  trellis trell(&hmm, seqs);
  trell.Viterbi();
  TracebackPath path(&hmm);
  trell.Traceback(path);
  cout << "viterbi path (log prob " << trell.ending_viterbi_log_prob() << "):" << endl;
  for(unsigned iseq=0; iseq<seqs.n_seqs(); ++iseq) {
    cout << "  sequence: ";
    seqs[iseq].Print();
  }

  path.abbreviate();
  cout << "  path:     ";
  cout << path;

  trell.Forward();
  cout << "\nforward log prob: " << trell.ending_forward_log_prob() << endl;

  if(outfile_arg.getValue().length() > 0) {
    ofstream ofs;
    ofs.open(outfile_arg.getValue());
    assert(ofs.is_open());
    ofs << trell.ending_forward_log_prob() << "\t" << path << endl;
    ofs.close();
  }
}
