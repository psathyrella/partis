#include <iostream>
#include <fstream>

#include "model.h"
#include "trellis.h"
#include "tclap/CmdLine.h"

using namespace ham;
using namespace TCLAP;
using namespace std;

int main(int argc, const char *argv[]) {
  // set up command line arguments
  ValueArg<string> hmmfname_arg("f", "hmmfname", "hmm (.yaml) model file", true, "", "string");
  ValueArg<string> seq_arg("s", "seq", "sequence of symbols", true, "", "string");
  ValueArg<string> seq2_arg("t", "seq2", "second sequence", false, "", "string");
  ValueArg<string> outfile_arg("o", "outfile", "output text file", false, "", "string");
  SwitchArg pair_arg("p", "pair", "is this a pair hmm?", false);
  try {
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmfname_arg);
    cmd.add(seq_arg);
    cmd.add(seq2_arg);
    cmd.add(outfile_arg);
    cmd.add(pair_arg);
    cmd.parse(argc, argv);
    if(pair_arg.getValue())
      assert(seq2_arg.getValue() != "");
    if(seq2_arg.getValue() != "")
      assert(pair_arg.getValue());
  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  // read hmm model file
  Model hmm;
  hmm.Parse(hmmfname_arg.getValue());

  // create sequences from command line
  Sequences seqs;
  Sequence seq(hmm.track(0), "seq", seq_arg.getValue());
  seqs.AddSeq(seq);
  if(pair_arg.getValue()) {
    Sequence seq2(hmm.track(0), "seq2", seq2_arg.getValue());
    seqs.AddSeq(seq2);
  }

  // make trellis and run
  trellis trell(&hmm, seqs);
  trell.Viterbi();
  TracebackPath path(&hmm);
  trell.Traceback(path);
  cout << "viterbi path (log prob " << trell.ending_viterbi_log_prob() << "):" << endl;
  cout << "  sequence: ";
  seq.Print();
  if(pair_arg.getValue()) {
    cout << "         2: ";
    seqs[1].Print();
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
