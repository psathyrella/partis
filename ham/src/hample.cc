#include <iostream>

#include "model.h"
#include "trellis.h"
#include "tclap/CmdLine.h"

using namespace ham;
using namespace TCLAP;
using namespace std;

int main(int argc, const char *argv[]) {
  ValueArg<string> hmmfname_arg("f", "hmmfname", "hmm (.yaml) model file", true, "", "string");
  ValueArg<string> seq_arg("s", "seq", "sequence of dice rolls", true, "", "string");
  try {
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmfname_arg);
    cmd.add(seq_arg);
    cmd.parse(argc, argv);
  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  Model hmm;
  hmm.Parse(hmmfname_arg.getValue());
  // Track track("dice", 1, {"1", "2", "3", "4", "5", "6"});
  Sequence seq("rolls", seq_arg.getValue(), hmm.track(0));
  trellis trell(&hmm, &seq);
  trell.Viterbi();
  TracebackPath path(&hmm);
  trell.Traceback(path);
  cout << "viterbi path (log prob " << trell.ending_viterbi_log_prob() << "):" << endl;
  cout << "  sequence: ";
  seq.Print();
  cout << "  path:     ";
  path.print_labels();
  trell.Forward();
  cout << "\nforward log prob: " << trell.forward_log_prob() << endl;
}
