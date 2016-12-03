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
void CheckChunkCaching(Model &hmm, Trellis &trellis, Sequences seqs);  // for checking with scons test, ignore if you're not scons

// ----------------------------------------------------------------------------------------
int main(int argc, const char *argv[]) {

  // set up command line arguments
  ValueArg<string> hmmfname_arg("f", "hmmfname", "hmm (.yaml) model file", true, "", "string");
  ValueArg<string> seqs_arg("s", "seqs", "colon-separated list of sequences", true, "", "string");
  ValueArg<string> outfile_arg("o", "outfile", "output text file", false, "", "string");
  try {
    CmdLine cmd("ham -- the fantabulous HMM compiler", ' ', "");
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
  vector<string> seqstrs = SplitString(seqs_arg.getValue(), ":");

  // create sequences from command line
  Sequences seqs;
  for(auto & seqstr : seqstrs)
    seqs.AddSeq(Sequence(hmm.track(), "seq", seqstr));

  // make the trellis, a wrapper for holding the DP tables and running the algorithms
  Trellis trell(&hmm, seqs);
  trell.Viterbi();
  TracebackPath path(&hmm);
  trell.Traceback(path);
  cout << "viterbi path (log prob " << trell.ending_viterbi_log_prob() << "):" << endl;
  for(unsigned iseq = 0; iseq < seqs.n_seqs(); ++iseq) {
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
  CheckChunkCaching(hmm, trell, seqs);
}

// ----------------------------------------------------------------------------------------
// check dp table chunk caching (just for use by `scons test`)
void CheckChunkCaching(Model &hmm, Trellis &trell, Sequences seqs) {
  for(size_t length = 1; length < seqs.GetSequenceLength(); ++length) {
    // first make a new trellis on the substring of length <length> using chunk caching
    Sequences subseqs(seqs, 0, length);
    Trellis subtrell(&hmm, subseqs, &trell);
    subtrell.Viterbi();
    TracebackPath subpath(&hmm);
    subtrell.Traceback(subpath);
    subtrell.Forward();

    // then make another one on the same sequence from scratch
    Trellis checktrell(&hmm, subseqs);
    checktrell.Viterbi();
    TracebackPath checkpath(&hmm);
    checktrell.Traceback(checkpath);
    checktrell.Forward();

    if(subtrell.ending_forward_log_prob() == -INFINITY) {
      cout << "  no valid path for subsequence of length " << length << endl;
      continue;
    }

    // and finally, make sure they got the same answer
    assert(checkpath.size() == subpath.size());
    assert(checkpath.size() == length);
    for(size_t ipos = 0; ipos < length; ++ipos) {
      if(checkpath[ipos] != subpath[ipos])
        throw runtime_error("ERROR dp table chunk caching failed -- didn't give the same viterbi path");
    }
    double eps(1e-10);
    if(fabs(checktrell.ending_viterbi_log_prob() - subtrell.ending_viterbi_log_prob()) > eps)
      throw runtime_error("ERROR dp table chunk caching failed -- didn't give the same viterbi log prob " + to_string(checktrell.ending_viterbi_log_prob()) + " " + to_string(subtrell.ending_viterbi_log_prob()));
    if(fabs(checktrell.ending_forward_log_prob() - subtrell.ending_forward_log_prob()) > eps)
      throw runtime_error("ERROR dp table chunk caching failed -- didn't give the same forward log prob: " + to_string(checktrell.ending_forward_log_prob()) + " " + to_string(subtrell.ending_forward_log_prob()));
  }
  cout << "caching ok!" << endl;
}
