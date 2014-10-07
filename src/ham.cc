#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "jobholder.h"
#include "germlines.h"
#include "tclap/CmdLine.h"
using namespace TCLAP;
using namespace ham;
using namespace std;

// ----------------------------------------------------------------------------------------
// class for reading csv input file
class Args {
public:
  Args(int argc, const char * argv[]);
  string hmmdir() { return hmmdir_arg_.getValue(); }
  string datadir() { return datadir_arg_.getValue(); }
  string infile() { return infile_arg_.getValue(); }
  string outfile() { return outfile_arg_.getValue(); }
  string algorithm() { return algorithm_arg_.getValue(); }
  int debug() { return debug_arg_.getValue(); }
  int n_best_events() { return n_best_events_arg_.getValue(); }
  bool pair() { return pair_arg_.getValue(); }

  // command line arguments
  vector<string> algo_strings_;
  vector<int> debug_ints_;
  ValuesConstraint<string> algo_vals_;
  ValuesConstraint<int> debug_vals_;
  ValueArg<string> hmmdir_arg_, datadir_arg_, infile_arg_, outfile_arg_, algorithm_arg_;
  ValueArg<int> debug_arg_, n_best_events_arg_;
  SwitchArg pair_arg_;

  // arguments read from csv input file
  string all_only_genes_;
  map<string, vector<string> > strings_;
  map<string, vector<int> > integers_;
  set<string> str_headers_, int_headers_;
};

// ----------------------------------------------------------------------------------------
Args::Args(int argc, const char * argv[]):
  algo_strings_{"viterbi", "forward"},
  debug_ints_{0,1,2},
  algo_vals_(algo_strings_),
  debug_vals_(debug_ints_),
  hmmdir_arg_("m", "hmmdir", "directory in which to look for hmm model files", true, "", "string"),
  datadir_arg_("d", "datadir", "directory in which to look for non-sample-specific data (eg human germline seqs)", true, "", "string"),
  infile_arg_("i", "infile", "input (whitespace-separated) file", true, "", "string"),
  outfile_arg_("o", "outfile", "output csv file", true, "", "string"),
  algorithm_arg_("a", "algorithm", "algorithm to run", true, "", &algo_vals_),
  debug_arg_("g", "debug", "debug level", false, 0, &debug_vals_),
  n_best_events_arg_("n", "n_best_events", "number of candidate recombination events to write to file", true, -1, "int"),
  pair_arg_("p","pair","is this a pair hmm?", false),
  str_headers_{"only_genes", "name", "seq", "second_name", "second_seq"},
  int_headers_{"k_v_min", "k_v_max", "k_d_min", "k_d_max"}
{
  try {  
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmdir_arg_);
    cmd.add(datadir_arg_);
    cmd.add(infile_arg_);
    cmd.add(outfile_arg_);
    cmd.add(algorithm_arg_);
    cmd.add(debug_arg_);
    cmd.add(n_best_events_arg_);
    cmd.add(pair_arg_);
    cmd.parse(argc, argv);
  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }

  for (auto &head: str_headers_)
    strings_[head] = vector<string>();
  for (auto &head: int_headers_)
    integers_[head] = vector<int>();
  
  ifstream ifs(infile());
  assert(ifs.is_open());
  string line;
  // get header line
  getline(ifs,line);
  stringstream ss(line);
  vector<string> headers;  // keep track of the file's column order
  while (!ss.eof()) {
    string head;
    ss >> head;
    headers.push_back(head);
  }
  while (getline(ifs,line)) {
    stringstream ss(line);
    string tmpstr;
    int tmpint;
    for (auto &head: headers) {
      if (str_headers_.find(head) != str_headers_.end()) {
	ss >> tmpstr;
	strings_[head].push_back(tmpstr);
	if (head=="only_genes") {
	  if (all_only_genes_.size() == 0)
	    all_only_genes_ = tmpstr;
	  else
	    all_only_genes_ += ":" + tmpstr;
	}
      } else if (int_headers_.find(head) != int_headers_.end()) {
	ss >> tmpint;
	integers_[head].push_back(tmpint);
      } else {
	assert(0);
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<Sequences*> GetSeqs(Args &args, Track *trk) {
  vector<Sequences*> all_seqs;
  for (size_t iseq=0; iseq<args.strings_["seq"].size(); ++iseq) {
    Sequences *seqs = new Sequences;
    Sequence *sq = new Sequence(args.strings_["name"][iseq], args.strings_["seq"][iseq], trk);
    seqs->AddSeq(sq);
    if (trk->n_seqs == 2) {
      assert(args.strings_["second_seq"][iseq].size() > 0);
      assert(args.strings_["second_seq"][iseq] != "x");
      Sequence *second_sq = new Sequence(args.strings_["second_name"][iseq], args.strings_["second_seq"][iseq], trk);
      seqs->AddSeq(second_sq);
    } else {
      assert(args.strings_["second_seq"][iseq] == "x");  // er, not really necessary, I suppose...
    }

    all_seqs.push_back(seqs);
  }
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score);
// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  Args args(argc, argv);
  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  assert(ofs.is_open());
  if (args.algorithm() == "viterbi")
    ofs << "unique_id,second_unique_id,v_gene,d_gene,j_gene,vd_insertion,dj_insertion,v_3p_del,d_5p_del,d_3p_del,j_5p_del,score,seq,second_seq,errors" << endl;
  else
    ofs << "unique_id,second_unique_id,score,errors" << endl;

  // init some infrastructure
  size_t n_seqs_per_track(args.pair() ? 2 : 1);
  vector<string> characters{"A","C","G","T"};
  Track trk("NUKES", n_seqs_per_track, characters);
  vector<Sequences*> seqs(GetSeqs(args, &trk));
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), n_seqs_per_track, gl);

  assert(seqs.size() == args.strings_["name"].size());
  for (size_t is=0; is<seqs.size(); is++) {
    if (args.debug()) cout << "  ---------" << endl;
    if (args.pair()) assert(seqs[is]->n_seqs() == 2);
    KSet kmin(args.integers_["k_v_min"][is], args.integers_["k_d_min"][is]);
    KSet kmax(args.integers_["k_v_max"][is], args.integers_["k_d_max"][is]);
    KBounds kbounds(kmin, kmax);

    JobHolder jh(gl, hmms, args.algorithm(), args.strings_["only_genes"][is]);
    jh.SetDebug(args.debug());
    jh.SetNBestEvents(args.n_best_events());

    Result result(kbounds);
    do {
      result = jh.Run(*seqs[is], kbounds);
      if (result.could_not_expand()) cout << "WARNING " << (*seqs[is])[0].name() << ((*seqs[is]).n_seqs()==2 ? (*seqs[is])[1].name() : "") << " couldn't expand k bounds for " << kbounds.stringify() << endl;
      kbounds = result.better_kbounds();
    } while (result.boundary_error() && !result.could_not_expand());

    double score(result.total_score());
    if (args.algorithm() == "forward" && args.pair()) {
      Result result_a = jh.Run((*seqs[is])[0], kbounds);  // denominator in P(A,B) / (P(A) P(B))
      Result result_b = jh.Run((*seqs[is])[1], kbounds);
      
      if (result_a.boundary_error() || result_b.boundary_error()) cout << "WARNING boundary errors for " << (*seqs[is])[0].name() << " " << (*seqs[is])[1].name() << endl;
      if (args.debug()) printf("%70s %8.2f - %8.2f - %8.2f = %8.3f\n", "", score, result_a.total_score(), result_b.total_score(), score - result_a.total_score() - result_b.total_score());
      if (args.debug()) printf("%70s %8.1e / %8.1e / %8.1e = %8.1f\n", "", exp(score), exp(result_a.total_score()), exp(result_b.total_score()), exp(score) / (exp(result_a.total_score())*exp(result_b.total_score())));
      score = score - result_a.total_score() - result_b.total_score();
    }

    if (args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {  // if we were asked for more events than we found
      if (result.events_.size() > 0)
	cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
	assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args, result.events_, *seqs[is], score);
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score) {
  if (args.algorithm() == "viterbi") {
    // assert(args.n_best_events() <= events.size());  // make sure we found at least as many valid events as were asked for with <--n_best_events>
    size_t n_max = min(size_t(args.n_best_events()), events.size());
    for (size_t ievt=0; ievt<n_max; ++ievt) {
      RecoEvent *event = &events[ievt];
      string second_seq_name,second_seq;
      if (args.pair()) {
	second_seq_name = event->second_seq_name_;
	second_seq = event->second_seq_;
      }
      ofs
	<< event->seq_name_
	<< "," << second_seq_name
	<< "," << event->genes_["v"]
	<< "," << event->genes_["d"]
	<< "," << event->genes_["j"]
	<< "," << event->insertions_["vd"]
	<< "," << event->insertions_["dj"]
	<< "," << event->deletions_["v_3p"]
	<< "," << event->deletions_["d_5p"]
	<< "," << event->deletions_["d_3p"]
	<< "," << event->deletions_["j_5p"]
	<< "," << event->score_
	<< "," << event->seq_
	<< "," << second_seq
	<< "," << ""  //errors
	<< endl;
    }
  } else {
    assert(seqs.n_seqs() == 2);  // er, at least for the moment
    ofs
      << seqs[0].name()
      << "," << seqs[1].name()
      << "," << total_score
      << "," << ""  //errors
      << endl;
  }
}
