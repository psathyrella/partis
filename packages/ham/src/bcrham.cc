#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iomanip>
// #include <ctime>
#include <fstream>
#include <cfenv>
#include "jobholder.h"
#include "germlines.h"
#include "tclap/CmdLine.h"
using namespace TCLAP;
using namespace ham;
using namespace std;

// ----------------------------------------------------------------------------------------
// input processing class
// NOTE some input is passed on the command line (global configuration), while some is passed in a csv file (stuff that depends on each (pair of) sequence(s)).
class Args {
public:
  Args(int argc, const char * argv[]);
  vector<string> Split(string arglist);  // str.split(':'), but all hackey-like 'cause it's c++
  string hmmdir() { return hmmdir_arg_.getValue(); }
  string datadir() { return datadir_arg_.getValue(); }
  string infile() { return infile_arg_.getValue(); }
  string outfile() { return outfile_arg_.getValue(); }
  string algorithm() { return algorithm_arg_.getValue(); }
  // string algorithm() { return algorithm_; }
  int debug() { return debug_arg_.getValue(); }
  // int debug() { return debug_; }
  int n_best_events() { return n_best_events_arg_.getValue(); }
  bool pair() { return pair_arg_.getValue(); }
  bool chunk_cache() { return chunk_cache_arg_.getValue(); }

  // command line arguments
  vector<string> algo_strings_;
  vector<int> debug_ints_;
  ValuesConstraint<string> algo_vals_;
  ValuesConstraint<int> debug_vals_;
  ValueArg<string> hmmdir_arg_, datadir_arg_, infile_arg_, outfile_arg_, algorithm_arg_;
  ValueArg<int> debug_arg_, n_best_events_arg_;
  SwitchArg pair_arg_;
  SwitchArg chunk_cache_arg_;

  // arguments read from csv input file
  map<string, vector<string> > strings_;
  map<string, vector<int> > integers_;
  map<string, vector<vector<string> > > str_lists_;
  set<string> str_headers_, int_headers_, str_list_headers_;

  // // extra values to cache command line args (TCLAP calls to ValuesConstraint::check() seem to be really slow
  // UPDATE hmm, didn't seem to help. leave it for the moment
  // string algorithm_;
  // int debug_;
};

// ----------------------------------------------------------------------------------------
Args::Args(int argc, const char * argv[]):
  algo_strings_ {"viterbi", "forward"},
  debug_ints_ {0, 1, 2},
  algo_vals_(algo_strings_),
  debug_vals_(debug_ints_),
  hmmdir_arg_("m", "hmmdir", "directory in which to look for hmm model files", true, "", "string"),
  datadir_arg_("d", "datadir", "directory in which to look for non-sample-specific data (eg human germline seqs)", true, "", "string"),
  infile_arg_("i", "infile", "input (whitespace-separated) file", true, "", "string"),
  outfile_arg_("o", "outfile", "output csv file", true, "", "string"),
  algorithm_arg_("a", "algorithm", "algorithm to run", true, "", &algo_vals_),
  debug_arg_("g", "debug", "debug level", false, 0, &debug_vals_),
  n_best_events_arg_("n", "n_best_events", "number of candidate recombination events to write to file", true, -1, "int"),
  pair_arg_("p", "pair", "is this a pair hmm?", false),
  chunk_cache_arg_("c", "chunk-cache", "perform chunk caching?", false),
  str_headers_ {},
  int_headers_ {"k_v_min", "k_v_max", "k_d_min", "k_d_max"},
  str_list_headers_{"names", "seqs", "only_genes"}  // args that are passed as colon-separated lists
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
    cmd.add(chunk_cache_arg_);

    cmd.parse(argc, argv);

    // algorithm_ = algorithm_arg_.getValue();
    // debug_ = debug_arg_.getValue();
  } catch(ArgException &e) {
    cerr << "ERROR: " << e.error() << " for argument " << e.argId() << endl;
    throw;
  }

  for(auto & head : str_headers_)
    strings_[head] = vector<string>();
  for(auto & head : int_headers_)
    integers_[head] = vector<int>();

  ifstream ifs(infile());
  assert(ifs.is_open());
  string line;
  // get header line
  getline(ifs, line);
  stringstream ss(line);
  vector<string> headers;  // keep track of the file's column order
  while(!ss.eof()) {
    string head;
    ss >> head;
    if(head != "")   // make sure not to add a blank header at the end of the file
      headers.push_back(head);
  }
  while(getline(ifs, line)) {
    stringstream ss(line);
    string tmpstr;
    int tmpint;
    for(auto & head : headers) {
      if(str_headers_.find(head) != str_headers_.end()) {
        ss >> tmpstr;
        strings_[head].push_back(tmpstr);
      } else if(str_list_headers_.find(head) != str_list_headers_.end()) {
	ss >> tmpstr;
	str_lists_[head].push_back(Split(tmpstr));
      } else if(int_headers_.find(head) != int_headers_.end()) {
        ss >> tmpint;
        integers_[head].push_back(tmpint);
      } else {
        throw runtime_error("ERROR header " + head + "' not found");
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
vector<string> Args::Split(string argstr) {
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
// read input sequences from file and return as vector of sequences
vector<Sequences> GetSeqs(Args &args, Track *trk) {
  vector<Sequences> all_seqs;
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    Sequences seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq], trk);
      seqs.AddSeq(sq);
    }
    all_seqs.push_back(seqs);
  }
  assert(all_seqs.size() == args.str_lists_["names"].size());
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score);
void print_forward_scores(double ab_score, vector<double> single_scores);
// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  Args args(argc, argv);
  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  assert(ofs.is_open());
  if(args.algorithm() == "viterbi")
    ofs << "unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion,v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,score,seqs,errors" << endl;
  else
    ofs << "unique_ids,score,errors" << endl;

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track trk("NUKES", characters);
  vector<Sequences> seqs(GetSeqs(args, &trk));
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl);
  // clock_t cache_start(clock());
  // hmms.CacheAll();
  // cout << "t " << ((clock() - cache_start) / (double)CLOCKS_PER_SEC) << endl;
  // clock_t run_start(clock());

  for(size_t is = 0; is < seqs.size(); is++) {
    if(args.debug()) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][is], args.integers_["k_d_min"][is]);
    KSet kmax(args.integers_["k_v_max"][is], args.integers_["k_d_max"][is]);
    KBounds kbounds(kmin, kmax);

    JobHolder jh(gl, hmms, args.algorithm(), args.str_lists_["only_genes"][is]);
    jh.SetDebug(args.debug());
    jh.SetChunkCache(args.chunk_cache());
    jh.SetNBestEvents(args.n_best_events());

    Result result(kbounds);
    do {
      result = jh.Run(seqs[is], kbounds);
      if(result.could_not_expand()) {
	cout << "WARNING " << seqs[is].name_str() << " couldn't expand k bounds for " << kbounds.stringify() << endl;
      }
      kbounds = result.better_kbounds();
    } while(result.boundary_error() && !result.could_not_expand());

    double numerator(result.total_score()),total_score(result.total_score());;
    if(args.algorithm() == "forward" && seqs[is].n_seqs() > 1) {  // denominator in P(A,B) / (P(A) P(B)). NOTE now for C, D, E, F, ...
      vector<double> single_scores;  // NOTE log probs, not scores, but I haven't managed to finish switching over to the new terminology
      for (size_t iseq = 0; iseq < seqs[is].n_seqs(); ++iseq) {  // NOTE kind of confusing, but <is> is looping over queries, while <iseq> is looping over the sequences within the <is>th query.
	Result single_result = jh.Run(seqs[is][iseq], kbounds);  // result for single sequence
	total_score -= single_result.total_score();

	single_scores.push_back(single_result.total_score());
	if(single_result.boundary_error())
	  if (args.debug())  // boundary errors are bad for single sequences, or for pairs of sequences that are clonally related... but are totally expected for unrelated pairs
	    cout << "WARNING boundary errors for " << seqs[is][iseq].name() << " when together with " << seqs[is].name_str() << endl;
      }
      if(args.debug())
	print_forward_scores(numerator, single_scores);
    }

    if(args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {   // if we were asked for more events than we found
      if(result.events_.size() > 0)
        cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
        assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args, result.events_, seqs[is], total_score);
  }

  // cout << "t " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score) {
  if(args.algorithm() == "viterbi") {
    // assert(args.n_best_events() <= events.size());  // make sure we found at least as many valid events as were asked for with <--n_best_events>
    size_t n_max = min(size_t(args.n_best_events()), events.size());
    for(size_t ievt = 0; ievt < n_max; ++ievt) {
      RecoEvent *event = &events[ievt];
      string second_seq_name, second_seq;
      ofs  // be very, very careful to change this *and* the csv header above at the same time
	<< seqs.name_str(":")
	<< "," << event->genes_["v"]
	<< "," << event->genes_["d"]
	<< "," << event->genes_["j"]
	<< "," << event->insertions_["fv"]
	<< "," << event->insertions_["vd"]
	<< "," << event->insertions_["dj"]
	<< "," << event->insertions_["jf"]
	<< "," << event->deletions_["v_5p"]
	<< "," << event->deletions_["v_3p"]
	<< "," << event->deletions_["d_5p"]
	<< "," << event->deletions_["d_3p"]
	<< "," << event->deletions_["j_5p"]
	<< "," << event->deletions_["j_3p"]
	<< "," << event->score_
	<< "," << seqs.seq_str(":")
	<< "," << ""  //errors
	<< endl;
    }
  } else {
    ofs
        << seqs.name_str(":")
        << "," << total_score
        << "," << ""  //errors
        << endl;
  }
}
// ----------------------------------------------------------------------------------------
void print_forward_scores(double ab_score, vector<double> single_scores) {
  // NOTE need to update to work with k-hmms for k > 2
  printf("%70s %8.2f - %8.2f - %8.2f = %8.3f\n", "", ab_score, single_scores[0], single_scores[1], ab_score - single_scores[0] - single_scores[1]);
  feclearexcept(FE_UNDERFLOW | FE_OVERFLOW);
  double ab_prob(exp(ab_score));
  double a_prob(exp(single_scores[0]));
  double b_prob(exp(single_scores[1]));
  double bayes_factor(ab_prob / (a_prob*b_prob));
  if (fetestexcept(FE_UNDERFLOW | FE_OVERFLOW))
    printf("%70s under/overflow when leaving log space\n", "");
  else
    printf("%70s %8.1e / (%8.1e * %8.1e) = %8.1f\n", "", ab_prob, a_prob, b_prob, bayes_factor);
}

