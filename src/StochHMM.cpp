#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <map>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "StochHMMlib.h"
#include "job_holder.h"
#include "germlines.h"

#include "StochHMM_usage.h"
using namespace StochHMM;
using namespace std;
// NOTE use this to compile: export CPLUS_INCLUDE_PATH=/home/dralph/Dropbox/include/eigen-eigen-6b38706d90a9; make
#define STATE_MAX 1024
void run_casino();
double run_job(model *hmm, sequences *seqs);
void print_output(multiTraceback*, string&);
void print_output(traceback_path*, string&);
void print_posterior(trellis&);

// init command-line options
opt_parameters commandline[] = {
  // {"-model:-m"     ,OPT_STRING     ,false  ,"",    {}},
  // {"-seq"          ,OPT_STRING     ,false   ,"",    {}},
  {"-hmmtype"      ,OPT_STRING     ,false  ,"single", {"single", "pair"}},
  // {"-only_genes"   ,OPT_STRING     ,false  ,"",   {}},
  {"-hmmdir"       ,OPT_STRING     ,false  ,"",   {}},
  {"-infile"       ,OPT_STRING     ,true  ,"",   {}},
  {"-debug"        ,OPT_INT        ,false  ,"",    {}},
  // {"-k_v_guess"    ,OPT_INT        ,false   ,"",    {}},
  // {"-k_d_guess"    ,OPT_INT        ,false   ,"",    {}},
  // {"-v_fuzz"       ,OPT_INT        ,false  ,"3",    {}},
  // {"-d_fuzz"       ,OPT_INT        ,false  ,"3",    {}},
  //Non-Stochastic Decoding
  {"-viterbi"      ,OPT_NONE       ,false  ,"",    {}},
  {"-forward"      ,OPT_NONE       ,false  ,"",    {}},
  {"-posterior"    ,OPT_NONE       ,false  ,"",    {}},
  //Stochastic Decoding
  {"-stochastic"   ,OPT_FLAG       ,false  ,"",    {"viterbi"}},
  {"-repetitions:-rep",OPT_INT     ,false  ,"1000",{}},
  //Output Files and Formats
  {"-path:-p"      ,OPT_STRING     ,false  ,"",    {}},
  {"-label:-l"     ,OPT_STRING     ,false  ,"",    {}},
  {"-hits"         ,OPT_STRING     ,false  ,"",    {}},
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);  //Stores the number of options in opt
options opt;  //Global options for parsed command-line options

// ----------------------------------------------------------------------------------------
vector<sequences*> GetSeqs(string seqfname, track *trk) {
  vector<sequences*> all_seqs;
  ifstream ifs(seqfname);
  assert(ifs.is_open());
  while(!ifs.eof()) {
    sequences *seqs = new sequences;
    for (size_t iseq=0; iseq<trk->n_seqs; iseq++) {  // for pair hmm, we push back two seqs for each track
      sequence *sq = new(nothrow) sequence(false);
      assert(sq);
      sq->getFasta(ifs, trk);
      seqs->addSeq(sq);
    }
    all_seqs.push_back(seqs);
  }
  ifs.close();
  return all_seqs;
}
// ----------------------------------------------------------------------------------------
class Args {
public:
  Args(string fname);
  map<string, vector<string> > strings_;
  map<string, vector<int> > integers_;
  set<string> str_headers_, int_headers_;
  // vector<string> names, second_names, seqs, second_seqs, only_geneses;
  // vector<int> k_v_guesses, k_d_guesses, v_fuzzes, d_fuzzes;
};
Args::Args(string fname):
  str_headers_{"only_genes", "name", "seq", "second_name", "second_seq"},
  int_headers_{"k_v_guess", "k_d_guess", "v_fuzz", "d_fuzz"}
{
  for (auto &head: str_headers_)
    strings_[head] = vector<string>();
  for (auto &head: int_headers_)
    integers_[head] = vector<int>();
  
  ifstream ifs(fname);
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
      } else if (int_headers_.find(head) != int_headers_.end()) {
	ss >> tmpint;
	integers_[head].push_back(tmpint);
      } else {
	assert(0);
      }
    }
  }
  // for (size_t i=0; i<strings_["seq"].size(); ++i) {
  //   for (auto &head: str_headers_) {
  //     cout << head << "  " << strings_[head][i] << endl;
  //   }
  //   for (auto &head: int_headers_)
  //     cout << head << "  " << integers_[head][i] << endl;
  //   cout << endl;
  // }
}
// ----------------------------------------------------------------------------------------
vector<sequences*> GetSeqs(Args &args, track *trk) {
  vector<sequences*> all_seqs;
  for (size_t iseq=0; iseq<args.strings_["seq"].size(); ++iseq) {
    sequences *seqs = new sequences;
    sequence *sq = new(nothrow) sequence(args.strings_["seq"][iseq], trk, args.strings_["name"][iseq]);
    assert(sq);
    seqs->addSeq(sq);
    if (trk->n_seqs == 2) {
      assert(args.strings_["second_seq"][iseq].size() > 0);
      assert(args.strings_["second_seq"][iseq] != "x");
      sequence *second_sq = new(nothrow) sequence(args.strings_["second_seq"][iseq], trk, args.strings_["second_name"][iseq]);
      assert(second_sq);
      seqs->addSeq(second_sq);
    } else {
      assert(args.strings_["second_seq"][iseq] == "x");  // er, not really necessary, I suppose...
    }

    all_seqs.push_back(seqs);
  }
  return all_seqs;
}
// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  GermLines gl;
  srand(time(NULL));
  opt.set_parameters(commandline, opt_size, usage);
  opt.parse_commandline(argc,argv);
  string algorithm;
  if (opt.isSet("-viterbi"))
    algorithm = "viterbi";
  else if (opt.isSet("-forward"))
    algorithm = "forward";

  Args args(opt.sopt("-infile"));

  size_t n_seqs_per_track(opt.sopt("-hmmtype")=="pair" ? 2 : 1);
  vector<string> characters{"A","C","G","T"};
  track trk("NUKES", n_seqs_per_track, characters);
  vector<sequences*> seqs(GetSeqs(args, &trk));
  HMMHolder hmms(opt.isSet("-hmmdir") ? opt.sopt("-hmmdir") : "./bcell", n_seqs_per_track);

  assert(seqs.size() == 1);  // for the moment this makes the most sense.
  assert(seqs.size() == args.strings_["name"].size());
  for (size_t is=0; is<seqs.size(); is++) {
    JobHolder jh(n_seqs_per_track, algorithm, seqs[is], &hmms, args.strings_["only_genes"][is]);
    if (opt.isSet("-debug"))
      jh.SetDebug(opt.iopt("-debug"));
    int k_v_guess = args.integers_["k_v_guess"][is];
    int k_d_guess = args.integers_["k_d_guess"][is];
    int v_fuzz = args.integers_["v_fuzz"][is];
    int d_fuzz = args.integers_["d_fuzz"][is];
    jh.Run(max(k_v_guess - v_fuzz, 1), 2*v_fuzz, max(k_d_guess - d_fuzz, 1), 2*d_fuzz);
  }
  return 0;
}

// ----------------------------------------------------------------------------------------
void run_casino() {
  model hmm;
  hmm.import("examples/Dice.hmm");
  ifstream ifs("examples/Dice_short.fa");
  assert(ifs.is_open());
  while (!ifs.eof()) {
    sequences seqs;
    sequence *seq = new sequence;
    seq->getFasta(ifs, hmm.getTrack((size_t)0));
    seqs.addSeq(seq);
    run_job(&hmm, &seqs);
  }
  ifs.close();
}

// ----------------------------------------------------------------------------------------
double run_job(model *hmm, sequences *seqs) {
  double score(-INFINITY);
  trellis trell(hmm, seqs);
  if (opt.isSet("-viterbi")) {
    trell.viterbi();
    traceback_path path(hmm);
    trell.traceback(path);
    print_output(&path, seqs->getHeader());
    score = path.getScore();
  } else if (opt.isSet("-forward")) {
    trell.forward();
    score = trell.getForwardProbability();
  } else if (opt.isSet("-posterior")) {
    trell.posterior();
    if (opt.isSet("-path") || opt.isSet("-label")) {  //If we need a posterior traceback b/c path or label is defined
      traceback_path path(hmm);
      trell.traceback_posterior(path);  // NOTE posterior score is not set with this option
      print_output(&path, seqs->getHeader());
    } else {
      score = trell.getForwardProbability();
      // std::cout << seqs->stringifyWOHeader();
      // cout << "  fwd score: " << trell.getBackwardProbability() << endl;
      // cout << "  bwd score: " << trell.getForwardProbability() << endl;
      // print_posterior(trell);
    }
  } else if (opt.isSet("-stochastic")) {
    assert(opt.isFlagSet("-stochastic", "viterbi"));  // removed the other ones for the time being
    trell.stochastic_viterbi();
    multiTraceback paths;
    trell.stochastic_traceback(paths, opt.iopt("-rep"));
    print_output(&paths, seqs->getHeader());
  } else {
    cout << "ERROR no option" << endl;
    assert(0);
  }
  return score;
}

// ----------------------------------------------------------------------------------------
void print_output(multiTraceback* tb, string& header) {
  tb->finalize();
  bool previous(true);
  if (opt.isSet("-hits")){
    tb->print_hits();
    previous=false;
  }
  if (opt.isSet("-label")){
    tb->print_label();
    previous=false;
  }
  if (opt.isSet("-path") || previous) {  //Print path by default if nothing else is set
    tb->print_path();
  }
}

// ----------------------------------------------------------------------------------------
void print_output(traceback_path* tb, string& header){
  bool previous(true);
  if (opt.isSet("-label")) {
    cout << ">" << header ;
    cout << "\tScore: " << tb->getScore() << endl;
    tb->print_label();
    previous=false;
  }
  if (opt.isSet("-path") || previous){
    cout << ">" << header ;
    cout << "\tScore: " << tb->getScore() << endl;
    tb->print_path();
  }
  return;
}

// ----------------------------------------------------------------------------------------
//Print the posterior probabilities for each state at each position
//Each state is in separate column
//Each row is on different row
void print_posterior(trellis& trell) {
  model* hmm = trell.getModel();
  double_2D* table = trell.getPosteriorTable();
  size_t state_size = hmm->state_size();
  char cstr[200];
  string output;
  output += "Posterior Probabilities Table\n";
  output += "Model:\t" + hmm->getName() + "\n";
  output += "Sequence:\t" + trell.getSeq()->getHeader() + "\n";
  sprintf(cstr, "Probability of Sequence from Forward: Natural Log'd\t%f\n", trell.getForwardProbability());
  output += cstr;
  sprintf(cstr, "Probability of Sequence from Backward:Natural Log'd\t%f\n", trell.getBackwardProbability());
  output += cstr;
  output+= "Position";
  for(size_t i=0; i<state_size; ++i) { // print each state name
    output += "\t" + hmm->getStateName(i);
  }
  output += "\n";
  cout <<  output;

  for(size_t position=0; position<table->size(); ++position) {
    sprintf(cstr, "%ld", position+1);
    output = cstr;
    for (size_t st=0; st<state_size; st++) {
      float val = exp((*table)[position][st]);
      if (val <= 0.001){
	output += "\t0";
      } else if (val == 1.0) {
	output += "\t1";
      } else {
	sprintf(cstr,"\t%.3f", exp((*table)[position][st]));
	output+= cstr;
      }
    }
    output += "\n";
    cout << output;
  }
  cout << endl;
  return;
}
