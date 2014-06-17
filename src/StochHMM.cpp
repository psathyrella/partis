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

#include "StochHMM_usage.h"
using namespace StochHMM;
using namespace std;

#define STATE_MAX 1024

void run_casino();
double run(JobHolder &jh, size_t k_v, size_t k_d);
double run_region(model *hmm, sequences *seqs);

// ----------------------------------------------------------------------------------------
void print_output(multiTraceback*, string&);
// void print_output(vector<traceback_path>&, string&);
void print_output(traceback_path*, string&);
void print_posterior(trellis&);

// init command-line options
opt_parameters commandline[]={
  {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
  //Required
  {"-model:-m"    ,OPT_STRING     ,false   ,"",    {}},
  {"-seq:-s:-track",OPT_STRING    ,true   ,"",    {}},
  //Non-Stochastic Decoding
  {"-viterbi"     ,OPT_NONE       ,false  ,"",    {}},
  {"-forward"     ,OPT_NONE       ,false  ,"",    {}},
  {"-posterior"   ,OPT_NONE       ,false  ,"",    {}},
  //Stochastic Decoding
  {"-stochastic"  ,OPT_FLAG       ,false  ,"",    {"viterbi"}},
  {"-repetitions:-rep",OPT_INT    ,false  ,"1000",{}},
  //Output Files and Formats
  {"-path:-p"     ,OPT_STRING     ,false  ,"",    {}},
  {"-label:-l"    ,OPT_STRING     ,false  ,"",    {}},
  {"-hits"        ,OPT_STRING     ,false  ,"",    {}},
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);  //Stores the number of options in opt
options opt;  //Global options for parsed command-line options
StateFuncs default_functions;  // this will automatically initialize all the Univariate and Multivariate PDFs
vector<string> regions{"v","d","j"};

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  opt.set_parameters(commandline, opt_size, usage);
  opt.parse_commandline(argc,argv);
  assert(opt.isSet("-seq"));

  if (opt.isSet("-model")) {
    assert(opt.sopt("-model") == "dice");
    run_casino();
    return 0;
  }
  JobHolder jh("bcell", regions, opt.sopt("-seq"));
  double best_score(-INFINITY);
  size_t best_k_v,best_k_d;
  for (size_t k_v=295; k_v<297; ++k_v) {
    for (size_t k_d=16; k_d<18; ++k_d) {
      // k_v = 296; k_d = 17;
      double score = run(jh, k_v, k_d);
      cout
	<< setw(12) << k_v
	<< setw(12) << k_d
	<< setw(12) << score
	<< endl;
      if (score > best_score) {
	best_score = score;
	best_k_v = k_v;
	best_k_d = k_d;
      }
    }
  }
  cout
    << "best: "
    << setw(12) << best_k_v
    << setw(12) << best_k_d
    << setw(12) << best_score
    << endl;
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
    run_region(&hmm, &seqs);
  }
  ifs.close();
}

// ----------------------------------------------------------------------------------------
double run(JobHolder &jh, size_t k_v, size_t k_d) {
  double score(-INFINITY);
  map<string,sequences> subseqs = jh.GetSubSeqs(k_v, k_d);
  for (auto &region : regions) {
    double tmp_score = run_region(jh.hmm(region), &subseqs[region]);
    if (score==-INFINITY) {
      score = tmp_score;
    } else {
      score += tmp_score;
    }
  }
  return score;
}

// ----------------------------------------------------------------------------------------
double run_region(model *hmm, sequences *seqs) {
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
