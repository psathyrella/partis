//Copyright (c) 2007-2012 Paul C Lott
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iomanip>
#include <time.h>
#include <fstream>
#include "StochHMMlib.h"

#include "StochHMM_usage.h"
using namespace StochHMM;

#define STATE_MAX 1024


void import_model(model&);
void import_sequence(model&);

void print_output(multiTraceback*, std::string&);
void print_output(std::vector<traceback_path>&, std::string&);
void print_output(traceback_path*, std::string&);
void print_posterior(trellis&);

// set command-line options
opt_parameters commandline[]={
  {"-help:-h"     ,OPT_NONE       ,false  ,"",    {}},
  //Required
  {"-model:-m"    ,OPT_STRING     ,true   ,"",    {}},
  {"-seq:-s:-track",OPT_STRING    ,true   ,"",    {}},
  //Non-Stochastic Decoding
  {"-viterbi"     ,OPT_NONE       ,false  ,"",    {}},
  {"-posterior"   ,OPT_STRING         ,false  ,"",    {}},
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
seqTracks jobs;  //seqTracks stores multiple jobs(model and multiple sequences)
StateFuncs default_functions;  // this will automatically initialize all the Univariate and Multivariate PDFs

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  opt.set_parameters(commandline, opt_size, usage);
  opt.parse_commandline(argc,argv);
  assert(opt.isSet("-model"));
  assert(opt.isSet("-seq"));

  model hmm;
  hmm.import(opt.sopt("-model"), &default_functions);
  jobs.loadSeqs(hmm, opt.sopt("-seq"), FASTA);  // calls importJobs, which calls getNext() *once* but no more.

  seqJob *job = jobs.getJob();  // job consists of a model and associated sequences
  while (job) {  // For each job(sequence) perform the analysis
    sequences *seqs(job->getSeqs());
    trellis trell(&hmm, seqs);
    if (opt.isSet("-posterior")) {
      trell.posterior();
      if (opt.isSet("-path") || opt.isSet("-label")) {  //If we need a posterior traceback b/c path,label,or GFF is defined
	traceback_path path(&hmm);
	trell.traceback_posterior(path);
	print_output(&path, seqs->getHeader());
      } else {
	print_posterior(trell);
      }
    } else if (opt.isSet("-viterbi")) {
      trell.viterbi();
      traceback_path path(&hmm);
      trell.traceback(path);
      print_output(&path, seqs->getHeader());
    } else if (opt.isSet("-stochastic")) {
      assert(opt.isFlagSet("-stochastic", "viterbi"));  // removed the other ones for the time being
      trell.stochastic_viterbi();
      multiTraceback paths;
      trell.stochastic_traceback(paths, opt.iopt("-rep"));
      print_output(&paths, seqs->getHeader());
    }
    job = jobs.getJob();  // get next job
  }
}

// ----------------------------------------------------------------------------------------
void print_output(multiTraceback* tb, std::string& header){
    
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
    
  //Print path by default if nothing else is set
  if (opt.isSet("-path") || previous){
    tb->print_path();
  }
}

// ----------------------------------------------------------------------------------------
void print_output(traceback_path* tb, std::string& header){
    
  bool previous(true);
    
  if (opt.isSet("-label")){
    std::cout << ">" << header ;
    std::cout << "\tScore: " << tb->getScore() << std::endl;
    tb->print_label();
    previous=false;
  }
    
  if (opt.isSet("-path") || previous){
    std::cout << ">" << header ;
    std::cout << "\tScore: " << tb->getScore() << std::endl;
    tb->print_path();
  }
    
  return;
}


// ----------------------------------------------------------------------------------------
//Print the posterior probabilities for each state at each position
//Each state is in separate column
//Each row is on different row
void print_posterior(trellis& trell){
  model* hmm = trell.getModel();
  double_2D* table = trell.getPosteriorTable();
  size_t state_size = hmm->state_size();
  char cstr[200];
        
  std::string output;
  output+="Posterior Probabilities Table\n";
  output+="Model:\t" + hmm->getName() + "\n";
  output+="Sequence:\t" + trell.getSeq()->getHeader() + "\n";
  sprintf(cstr, "Probability of Sequence from Forward: Natural Log'd\t%f\n",trell.getForwardProbability());
  output+= cstr;
  sprintf(cstr, "Probability of Sequence from Backward:Natural Log'd\t%f\n",trell.getBackwardProbability());
  output+= cstr;
  output+= "Position";
  for(size_t i=0;i< state_size; ++i){
    output+= "\t" + hmm->getStateName(i);
  }
  output+="\n";
        
  std::cout <<  output;
        

  for(size_t position = 0; position < table->size(); ++position){
    sprintf(cstr, "%ld", position+1);
    output= cstr;
    for (size_t st = 0 ; st < state_size ; st++){
      float val  = exp((*table)[position][st]);
      if (val<= 0.001){
	output+="\t0";
      }
      else if (val == 1.0){
	output+="\t1";
      }
      else{
	sprintf(cstr,"\t%.3f", exp((*table)[position][st]));
	output+= cstr;
      }

    }
    output+="\n";
    std::cout << output;
  }

  std::cout << std::endl;
        
  return;
        
}

