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

// #include "tclap/CmdLine.h"
// using namespace TCLAP;
// 
#include "StochHMM_usage.h"
using namespace StochHMM;
using namespace std;
#define STATE_MAX 1024  // TODO reduce the state max to something reasonable?
void StreamOutput(ofstream &ofs, options &opt, vector<RecoEvent> &events, sequences &seqs, double total_score);

// ----------------------------------------------------------------------------------------
// init command-line options
opt_parameters commandline[] = {
  {"--hmmdir",            OPT_STRING, false, "",   {}},
  {"--datadir",           OPT_STRING, false, "",   {}},
  {"--infile",            OPT_STRING, true,  "",   {}},
  {"--outfile",           OPT_STRING, true,  "",   {}},
  {"--algorithm",         OPT_STRING, true,  "",   {}},
  {"--debug",             OPT_INT,    false, "0",  {}},
  {"--n_best_events",     OPT_INT,    true,  "",   {}},
  {"--pair",              OPT_INT,    false, "0",  {}},  // holy crap why is flag not a flag?
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);  //Stores the number of options in opt
options opt;  //Global options for parsed command-line options

// ----------------------------------------------------------------------------------------
// class for reading csv input file
class Args {
public:
  Args(string fname);

  string all_only_genes_;
  map<string, vector<string> > strings_;
  map<string, vector<int> > integers_;
  set<string> str_headers_, int_headers_;
};

// ----------------------------------------------------------------------------------------
Args::Args(string fname):
  str_headers_{"only_genes", "name", "seq", "second_name", "second_seq"},
  int_headers_{"k_v_min", "k_v_max", "k_d_min", "k_d_max"}
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
// // ----------------------------------------------------------------------------------------
// 	// Wrap everything in a try block.  Do this every time, 
// 	// because exceptions will be thrown for problems. 
// 	try {  

// 	// Define the command line object.
// 	CmdLine cmd("Command description message", ' ', "0.9");

// 	// Define a value argument and add it to the command line.
// 	ValueArg<string> nameArg("n","name","Name to print",true,"homer","string");
// 	cmd.add( nameArg );

// 	// Define a switch and add it to the command line.
// 	SwitchArg reverseSwitch("r","reverse","Print name backwards", false);
// 	cmd.add( reverseSwitch );

// 	// Parse the args.
// 	cmd.parse( argc, argv );

// 	// Get the value parsed by each arg. 
// 	string name = nameArg.getValue();
// 	bool reverseName = reverseSwitch.getValue();

// 	// Do what you intend too...
// 	if ( reverseName )
// 	{
// 		reverse(name.begin(),name.end());
// 		cout << "My name (spelled backwards) is: " << name << endl;
// 	}
// 	else
// 		cout << "My name is: " << name << endl;


// 	} catch (ArgException &e)  // catch any exceptions
// 	{ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

// // ----------------------------------------------------------------------------------------
  srand(time(NULL));
  opt.set_parameters(commandline, opt_size, "");
  opt.parse_commandline(argc,argv);
  assert(opt.sopt("--algorithm") == "viterbi" || opt.sopt("--algorithm") == "forward");
  assert(opt.iopt("--pair") == 0 || opt.iopt("--pair") == 1);
  assert(opt.iopt("--debug") >=0 && opt.iopt("--debug") < 3);
  Args args(opt.sopt("--infile"));

  // write csv output headers
  ofstream ofs;
  ofs.open(opt.sopt("--outfile"));
  assert(ofs.is_open());
  if (opt.sopt("--algorithm") == "viterbi")
    ofs << "unique_id,second_unique_id,v_gene,d_gene,j_gene,vd_insertion,dj_insertion,v_3p_del,d_5p_del,d_3p_del,j_5p_del,score,seq,second_seq,errors" << endl;
  else
    ofs << "unique_id,second_unique_id,score,errors" << endl;

  // init some stochhmm infrastructure
  size_t n_seqs_per_track(opt.iopt("--pair") ? 2 : 1);
  vector<string> characters{"A","C","G","T"};
  track trk("NUKES", n_seqs_per_track, characters);
  vector<sequences*> seqs(GetSeqs(args, &trk));
  GermLines gl(opt.sopt("--datadir"));
  HMMHolder hmms(opt.sopt("--hmmdir"), n_seqs_per_track, gl);

  assert(seqs.size() == args.strings_["name"].size());
  for (size_t is=0; is<seqs.size(); is++) {
    if (opt.iopt("--debug")) cout << "  ---------" << endl;
    if (opt.iopt("--pair")) assert(seqs[is]->n_seqs() == 2);
    KSet kmin(args.integers_["k_v_min"][is], args.integers_["k_d_min"][is]);
    KSet kmax(args.integers_["k_v_max"][is], args.integers_["k_d_max"][is]);
    KBounds kbounds(kmin, kmax);

    JobHolder jh(gl, hmms, opt.sopt("--algorithm"), args.strings_["only_genes"][is]);
    jh.SetDebug(opt.iopt("--debug"));
    jh.SetNBestEvents(opt.iopt("--n_best_events"));

    Result result(kbounds);
    do {
      result = jh.Run(*seqs[is], kbounds);
      if (result.could_not_expand()) cout << "WARNING " << (*seqs[is])[0].name() << ((*seqs[is]).n_seqs()==2 ? (*seqs[is])[1].name() : "") << " couldn't expand k bounds for " << kbounds.stringify() << endl;
      kbounds = result.better_kbounds();
    } while (result.boundary_error() && !result.could_not_expand());

    double score(result.total_score());
    if (opt.sopt("--algorithm") == "forward" && opt.iopt("--pair")) {
      Result result_a = jh.Run((*seqs[is])[0], kbounds);  // denominator in P(A,B) / (P(A) P(B))
      Result result_b = jh.Run((*seqs[is])[1], kbounds);
      
      if (result_a.boundary_error() || result_b.boundary_error()) cout << "WARNING boundary errors for " << (*seqs[is])[0].name() << " " << (*seqs[is])[1].name() << endl;
      if (opt.iopt("--debug")) printf("%70s %8.2f - %8.2f - %8.2f = %8.3f\n", "", score, result_a.total_score(), result_b.total_score(), score - result_a.total_score() - result_b.total_score());
      if (opt.iopt("--debug")) printf("%70s %8.1e / %8.1e / %8.1e = %8.1f\n", "", exp(score), exp(result_a.total_score()), exp(result_b.total_score()), exp(score) / (exp(result_a.total_score())*exp(result_b.total_score())));
      score = score - result_a.total_score() - result_b.total_score();
    }

    if (opt.sopt("--algorithm") == "viterbi" && size_t(opt.iopt("--n_best_events")) > result.events_.size()) {  // if we were asked for more events than we found
      if (result.events_.size() > 0)
	cout << "WARNING asked for " << opt.iopt("--n_best_events") << " events but only found " << result.events_.size() << endl;
      else
	assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, opt, result.events_, *seqs[is], score);
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, options &opt, vector<RecoEvent> &events, sequences &seqs, double total_score) {
  if (opt.sopt("--algorithm") == "viterbi") {
    // assert(opt.iopt("--n_best_events") <= events.size());  // make sure we found at least as many valid events as were asked for with <--n_best_events>
    size_t n_max = min(size_t(opt.iopt("--n_best_events")), events.size());
    for (size_t ievt=0; ievt<n_max; ++ievt) {
      RecoEvent *event = &events[ievt];
      string second_seq_name,second_seq;
      if (opt.iopt("--pair")) {
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
