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
#define STATE_MAX 1024  // TODO reduce the state max to something reasonable?
void StreamOutput(ofstream &ofs, options &opt, vector<RecoEvent> &events, sequences &seqs, double total_score);

// ----------------------------------------------------------------------------------------
// init command-line options
opt_parameters commandline[] = {
  {"--hmmdir",            OPT_STRING, false, "./bcell", {}},
  {"--infile",            OPT_STRING, true,  "",        {}},
  {"--outfile",           OPT_STRING, true,  "",        {}},
  {"--algorithm",         OPT_STRING, true,  "viterbi", {}},
  {"--debug",             OPT_INT,    false, "0",       {}},
  {"--n_best_events",     OPT_INT,    true,  "",        {}},
  {"--pair",              OPT_INT,    false, "0",       {}},  // holy crap why is flag not a flag?
  {"--single_gene_probs", OPT_INT,    false, "0",       {}},  // get the probability that each single gene version occurs in each query sequence and write to file (for preclustering)
};

int opt_size=sizeof(commandline)/sizeof(commandline[0]);  //Stores the number of options in opt
options opt;  //Global options for parsed command-line options

// ----------------------------------------------------------------------------------------
// class for reading csv input info
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
  srand(time(NULL));
  opt.set_parameters(commandline, opt_size, "");
  opt.parse_commandline(argc,argv);
  assert(opt.sopt("--algorithm") == "viterbi" || opt.sopt("--algorithm") == "forward");
  assert(opt.iopt("--pair") == 0 || opt.iopt("--pair") == 1);
  assert(opt.iopt("--debug") >=0 && opt.iopt("--debug") < 3);
  if (opt.iopt("--single_gene_probs"))
    assert(!opt.iopt("--pair"));
  Args args(opt.sopt("--infile"));

  // write csv output headers
  ofstream ofs;
  ofs.open(opt.sopt("--outfile"));
  assert(ofs.is_open());
  if (opt.iopt("--single_gene_probs")) {
    ofs << "unique_id,scores" << endl;
  } else {
    if (opt.sopt("--algorithm") == "viterbi")
      ofs << "unique_id,second_unique_id,v_gene,d_gene,j_gene,vd_insertion,dj_insertion,v_3p_del,d_5p_del,d_3p_del,j_5p_del,score,seq,second_seq" << endl;
    else
      ofs << "unique_id,second_unique_id,score" << endl;
  }

  // init some stochhmm infrastructure
  size_t n_seqs_per_track(opt.iopt("--pair") ? 2 : 1);
  vector<string> characters{"A","C","G","T"};
  track trk("NUKES", n_seqs_per_track, characters);
  vector<sequences*> seqs(GetSeqs(args, &trk));
  GermLines gl;
  HMMHolder hmms(opt.sopt("--hmmdir"), n_seqs_per_track, gl);

  assert(seqs.size() == args.strings_["name"].size());
  for (size_t is=0; is<seqs.size(); is++) {
    int k_start = max(1, args.integers_["k_v_guess"][is] - args.integers_["v_fuzz"][is]);
    int d_start = max(1, args.integers_["k_d_guess"][is] - args.integers_["d_fuzz"][is]);
    int n_k_v = 2 * args.integers_["v_fuzz"][is];  // NOTE this is the *maximum* fuzz allowed -- job_holder::Run is allowed to skip ksets that don't make sense
    int n_k_d = 2 * args.integers_["d_fuzz"][is];  // TODO oh wait shouldn't k_d be allowed to be zero?

    JobHolder jh(gl, hmms, opt.sopt("--algorithm"), args.strings_["only_genes"][is]);
    jh.SetDebug(opt.iopt("--debug"));
    jh.SetNBestEvents(opt.iopt("--n_best_events"));

    Result result = jh.Run(*seqs[is], k_start, n_k_v, d_start, n_k_d);
    double score(result.total_score_);
    if (opt.sopt("--algorithm") == "forward" && opt.iopt("--pair") && !opt.iopt("--single_gene_probs")) {
      assert(seqs[is]->size() == 2);
      Result result_a = jh.Run((*seqs[is])[0], k_start, n_k_v, d_start, n_k_d);
      Result result_b = jh.Run((*seqs[is])[1], k_start, n_k_v, d_start, n_k_d);
      score = score - result_a.total_score_ - result_b.total_score_;
    }

    if (opt.iopt("--single_gene_probs")) {
      assert(seqs[is]->size() == 1);
      jh.WriteBestGeneProbs(ofs, (*seqs[is])[0].name_);
    } else {
      StreamOutput(ofs, opt, result.events_, *seqs[is], score);
    }
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, options &opt, vector<RecoEvent> &events, sequences &seqs, double total_score) {
  if (opt.sopt("--algorithm") == "viterbi") {
    assert(opt.iopt("--n_best_events") <= events.size());
    for (size_t ievt=0; ievt<opt.iopt("--n_best_events"); ++ievt) {
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
	<< endl;
    }
  } else {
    assert(seqs.size() == 2);  // er, at least for the moment
    ofs
      << seqs[0].name_
      << "," << seqs[1].name_
      << "," << total_score << endl;
  }
}
