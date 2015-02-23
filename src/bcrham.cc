#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <map>
#include <iomanip>
// #include <ctime>
#include <fstream>
#include <cfenv>
#include "jobholder.h"
#include "germlines.h"
#include "text.h"
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
  string hmmdir() { return hmmdir_arg_.getValue(); }
  string datadir() { return datadir_arg_.getValue(); }
  string infile() { return infile_arg_.getValue(); }
  string outfile() { return outfile_arg_.getValue(); }
  float hamming_fraction_cutoff() { return hamming_fraction_cutoff_arg_.getValue(); }
  string algorithm() { return algorithm_arg_.getValue(); }
  // string algorithm() { return algorithm_; }
  int debug() { return debug_arg_.getValue(); }
  // int debug() { return debug_; }
  int n_best_events() { return n_best_events_arg_.getValue(); }
  bool chunk_cache() { return chunk_cache_arg_.getValue(); }
  bool hagglom() { return hagglom_arg_.getValue(); }

  // command line arguments
  vector<string> algo_strings_;
  vector<int> debug_ints_;
  ValuesConstraint<string> algo_vals_;
  ValuesConstraint<int> debug_vals_;
  ValueArg<string> hmmdir_arg_, datadir_arg_, infile_arg_, outfile_arg_, algorithm_arg_;
  ValueArg<float> hamming_fraction_cutoff_arg_;
  ValueArg<int> debug_arg_, n_best_events_arg_;
  SwitchArg chunk_cache_arg_, hagglom_arg_;

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
              hamming_fraction_cutoff_arg_("f", "hamming-fraction-cutoff", "hamming fraction cutoff for clustering", true, -1.0, "float"),
              debug_arg_("g", "debug", "debug level", false, 0, &debug_vals_),
              n_best_events_arg_("n", "n_best_events", "number of candidate recombination events to write to file", true, -1, "int"),
              chunk_cache_arg_("c", "chunk-cache", "perform chunk caching?", false),
              hagglom_arg_("z", "hierarch-agglom", "", false),
              str_headers_ {},
              int_headers_ {"k_v_min", "k_v_max", "k_d_min", "k_d_max"},
str_list_headers_ {"names", "seqs", "only_genes"} { // args that are passed as colon-separated lists
  try {
    CmdLine cmd("ham -- the fantastic HMM compiler", ' ', "");
    cmd.add(hmmdir_arg_);
    cmd.add(datadir_arg_);
    cmd.add(infile_arg_);
    cmd.add(outfile_arg_);
    cmd.add(hamming_fraction_cutoff_arg_);
    cmd.add(algorithm_arg_);
    cmd.add(debug_arg_);
    cmd.add(n_best_events_arg_);
    cmd.add(chunk_cache_arg_);
    cmd.add(hagglom_arg_);

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
        str_lists_[head].push_back(SplitString(tmpstr, ":"));
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
// read input sequences from file and return as vector of sequences
vector<Sequences> GetSeqs(Args &args, Track *trk) {
  vector<Sequences> all_seqs;
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    Sequences seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(trk, args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq]);
      seqs.AddSeq(sq);
    }
    all_seqs.push_back(seqs);
  }
  assert(all_seqs.size() == args.str_lists_["names"].size());
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score, string errors);
void print_forward_scores(double numerator, vector<double> single_scores, double bayes_factor);
void hierarch_agglom(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs);
void glomerate(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs, map<string, Sequences> &info, map<string, KBounds> &kbinfo, map<string, vector<string> > &only_genes);
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
  vector<Sequences> qry_seq_list(GetSeqs(args, &trk));
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl);
  // hmms.CacheAll();

  if(args.hagglom()) {
    hierarch_agglom(hmms, gl, qry_seq_list, args, ofs);
    ofs.close();
    return 0;
  }
  
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    if(args.debug()) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kbounds(kmin, kmax);
    Sequences qry_seqs(qry_seq_list[iqry]);

    JobHolder jh(gl, hmms, args.algorithm(), args.str_lists_["only_genes"][iqry]);
    jh.SetDebug(args.debug());
    jh.SetChunkCache(args.chunk_cache());
    jh.SetNBestEvents(args.n_best_events());

    Result result(kbounds);
    vector<Result> denom_results(qry_seqs.n_seqs(), result);  // only used for forward if n_seqs > 1
    double numerator(-INFINITY);  // numerator in P(A,B,C,...) / (P(A)P(B)P(C)...)
    double bayes_factor(-INFINITY); // final result
    vector<double> single_scores(qry_seqs.n_seqs(), -INFINITY);  // NOTE log probs, not scores, but I haven't managed to finish switching over to the new terminology
    bool stop(false);
    string errors;
    do {
      errors = "";
      // clock_t run_start(clock());
      if(args.debug()) cout << "       ----" << endl;
      result = jh.Run(qry_seqs, kbounds);
      numerator = result.total_score();
      bayes_factor = numerator;
      if(args.algorithm() == "forward" && qry_seqs.n_seqs() > 1) {  // calculate factors for denominator
        for(size_t iseq = 0; iseq < qry_seqs.n_seqs(); ++iseq) {
          denom_results[iseq] = jh.Run(qry_seqs[iseq], kbounds);  // result for a single sequence
          single_scores[iseq] = denom_results[iseq].total_score();
          bayes_factor -= single_scores[iseq];
        }
      }

      kbounds = result.better_kbounds();
      for(auto & res : denom_results)
        kbounds = kbounds.LogicalOr(res.better_kbounds());

      stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
      for(auto & res : denom_results)
        stop &= !res.boundary_error() || res.could_not_expand();
      if(args.debug() && !stop)
        cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
      // cout << "      time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;
      if(result.boundary_error())
        errors = "boundary";
      for(auto & res : denom_results) // NOTE <errors> will still just have one "boundary" in it even if multiple results had boundary errors
        if(res.boundary_error())
          errors = "boundary";
    } while(!stop);

    // if(result.could_not_expand())
    //   cout << "WARNING " << qry_seqs.name_str() << " couldn't expand k bounds for " << kbounds.stringify() << endl;
    // if(single_result.boundary_error())
    //   cout << "WARNING boundary errors for " << qry_seqs[iseq].name() << " when together with " << qry_seqs.name_str() << endl;

    if(args.debug() && args.algorithm() == "forward" && qry_seqs.n_seqs() > 1)
      print_forward_scores(numerator, single_scores, bayes_factor);

    if(args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {   // if we were asked for more events than we found
      if(result.events_.size() > 0)
        cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
        assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args, result.events_, qry_seqs, bayes_factor, errors);
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, Sequences &seqs, double total_score, string errors) {
  if(args.algorithm() == "viterbi") {
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
          << "," << errors
          << endl;
    }
  } else {
    ofs
        << seqs.name_str(":")
        << "," << total_score
        << "," << errors
        << endl;
  }
}
// ----------------------------------------------------------------------------------------
void print_forward_scores(double numerator, vector<double> single_scores, double bayes_factor) {
  printf("   %8.3f = ", bayes_factor);
  printf("%2s %8.2f", "", numerator);
  for(auto & score : single_scores)
    printf(" - %8.2f", score);
  printf("\n");
}
// ----------------------------------------------------------------------------------------
int minimal_hamming_distance(Sequences &seqs_a, Sequences &seqs_b) {
  // Minimal hamming distance between any sequence in <seqs_a> and any sequence in <seqs_b>
  // NOTE for now, we require sequences are the same length (otherwise we have to deal with alignming them which is what we would call a CAN OF WORMS.
  assert(seqs_a.n_seqs() > 0 && seqs_b.n_seqs() > 0);
  int min_distance(seqs_a[0].size());
  for(size_t is=0; is<seqs_a.n_seqs(); ++is) {
    for(size_t js=0; js<seqs_b.n_seqs(); ++js) {
      Sequence seq_a = seqs_a[is];
      Sequence seq_b = seqs_b[js];
      assert(seq_a.size() == seq_b.size());
      int distance(0);
      for(size_t ic=0; ic<seq_a.size(); ++ic) {
	if(seq_a[ic] != seq_b[ic])
	  ++distance;
      }
      if(distance < min_distance)
	min_distance = distance;
    }
  }
  return min_distance;
}
// ----------------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------------
void hierarch_agglom(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs) {  // reminder: <qry_seq_list> is a list of lists
  map<string, Sequences> info;  // convert the vector to a map
  map<string, vector<string> > only_genes;
  map<string, KBounds> kbinfo;
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    info[qry_seq_list[iqry].name_str(":")] = qry_seq_list[iqry];  // NOTE this is probably kind of inefficient to remake the Sequences all the time
    only_genes[qry_seq_list[iqry].name_str(":")] = args.str_lists_["only_genes"][iqry];

    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kb(kmin, kmax);
    kbinfo[qry_seq_list[iqry].name_str(":")] = kb;
  }
  // for(int i=0; i<3; ++i) {
  while(info.size() > 0) {
    cout << "===================================" << endl;
    glomerate(hmms, gl, qry_seq_list, args, ofs, info, kbinfo, only_genes);
  }
}
// ----------------------------------------------------------------------------------------
void glomerate(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args &args, ofstream &ofs, map<string, Sequences> &info,
				 map<string, KBounds> &kbinfo, map<string, vector<string> > &only_genes) {  // reminder: <qry_seq_list> is a list of lists
  double max_log_prob(-INFINITY);
  pair<string, string> max_pair; // pair of clusters with largest log prob (i.e. the ratio of their prob together to their prob apart is largest)
  // Sequences max_seqs;
  vector<string> max_only_genes;
  KBounds max_kbounds;

  set<string> finished;
  for(auto &kv_a : info) {
    for(auto &kv_b : info) {
      if(kv_a.first == kv_b.first) continue;
      vector<string> names{kv_a.first, kv_b.first};
      sort(names.begin(), names.end());
      string namekey(names[0] + "---" + names[1]);
      if(finished.count(namekey)) {
	continue;
      } else {
	finished.insert(namekey);
      }
  // for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
  //   // NOTE this is really, really, really wastefull to recalculate the individual-cluster (denomiator) log probs every time through
  //   for(size_t jqry = iqry+1; jqry < qry_seq_list.size(); jqry++) {
      // if(iqry == jqry) continue;
      if(args.debug()) cout << "  " << kv_a.first << " " << kv_b.first << "  ---------" << endl;
      // Sequences qry_seqs_a(qry_seq_list[iqry]), qry_seqs_b(qry_seq_list[jqry]);
      Sequences qry_seqs_a(kv_a.second), qry_seqs_b(kv_b.second);
      double hamming_fraction = float(minimal_hamming_distance(qry_seqs_a, qry_seqs_b)) / qry_seqs_a[0].size();  // minimal_hamming_distance() will fail if the seqs aren't all the same length
      cout << hamming_fraction << " " << args.hamming_fraction_cutoff() << endl;
      if(hamming_fraction > args.hamming_fraction_cutoff()) {  // if all sequences in a are too far away from all sequences in b
	if(args.debug()) cout << "    skipping hamming: " << hamming_fraction << endl;
	continue;
      }

      Sequences union_seqs;
      for(size_t is=0; is<qry_seqs_a.n_seqs(); ++is)
	union_seqs.AddSeq(qry_seqs_a[is]);
      for(size_t is=0; is<qry_seqs_b.n_seqs(); ++is)
	union_seqs.AddSeq(qry_seqs_b[is]);
      vector<string> union_only_genes = only_genes[kv_a.first];
      for(auto &g : only_genes[kv_b.first])  // NOTE this will add duplicates (that's no big deal, though)
	union_only_genes.push_back(g);

      JobHolder jh(gl, hmms, args.algorithm(), union_only_genes);
      jh.SetDebug(args.debug());
      jh.SetChunkCache(args.chunk_cache());
      jh.SetNBestEvents(args.n_best_events());

      KBounds union_kbounds = kbinfo[kv_a.first].LogicalOr(kbinfo[kv_b.first]);
      Result numer_result(union_kbounds), result_a(union_kbounds), result_b(union_kbounds);
      double numer_logprob(-INFINITY);  // numerator in P(A,B,C,...) / (P(A)P(B)P(C)...)
      double logprob_a(-INFINITY), logprob_b(-INFINITY);
      double bayes_factor(-INFINITY); // or maybe likelihood ratio?
      bool stop(false);
      string errors;
      do {
	errors = "";
	// clock_t run_start(clock());
	if(args.debug()) cout << "       ----" << endl;
	numer_result = jh.Run(union_seqs, union_kbounds);
	numer_logprob = numer_result.total_score();
	result_a = jh.Run(qry_seqs_a, union_kbounds);
	result_b = jh.Run(qry_seqs_b, union_kbounds);
	logprob_a = result_a.total_score();
	logprob_b = result_b.total_score();
	bayes_factor = numer_logprob - logprob_a - logprob_b;

	union_kbounds = numer_result.better_kbounds();
	union_kbounds = union_kbounds.LogicalOr(result_a.better_kbounds());
	union_kbounds = union_kbounds.LogicalOr(result_b.better_kbounds());

	stop = !numer_result.boundary_error() || numer_result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
	stop &= !result_a.boundary_error() || result_a.could_not_expand();
	stop &= !result_b.boundary_error() || result_b.could_not_expand();
	if(args.debug() && !stop)
	  cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
	// cout << "      time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;
	if(numer_result.boundary_error() || result_a.boundary_error() || result_b.boundary_error())
	  errors = "boundary";
      } while(!stop);

      if(args.debug())
	print_forward_scores(numer_logprob, {logprob_a, logprob_b}, bayes_factor);

      // StreamOutput(ofs, args, numer_result.events_, qry_seqs, bayes_factor, errors);
      if(bayes_factor > max_log_prob) {
	max_log_prob = bayes_factor;
	max_pair = pair<string, string>(kv_a.first, kv_b.first);
	// max_seqs = union_seqs;
	max_only_genes = union_only_genes;
	max_kbounds = union_kbounds;
      }
    }
  }

  // then merge the two best clusters
  if(max_log_prob == -INFINITY) {
    throw runtime_error("no more!\n");
  }
  // info[max_seqs.name_str(":")] = max_seqs;
  Sequences max_seqs;
  // info[max_seqs.name_str(":")] = Sequences();
  // for(size_t is=0; is<info[max_pair.first].n_seqs(); ++is)
  //   info[max_seqs.name_str(":")].AddSeq(info[max_pair.first][is]);
  // cout << "1---" << endl;
  // for(size_t is=0; is<info[max_pair.second].n_seqs(); ++is)
  //   info[max_seqs.name_str(":")].AddSeq(info[max_pair.second][is]);

  // NOTE malloc/free problems with the way of setting max_seqs above!!! should figure that out
  for(size_t is=0; is<info[max_pair.first].n_seqs(); ++is)
    max_seqs.AddSeq(info[max_pair.first][is]);
  for(size_t is=0; is<info[max_pair.second].n_seqs(); ++is)
    max_seqs.AddSeq(info[max_pair.second][is]);
  info[max_seqs.name_str(":")] = max_seqs;
  kbinfo[max_seqs.name_str(":")] = max_kbounds;
  only_genes[max_seqs.name_str(":")] = max_only_genes;

  info.erase(max_pair.first);
  info.erase(max_pair.second);
  kbinfo.erase(max_pair.first);
  kbinfo.erase(max_pair.second);
  only_genes.erase(max_pair.first);
  only_genes.erase(max_pair.second);

  if(args.debug()) {
    cout << "  clustering " << max_pair.first << " " << max_pair.second << endl;
    cout << "  to give" << endl;
    for(auto &kv : info) cout << kv.first << endl;
  }

}
