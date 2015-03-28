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
#include "bcrutils.h"
#include "text.h"
#include "args.h"
#include "glomerator.h"
#include "tclap/CmdLine.h"
using namespace TCLAP;
using namespace ham;
using namespace std;

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<Sequences> GetSeqs(Args &args, Track *trk) {
  vector<Sequences> all_seqs;
  set<string> all_names;  // keeps track which sequences we've added to make sure we weren't passed duplicates
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    Sequences seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(trk, args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq]);

      if(args.partition() && all_names.count(sq.name()))  // Not reall sure if we need this check, but I don't feel like thinking about it right now. In any case, we want to be able to add the same sequqence twice for run_algorithm
	throw runtime_error("ERROR tried to add sequence with name " + sq.name() + " twice in bcrham::GetSeqs");
      else
	all_names.insert(sq.name());

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

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  srand(time(NULL));
  Args args(argc, argv);
  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  if(!ofs.is_open())
    throw runtime_error("ERROR --outfile (" + args.outfile() + ") d.n.e.\n");
  if(args.algorithm() == "viterbi")
    ofs << "unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion,v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,score,seqs,errors" << endl;
  else if(args.algorithm() == "forward")
    ofs << "unique_ids,score,errors" << endl;

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track trk("NUKES", characters);
  vector<Sequences> qry_seq_list(GetSeqs(args, &trk));
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl);
  // hmms.CacheAll();

  if(args.partition()) {  // NOTE this is kind of hackey -- there's some code duplication between Glomerator and the loop below... but only a little, and they're doing fairly different things, so screw it for the time being
    Glomerator glom(hmms, gl, qry_seq_list, &args, &trk);
    glom.Cluster();
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
