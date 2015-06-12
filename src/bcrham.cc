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

#include "dphandler.h"
#include "bcrutils.h"
#include "text.h"
#include "args.h"
#include "glomerator.h"
#include "tclap/CmdLine.h"
#include "smctc.hh"

using namespace TCLAP;
using namespace ham;
using namespace std;

Glomerator *stupid_global_glom;  // I *(#*$$!*ING HATE GLOBALS

// ----------------------------------------------------------------------------------------
smc::particle<ClusterPath> SMCInit(smc::rng *rgen) {
  int initial_path_index(-1);
  double logweight;
  Partition initial_partition(stupid_global_glom->GetAnInitialPartition(initial_path_index, logweight));  // get the next initial partition (and increment the counter)
  double logprob = stupid_global_glom->LogProbOfPartition(initial_partition);
  ClusterPath thecp(initial_partition, logprob, logweight);
  thecp.initial_path_index_ = initial_path_index;
  return smc::particle<ClusterPath>(thecp, thecp.CurrentLogWeight());  //logprob);
}

// ----------------------------------------------------------------------------------------
void SMCMove(long time, smc::particle<ClusterPath> &ptl, smc::rng *rgen) {
  stupid_global_glom->Merge(ptl.GetValuePointer(), rgen);  // ...but I couldn't figure out a good way to do this without a global
  // ptl.SetLogWeight(ptl.GetValuePointer()->CurrentLogProb());
  ptl.SetLogWeight(ptl.GetValuePointer()->CurrentLogWeight());
}

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<vector<Sequence> > GetSeqs(Args &args, Track *trk) {
  vector<vector<Sequence> > all_seqs;
  assert(args.str_lists_["names"].size() == args.str_lists_["seqs"].size());
  assert(args.str_lists_["names"].size() == args.float_lists_["mute_freqs"].size());
  assert(args.str_lists_["names"].size() == args.int_lists_["cyst_positions"].size());
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    vector<Sequence> seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    assert(args.str_lists_["names"][iqry].size() == args.float_lists_["mute_freqs"][iqry].size());
    assert(args.str_lists_["names"][iqry].size() == args.int_lists_["cyst_positions"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(trk, args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq], args.int_lists_["cyst_positions"][iqry][iseq]);
      seqs.push_back(sq);
    }
    all_seqs.push_back(seqs);
  }
  assert(all_seqs.size() == args.str_lists_["names"].size());
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, vector<Sequence> &seqs, double total_score, string errors);
void print_forward_scores(double numerator, vector<double> single_scores, double lratio);

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
    ofs << "nth_best,unique_ids,v_gene,d_gene,j_gene,fv_insertion,vd_insertion,dj_insertion,jf_insertion,v_5p_del,v_3p_del,d_5p_del,d_3p_del,j_5p_del,j_3p_del,score,seqs,errors" << endl;
  else if(args.algorithm() == "forward")
    ofs << "unique_ids,score,errors" << endl;

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track track("NUKES", characters, "N");
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl, &track);
  // Track *track(hmms.track());
  vector<vector<Sequence> > qry_seq_list(GetSeqs(args, &track));
  // hmms.CacheAll();

  if(args.partition()) {  // NOTE this is kind of hackey -- there's some code duplication between Glomerator and the loop below... but only a little, and they're doing fairly different things, so screw it for the time being
    Glomerator glom(hmms, gl, qry_seq_list, &args, &track);
    if(args.smc_particles() == 1) {
      glom.Cluster();
    } else {
      if(args.debug()) cout << "   glomerating with smc" << endl;
      stupid_global_glom = &glom;
      smc::sampler<ClusterPath> smp(args.smc_particles(), SMC_HISTORY_NONE);
      smc::moveset<ClusterPath> mvs(SMCInit, SMCMove);
      smp.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
      smp.SetMoveSet(mvs);
      smp.Initialise();

      bool finished(true);
      do {
	smp.Iterate();
	finished = true;
	for(int ip=0; ip<args.smc_particles(); ++ip)
	  finished &= smp.GetParticleValue(ip).finished_;
      } while(!finished);

      vector<ClusterPath> paths;
      for(int ip=0; ip<args.smc_particles(); ++ip) {
	paths.push_back(smp.GetParticleValue(ip));
      }
      glom.WritePartitions(paths);
    }
    ofs.close();
    return 0;
  }

  DPHandler dph(args.algorithm(), &args, gl, hmms);

  // not partitioning
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    if(args.debug()) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kbounds(kmin, kmax);
    vector<Sequence> qry_seqs(qry_seq_list[iqry]);
    vector<double> mute_freqs(args.float_lists_["mute_freqs"][iqry]);
    double mean_mute_freq(avgVector(mute_freqs));

    vector<KBounds> kbvector(qry_seqs.size(), kbounds);
    if(args.truncate_seqs())
      TruncateSeqs(qry_seqs, kbvector);

    Result result(kbounds);
    vector<Result> denom_results(qry_seqs.size(), result);  // only used for forward if n_seqs > 1
    double numerator(-INFINITY);  // numerator in P(A,B,C,...) / (P(A)P(B)P(C)...)
    double lratio(-INFINITY); // final result
    vector<double> single_scores(qry_seqs.size(), -INFINITY);  // NOTE log probs, not scores, but I haven't managed to finish switching over to the new terminology
    bool stop(false);
    string errors;
    do {
      errors = "";
      // clock_t run_start(clock());
      if(args.debug()) cout << "       ----" << endl;
      result = dph.Run(qry_seqs, kbounds, args.str_lists_["only_genes"][iqry], mean_mute_freq);
      numerator = result.total_score();
      lratio = numerator;
      if(args.algorithm() == "forward" && qry_seqs.size() > 1) {  // calculate factors for denominator
        for(size_t iseq = 0; iseq < qry_seqs.size(); ++iseq) {
          denom_results[iseq] = dph.Run(qry_seqs[iseq], kbounds, args.str_lists_["only_genes"][iqry], mean_mute_freq);  // result for a single sequence  TODO hm, wait, should this be the individual mute freqs?
          single_scores[iseq] = denom_results[iseq].total_score();
          lratio -= single_scores[iseq];
        }
      }

      kbounds = result.better_kbounds();
      for(auto & res : denom_results)
        kbounds = kbounds.LogicalOr(res.better_kbounds());

      stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
      for(auto & res : denom_results)
        stop &= !res.boundary_error() || res.could_not_expand();
      if(args.debug() && !stop)
        cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
      // cout << "      time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << endl;
      if(result.boundary_error())
        errors = "boundary";
      for(auto & res : denom_results) // NOTE <errors> will still just have one "boundary" in it even if multiple results had boundary errors
        if(res.boundary_error())
          errors = "boundary";
    } while(!stop);

    if(args.debug() && args.algorithm() == "forward" && qry_seqs.size() > 1)
      print_forward_scores(numerator, single_scores, lratio);

    if(args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {   // if we were asked for more events than we found
      if(result.events_.size() > 0)
        cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
        assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args, result.events_, qry_seqs, lratio, errors);
  }

  ofs.close();
  return 0;
}

// ----------------------------------------------------------------------------------------
void StreamOutput(ofstream &ofs, Args &args, vector<RecoEvent> &events, vector<Sequence> &seqs, double total_score, string errors) {
  if(args.algorithm() == "viterbi") {
    size_t n_max = min(size_t(args.n_best_events()), events.size());
    for(size_t ievt = 0; ievt < n_max; ++ievt) {
      RecoEvent *event = &events[ievt];
      string second_seq_name, second_seq;
      ofs  // be very, very careful to change this *and* the csv header above at the same time
	<< ievt
	<< "," << SeqNameStr(seqs, ":")
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
	<< "," << SeqStr(seqs, ":")
	<< "," << errors
	<< endl;
}
  } else {
    ofs
        << SeqNameStr(seqs, ":")
        << "," << total_score
        << "," << errors
        << endl;
  }
}
// ----------------------------------------------------------------------------------------
void print_forward_scores(double numerator, vector<double> single_scores, double lratio) {
  printf("   %8.3f = ", lratio);
  printf("%2s %8.2f", "", numerator);
  for(auto & score : single_scores)
    printf(" - %8.2f", score);
  printf("\n");
}
