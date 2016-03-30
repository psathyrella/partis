#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <map>
#include <iomanip>
#include <ctime>
#include <fstream>
#include <cfenv>

#include "dphandler.h"
#include "bcrutils.h"
#include "text.h"
#include "args.h"
#include "glomerator.h"
#include "tclap/CmdLine.h"
// #include "smctc.hh"

using namespace TCLAP;
using namespace ham;
using namespace std;

// ----------------------------------------------------------------------------------------
vector<vector<Sequence> > GetSeqs(Args &args, Track *trk);
void run_algorithm(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args &args);

// ----------------------------------------------------------------------------------------
int main(int argc, const char * argv[]) {
  clock_t run_start(clock());
  Args args(argc, argv);
  if(args.smc_particles() > 1)
    assert(0);  // see commented code below
  srand(args.random_seed());

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track track("NUKES", characters, args.ambig_base());
  GermLines gl(args.datadir());
  HMMHolder hmms(args.hmmdir(), gl, &track);
  vector<vector<Sequence> > qry_seq_list(GetSeqs(args, &track));

  if(args.cache_naive_seqs()) {
    Glomerator glom(hmms, gl, qry_seq_list, &args, &track);
    glom.CacheNaiveSeqs();
  } else if(args.partition()) {  // NOTE this is kind of hackey -- there's some code duplication between Glomerator and the loop below... but only a little, and they're doing fairly different things, so screw it for the time being
    Glomerator glom(hmms, gl, qry_seq_list, &args, &track);
    glom.Cluster();
  } else {
    run_algorithm(hmms, gl, qry_seq_list, args);
  }

  printf("        bcrham time: %.1f\n", ((clock() - run_start) / (double)CLOCKS_PER_SEC));
  return 0;
}

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<vector<Sequence> > GetSeqs(Args &args, Track *trk) {
  vector<vector<Sequence> > all_seqs;
  assert(args.str_lists_["names"].size() == args.str_lists_["seqs"].size());
  assert(args.str_lists_["names"].size() == args.float_lists_["mute_freqs"].size());
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    vector<Sequence> seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
    assert(args.str_lists_["names"][iqry].size() == args.float_lists_["mute_freqs"][iqry].size());
    for(size_t iseq = 0; iseq < args.str_lists_["names"][iqry].size(); ++iseq) { // loop over each sequence in that query
      Sequence sq(trk, args.str_lists_["names"][iqry][iseq], args.str_lists_["seqs"][iqry][iseq]);
      seqs.push_back(sq);
    }
    all_seqs.push_back(seqs);
  }
  assert(all_seqs.size() == args.str_lists_["names"].size());
  return all_seqs;
}

// ----------------------------------------------------------------------------------------
void run_algorithm(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args &args) {

  DPHandler dph(args.algorithm(), &args, gl, hmms);

  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  if(!ofs.is_open())
    throw runtime_error("ERROR --outfile (" + args.outfile() + ") d.n.e.\n");
  StreamHeader(ofs, args.algorithm());

  int n_vtb_calculated(0), n_fwd_calculated(0);

  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    if(args.debug()) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kbounds(kmin, kmax);
    vector<Sequence> qry_seqs(qry_seq_list[iqry]);
    vector<double> mute_freqs(args.float_lists_["mute_freqs"][iqry]);
    double mean_mute_freq(avgVector(mute_freqs));

    vector<KBounds> kbvector(qry_seqs.size(), kbounds);

    Result result(kbounds);
    double logprob(-INFINITY);
    bool stop(false);
    string errors;
    do {
      errors = "";
      if(args.debug()) cout << "       ----" << endl;
      result = dph.Run(qry_seqs, kbounds, args.str_lists_["only_genes"][iqry], mean_mute_freq);
      logprob = result.total_score();
      kbounds = result.better_kbounds();
      stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
      if(args.debug() && !stop)
        cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
      if(result.boundary_error())
        errors = "boundary";
    } while(!stop);

    if(args.algorithm() == "viterbi" && size_t(args.n_best_events()) > result.events_.size()) {   // if we were asked for more events than we found
      if(result.events_.size() > 0)
        cout << "WARNING asked for " << args.n_best_events() << " events but only found " << result.events_.size() << endl;
      else
        assert(result.no_path_);  // if there's some *other* way we can end up with no events, I want to know about it
    }
    StreamOutput(ofs, args.algorithm(), min(size_t(args.n_best_events()), result.events_.size()), result.events_, qry_seqs, logprob, errors);
    if(args.algorithm() == "viterbi")
      ++n_vtb_calculated;
    else if(args.algorithm() == "forward")
      ++n_fwd_calculated;
  }
  printf("        calcd:   vtb %-4d  fwd %-4d\n", n_vtb_calculated, n_fwd_calculated);
  ofs.close();
}


// Glomerator *stupid_global_glom;  // I *(#*$$!*ING HATE GLOBALS

// } else {
//   if(args.debug()) cout << "   glomerating with smc" << endl;
//   stupid_global_glom = &glom;
//   smc::sampler<ClusterPath> smp(args.smc_particles(), SMC_HISTORY_NONE);
//   smc::moveset<ClusterPath> mvs(SMCInit, SMCMove);
//   smp.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
//   smp.SetMoveSet(mvs);
//   smp.Initialise();

//   bool finished(true);
//   do {
//   	smp.Iterate();
//   	finished = true;
//   	for(int ip=0; ip<args.smc_particles(); ++ip)
//   	  finished &= smp.GetParticleValue(ip).finished_;
//   } while(!finished);

//   vector<ClusterPath> paths;
//   for(int ip=0; ip<args.smc_particles(); ++ip) {
//   	paths.push_back(smp.GetParticleValue(ip));
//   }
//   glom.WritePartitions(paths);
// }
// // ----------------------------------------------------------------------------------------
// smc::particle<ClusterPath> SMCInit(smc::rng *rgen) {
//   int initial_path_index(-1);
//   throw runtime_error("I think this is ok, but I just changed the default path index in clusterpath from -1 to 0, and it should be double checked that this doesn't screw up anything in this fcn.\n");
//   double logweight;
//   Partition initial_partition(stupid_global_glom->GetAnInitialPartition(initial_path_index, logweight));  // get the next initial partition (and increment the counter)
//   double logprob = stupid_global_glom->LogProbOfPartition(initial_partition);
//   ClusterPath thecp(initial_partition, logprob, logweight);
//   thecp.initial_path_index_ = initial_path_index;
//   return smc::particle<ClusterPath>(thecp, thecp.CurrentLogWeight());  //logprob);
// }

// // ----------------------------------------------------------------------------------------
// void SMCMove(long time, smc::particle<ClusterPath> &ptl, smc::rng *rgen) {
//   stupid_global_glom->Merge(ptl.GetValuePointer(), rgen);  // ...but I couldn't figure out a good way to do this without a global
//   // ptl.SetLogWeight(ptl.GetValuePointer()->CurrentLogProb());
//   ptl.SetLogWeight(ptl.GetValuePointer()->CurrentLogWeight());
// }

