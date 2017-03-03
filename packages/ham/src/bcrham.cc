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
  srand(args.random_seed());

  // init some infrastructure
  vector<string> characters {"A", "C", "G", "T"};
  Track track("NUKES", characters, args.ambig_base());
  GermLines gl(args.datadir(), args.locus());
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

  printf("        time: bcrham %.1f\n", ((clock() - run_start) / (double)CLOCKS_PER_SEC));
  return 0;
}

// ----------------------------------------------------------------------------------------
// read input sequences from file and return as vector of sequences
vector<vector<Sequence> > GetSeqs(Args &args, Track *trk) {
  vector<vector<Sequence> > all_seqs;
  assert(args.str_lists_["names"].size() == args.str_lists_["seqs"].size());
  for(size_t iqry = 0; iqry < args.str_lists_["names"].size(); ++iqry) { // loop over queries, where each query can be composed of one, two, or k sequences
    vector<Sequence> seqs;
    assert(args.str_lists_["names"][iqry].size() == args.str_lists_["seqs"][iqry].size());
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

  // write csv output headers
  ofstream ofs;
  ofs.open(args.outfile());
  if(!ofs.is_open())
    throw runtime_error("ERROR --outfile (" + args.outfile() + ") d.n.e.\n");
  StreamHeader(ofs, args.algorithm());

  int n_vtb_calculated(0), n_fwd_calculated(0);

  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    if(args.debug() > 1) cout << "  ---------" << endl;
    KSet kmin(args.integers_["k_v_min"][iqry], args.integers_["k_d_min"][iqry]);
    KSet kmax(args.integers_["k_v_max"][iqry], args.integers_["k_d_max"][iqry]);
    KBounds kbounds(kmin, kmax);
    vector<Sequence> qry_seqs(qry_seq_list[iqry]);

    DPHandler dph(args.algorithm(), &args, gl, hmms);
    Result result = dph.Run(qry_seqs, kbounds, args.str_lists_["only_genes"][iqry], args.floats_["mut_freq"][iqry]);
    // if(FishyMultiSeqAnnotation(qry_seqs.size(), result.best_event()))
    //   dph.HandleFishyAnnotations(result, qry_seqs, kbounds, args.str_lists_["only_genes"][iqry], args.floats_["mut_freq"][iqry]);

    if(args.debug() > 1) cout << "       ----" << endl;

    if(result.no_path_)
      StreamErrorput(ofs, args.algorithm(), qry_seqs, "no_path");
    else if(args.algorithm() == "viterbi")
      StreamViterbiOutput(ofs, result.best_event(), qry_seqs, "");
    else if(args.algorithm() == "forward")
      StreamForwardOutput(ofs, qry_seqs, result.total_score(), "");
    else
      assert(0);

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

