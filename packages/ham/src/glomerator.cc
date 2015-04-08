#include "glomerator.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track):  // reminder: <qry_seq_list> is a list of lists
  hmms_(hmms),
  gl_(gl),
  args_(args),
  sampler_(args.smc_particles(), SMC_HISTORY_NONE),
  moveset_(SMCInit, SMCMove)
{
  ReadCachedLogProbs(track);

  sampler_.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
  sampler_.SetMoveSet(moveset_);
  sampler_.Initialise();

  // convert input vector to maps
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key(qry_seq_list[iqry].name_str(":"));
    KSet kmin(args->integers_["k_v_min"][iqry], args->integers_["k_d_min"][iqry]);
    KSet kmax(args->integers_["k_v_max"][iqry], args->integers_["k_d_max"][iqry]);
    KBounds kb(kmin, kmax);
    info_[key] = qry_seq_list[iqry];  // NOTE this is probably kind of inefficient to remake the Sequences all the time
    only_genes_[key] = args->str_lists_["only_genes"][iqry];
    kbinfo_[key] = kb;
  }

  initial_partition_ = GetClusterList(info_);
  initial_logprob_ = LogProbOfPartition(initial_partition_);
  if(args_->debug())
    PrintPartition(initial_partition, "initial");

  if(args.smc_particles() == 1)  // no smc, so push back one manually
    clusterpaths_.push_back(ClusterPath(initial_partition_, initial_logprob_));
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  if(args_->debug()) cout << "   glomerating" << endl;

  if(args.smc_particles() == 1) {  // don't do smc
    assert(clusterpaths_.size() == 1);
    do {
      Merge(clusterpaths_[0]);
    } while(!clusterpaths_[0].finished_);
  } else {
    do {
      sampler_.Iterate();
    } while(!AllFinished());
  }

  // if(args_->debug()) {
  //   cout << "  -----------------" << endl;  
  //   PrintPartition(best_partition_, "best");
  // }

  WritePartitions();
  WriteCachedLogProbs();
}

// ----------------------------------------------------------------------------------------
smc::particle<ClusterPath> Glomerator::SMCInit(smc::rng *rgen) {
  return smc::particle<ClusterPath>(ClusterPath(initial_partition_, initial_logprob_), initial_logprob_);
}

// ----------------------------------------------------------------------------------------
void Glomerator::SMCMove(long time, smc::particle<ClusterPath> &ptl, smc::rng *rgen) {
  
}

// ----------------------------------------------------------------------------------------
bool Glomerator::AllFinished() {
  bool all_finished(true);
  for(int ip = 0; ip < clusterpaths_.size(); ++ip)
    all_finished &= clusterpaths_[ip].finished_;
  return all_finished;
}

// ----------------------------------------------------------------------------------------
void Glomerator::ReadCachedLogProbs(Track *track) {
  ifstream ifs(args_->cachefile());
  if(!ifs.is_open()) {  // this means we don't have any cached results to start with, but we'll write out what we have at the end of the run to this file
    // throw runtime_error("ERROR cache file (" + args_->cachefile() + ") d.n.e.\n");
    return;
  }
  string line;

  // check the header is right (no cached info)
  if(!getline(ifs, line)) {
    return;  // return for zero length file
  }
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  vector<string> headstrs(SplitString(line, ","));
  assert(headstrs[0].find("unique_ids") == 0);
  assert(headstrs[1].find("score") == 0);
  assert(headstrs[2].find("naive-seq") == 0);

  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> column_list = SplitString(line, ",");
    assert(column_list.size() == 3);
    string unique_ids(column_list[0]);
    double logprob(stof(column_list[1]));
    string naive_seq(column_list[2]);
    log_probs_[unique_ids] = logprob;
    if(naive_seq.size() > 0)
      naive_seqs_[unique_ids] = naive_seq;
  }
  cout << "      read " << log_probs_.size() << " cached results" << endl;
}

// ----------------------------------------------------------------------------------------
// return a vector consisting of the keys in <partinfo>
Partitition Glomerator::GetClusterList(map<string, Sequences> &partinfo) {
  Partition clusters;
  for(auto &kv : partinfo)
    clusters.insert(kv.first);
  return clusters;
}

// ----------------------------------------------------------------------------------------
void Glomerator::GetSoloLogProb(string key) {
  // NOTE the only reason to have this separate from GetLogProb is the only_genes stuff
  DPHandler dph(args_, gl_, hmms_, only_genes_[key]);
  GetLogProb(dph, key, info_[key], kbinfo_[key]);
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(Partition &partition) {
  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  for(auto &key : partition) {
    if(log_probs_.count(key) == 0)
      GetSoloLogProb(key);  // should only happen for the initial partition
      // throw runtime_error("ERROR couldn't find key " + key + " in cached log probs\n");
    total_log_prob = AddWithMinusInfinities(total_log_prob, log_probs_[key]);
  }
  return total_log_prob;
}

// ----------------------------------------------------------------------------------------
void Glomerator::PrintPartition(Partition &partition, string extrastr) {
  const char *extra_cstr(extrastr.c_str());  // dammit I shouldn't need this line
  printf("    %-8.2f %s partition\n", LogProbOfPartition(partition), extra_cstr);
  for(auto &key : partition)
    cout << "          " << key << endl;
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteCachedLogProbs() {
  ofstream log_prob_ofs(args_->cachefile());
  if(!log_prob_ofs.is_open())
    throw runtime_error("ERROR cache file (" + args_->cachefile() + ") d.n.e.\n");
  log_prob_ofs << "unique_ids,score,naive-seq,errors" << endl;
  for(auto &kv : log_probs_) {
    if(naive_seqs_.count(kv.first) == 0)  // we end up having the log prob for a key, but not the naive sequence, if we calculated the merged prob for two clusters but didn't end up merging them
      continue;
    log_prob_ofs << kv.first << "," << kv.second << "," << naive_seqs_[kv.first] << "," << errors_[kv.first] << endl;  // NOTE if there's no errors, it just prints the empty string, which is fine
  }
  log_prob_ofs.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WritePartitions() {
  ofs_.open(args_->outfile());
  ofs_ << "path_index,partition,score" << endl;
  for(unsigned ipath=0; ipath<clusterpaths_.size(); ++ipath) {
    ClusterPath cp(clusterpaths_[ipath]);
    ofs_ << ipath << ",";
    for(unsigned ipart=0; ipart<cp.partitions_.size(); ++ipart) {
      for(auto &cluster : cp.partitions_[ipart]) {
	if(ic > 0)
	  ofs_ << ";";
	ofs_ << cluster;
      }
      ofs_ << "," << cp.log_probs_[ipart] << "\n";
    }
  }
  ofs_.close();
}

// ----------------------------------------------------------------------------------------
int Glomerator::HammingDistance(string seq_a, string seq_b) {
  if(seq_a.size() != seq_b.size())
    throw runtime_error("ERROR sequences different length in Glomerator::HammingDistance (" + seq_a + "," + seq_b + ")\n");
  int distance(0);
  for(size_t ic=0; ic<seq_a.size(); ++ic) {
    if(seq_a[ic] != seq_b[ic])
      ++distance;
  }
  return distance;
}

// ----------------------------------------------------------------------------------------
int Glomerator::HammingDistance(Sequence seq_a, Sequence seq_b) {
  assert(seq_a.size() == seq_b.size());
  int distance(0);
  for(size_t ic=0; ic<seq_a.size(); ++ic) {
    if(seq_a[ic] != seq_b[ic])
      ++distance;
  }
  return distance;
}

// ----------------------------------------------------------------------------------------
int Glomerator::NaiveHammingDistance(string key_a, string key_b) {
  GetNaiveSeq(key_a);
  GetNaiveSeq(key_b);
  return HammingDistance(naive_seqs_[key_a], naive_seqs_[key_b]);
}
// // ----------------------------------------------------------------------------------------
// int Glomerator::MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b) {
//   // Minimal hamming distance between any sequence in <seqs_a> and any sequence in <seqs_b>
//   // NOTE for now, we require sequences are the same length (otherwise we have to deal with alignming them which is what we would call a CAN OF WORMS.
//   assert(seqs_a.n_seqs() > 0 && seqs_b.n_seqs() > 0);
//   int min_distance(seqs_a[0].size());
//   for(size_t is=0; is<seqs_a.n_seqs(); ++is) {
//     for(size_t js=0; js<seqs_b.n_seqs(); ++js) {
//       int distance = HammingDistance(seqs_a[is], seqs_b[js]);
//       if(distance < min_distance)
// 	min_distance = distance;
//     }
//   }
//   return min_distance;
// }

// ----------------------------------------------------------------------------------------
// add log prob for <key>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetNaiveSeq(string key) {
  if(naive_seqs_.count(key))  // already did it
    return;

  DPHandler dph(args_, gl_, hmms_, only_genes_[key]);
  Result result(kbinfo_[key]);
  bool stop(false);
  do {
    result = dph.Run("viterbi", info_[key], kbinfo_[key]);
    kbinfo_[key] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.events_.size() < 1)
    throw runtime_error("ERROR no events for" + key + "\n");
  naive_seqs_[key] = result.events_[0].naive_seq_;  // NOTE it might be a bit wasteful to store these digitized? Or maybe it's better...
  if(result.boundary_error())
    errors_[key] = errors_[key] + ":boundary";
}

// ----------------------------------------------------------------------------------------
// add log prob for <name>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetLogProb(DPHandler &dph, string name, Sequences &seqs, KBounds &kbounds) {
  // NOTE that when this imporves the kbounds, that info doesn't get propagated to <kbinfo_>
  if(log_probs_.count(name))  // already did it
    return;
    
  Result result(kbounds);
  bool stop(false);
  do {
    result = dph.Run("forward", seqs, kbounds);
    kbounds = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  log_probs_[name] = result.total_score();
  if(result.boundary_error())
    errors_[name] = errors_[name] + ":boundary";
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinNames(string name1, string name2) {
  vector<string> names{name1, name2};
  sort(names.begin(), names.end());  // NOTE this doesn't sort *within* name1 or name2 when they're already comprised of several uids
  return names[0] + ":" + names[1];
}

// ----------------------------------------------------------------------------------------
// perform one merge step, i.e. find the two "nearest" clusters and merge 'em
void Glomerator::Merge(ClusterPath &path) {
  if(path.finished_)  // already finished this <path>, but we're still iterating 'cause some of the other paths aren't finished
    return;
  double max_log_prob(-INFINITY);
  pair<string, string> max_pair; // pair of clusters with largest log prob (i.e. the ratio of their prob together to their prob apart is largest)
  vector<string> max_only_genes;
  KBounds max_kbounds;

  int n_total_pairs(0), n_skipped_hamming(0);

  set<string> already_done;  // keeps track of which  a-b pairs we've already done
  for(auto &key_a : path.CurrentPartition()) {  // note that c++ maps are ordered
    for(auto &key_b : path.CurrentPartition()) {
      if(key_a == key_b) continue;
      string bothnamestr = JoinNames(key_a, key_b);
      if(already_done.count(bothnamestr))  // already did this pair
  	continue;
      else
  	already_done.insert(bothnamestr);
      Sequences a_seqs(info_[key_a]), b_seqs(info[key_b]);
  // for(unsigned i_a=0; i_a<path.CurrentPartition().size(); ++i_a) {
  //   for(unsigned i_b=i_a+1; i_b<path.CurrentPartition().size(); ++i_b) {
  //     string key_a(path.CurrentPartition()[i_a]), key_b(path.CurrentPartition()[i_b]);
  //     string bothnamestr = JoinNames(key_a, key_b);

      ++n_total_pairs;

      // NOTE it might help to also cache hamming fractions
      double hamming_fraction = float(NaiveHammingDistance(key_a, key_b)) / a_seqs[0].size();  // hamming distance fcn will fail if the seqs aren't the same length
      if(hamming_fraction > args_->hamming_fraction_cutoff()) {
	++n_skipped_hamming;
	continue;
      }

      // NOTE it's a bit wasteful to do this if we already have all three of 'em cached
      Sequences ab_seqs(a_seqs.Union(b_seqs));
      vector<string> ab_only_genes = only_genes_[key_a];
      for(auto &g : only_genes_[key_b])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
	ab_only_genes.push_back(g);
      KBounds ab_kbounds = kbinfo_[key_a].LogicalOr(kbinfo_[key_b]);

      DPHandler dph(args_, gl_, hmms_, ab_only_genes);  // NOTE it's an ok approximation to compare log probs for sequence sets that were run with different kbounds, but (I'm pretty sure) we do need to run them with the same set of genes. EDIT hm, well, maybe not. Anywa, it's ok for now

      // NOTE the error from using the single kbounds rather than the OR seems to be around a part in a thousand or less
      GetLogProb(dph, key_a, a_seqs, kbinfo_[key_a]);
      GetLogProb(dph, key_b, b_seqs, kbinfo_[key_b]);
      GetLogProb(dph, bothnamestr, ab_seqs, ab_kbounds);

      double bayes_factor(log_probs_[bothnamestr] - log_probs_[key_a] - log_probs_[key_b]);  // REMINDER a, b not necessarily same order as names[0], names[1]
      if(args_->debug()) {
	printf("       %8.3f = ", bayes_factor);
	printf("%2s %8.2f", "", log_probs_[bothnamestr]);
	printf(" - %8.2f - %8.2f", log_probs_[key_a], log_probs_[key_b]);
	printf("\n");
      }

      if(bayes_factor > max_log_prob) {
	max_log_prob = bayes_factor;
	max_pair = pair<string, string>(key_a, key_b);  // REMINDER not always same order as names[0], names[1]
	max_only_genes = ab_only_genes;
	max_kbounds = ab_kbounds;
      }
    }
  }

  // if <info_> only has one cluster, if hamming is too large between all remaining clusters/remaining bayes factors are -INFINITY
  if(max_log_prob == -INFINITY) {
    path.finished_ = true;
    return;
  }

  // then merge the two best clusters
  vector<string> max_names{max_pair.first, max_pair.second};
  sort(max_names.begin(), max_names.end());
  Sequences max_seqs(info_[max_names[0]].Union(info_[max_names[1]]));  // NOTE this will give the ordering {<first seqs>, <second seqs>}, which should be the same as in <max_name_str>. But I don't think it'd hurt anything if the sequences and names were in a different order
  string max_name_str = JoinNames(max_names[0], max_names[1]);  // NOTE the names[i] are not sorted *within* themselves, but <names> itself is sorted. This is ok, because we will never again encounter these sequences separately
  info_[max_name_str] = max_seqs;
  kbinfo_[max_name_str] = max_kbounds;
  only_genes_[max_name_str] = max_only_genes;

  GetNaiveSeq(max_name_str);

  Partition new_partition(path.CurrentPartition());
  new_partition.erase(max_pair.first);
  new_partition.erase(max_pair.second);
  new_partition.insert(max_name_str);

  path.AddPartition(new_partition, LogProbOfPartition(new_partition));

  if(args_->debug()) {
    printf("          hamming skipped %d / %d\n", n_skipped_hamming, n_total_pairs);
    printf("       merged %-8.2f %s and %s\n", max_log_prob, max_pair.first.c_str(), max_pair.second.c_str());
  }
}

}
