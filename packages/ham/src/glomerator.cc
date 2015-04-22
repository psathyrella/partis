#include "glomerator.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track) :
  hmms_(hmms),
  gl_(gl),
  args_(args),
  i_initial_partition_(0)
{
  ReadCachedLogProbs(track);

  Partition tmp_partition;
  int last_ipath(0);
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key(qry_seq_list[iqry].name_str(":"));

    int ipath(args_->integers_["path_index"][iqry]);
    // if((int)initial_partitions_.size() < ipath + 1) {  // if we need to push back a new initial partition (i.e. this is a new path/particle
    if(last_ipath != ipath) {  // if we need to push back a new initial partition (i.e. this is a new path/particle) NOTE I'm assuming the the first one will have <path_index> zero
      initial_partitions_.push_back(tmp_partition);
      initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));
      tmp_partition.clear();
      last_ipath = ipath;
    }
      
    tmp_partition.insert(key);

    if(info_.count(key) > 0)  // already added this cluster to the maps (but it'll appear more than once if it's in multiple paths/particles), so we only need to add the cluster info to <initial_partitions_>
      continue;

    info_[key] = qry_seq_list[iqry];
    only_genes_[key] = args_->str_lists_["only_genes"][iqry];

    KSet kmin(args_->integers_["k_v_min"][iqry], args_->integers_["k_d_min"][iqry]);
    KSet kmax(args_->integers_["k_v_max"][iqry], args_->integers_["k_d_max"][iqry]);
    KBounds kb(kmin, kmax);
    kbinfo_[key] = kb;

    vector<double> mute_freqs(args_->float_lists_["mute_freqs"][iqry]);
    mute_freqs_[key] = avgVector(mute_freqs);
  }
  // add the last initial partition (i.e. for the last path/particle)
  assert(tmp_partition.size() > 0);
  initial_partitions_.push_back(tmp_partition);
  initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));

  // the *first* time, we only get one path/partition from partis, so we push back a bunch of copies
  if((int)initial_partitions_.size() == 1 && args_->smc_particles() > 1)  {
    for(int ip=1; ip<args_->smc_particles(); ++ip) {
      initial_partitions_.push_back(tmp_partition);
      initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));
    }
  }

  if((int)initial_partitions_.size() != args_->smc_particles())
    throw runtime_error("wrong number of initial partitions " + to_string(initial_partitions_.size()) + " (should be " + to_string(args_->smc_particles()) + ")");

  // initial_partition_ = GetPartitionFromMap(info_);
  // initial_logprob_ = LogProbOfPartition(initial_partition_);
  if(args_->debug())
    for(auto &part : initial_partitions_)
      PrintPartition(part, "initial");
}

// ----------------------------------------------------------------------------------------
Glomerator::~Glomerator() {
  WriteCachedLogProbs();
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  if(args_->debug()) cout << "   glomerating" << endl;

  assert((int)initial_partitions_.size() == 1);
  ClusterPath cp(initial_partitions_[0], LogProbOfPartition(initial_partitions_[0]));
  do {
    Merge(&cp);
  } while(!cp.finished_);

  vector<ClusterPath> paths{cp};
  WritePartitions(paths);
}

// ----------------------------------------------------------------------------------------
Partition Glomerator::GetAnInitialPartition() {
  assert(i_initial_partition_ < (int)initial_partitions_.size());
  return initial_partitions_.at(i_initial_partition_++);
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
    assert(column_list.size() == 3 || column_list.size() == 4);  // if written by bcrham it'll have an error column, while if written by partitiondriver it won't (under normal circumstances we wouldn't see bcrham cachefiles right here, but it can be useful for testing)
    string unique_ids(column_list[0]);
    double logprob(stof(column_list[1]));
    string naive_seq(column_list[2]);
    log_probs_[unique_ids] = logprob;
    if(naive_seq.size() > 0)
      naive_seqs_[unique_ids] = naive_seq;
  }
  // cout << "      read " << log_probs_.size() << " cached results" << endl;
}

// ----------------------------------------------------------------------------------------
// return a vector consisting of the keys in <partinfo>
Partition Glomerator::GetPartitionFromMap(map<string, Sequences> &partinfo) {
  Partition clusters;
  for(auto &kv : partinfo)
    clusters.insert(kv.first);
  return clusters;
}

// ----------------------------------------------------------------------------------------
void Glomerator::GetSoloLogProb(string key) {
  // NOTE the only reason to have this separate from GetLogProb is the only_genes stuff
  DPHandler dph(args_, gl_, hmms_, only_genes_[key]);
  GetLogProb(dph, key, info_[key], kbinfo_[key], mute_freqs_[key]);
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
void Glomerator::WritePartitions(vector<ClusterPath> &paths) {
  ofs_.open(args_->outfile());
  ofs_ << "path_index,partition,score" << endl;
  int ipath(0);
  for(auto &cp : paths) {
    for(unsigned ipart=0; ipart<cp.partitions().size(); ++ipart) {
      ofs_ << ipath << ",";
      int ic(0);
      for(auto &cluster : cp.partitions()[ipart]) {
	if(ic > 0)
	  ofs_ << ";";
	ofs_ << cluster;
	++ic;
      }
      ofs_ << "," << cp.logprobs()[ipart] << endl;
    }
    ++ipath;
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
    result = dph.Run("viterbi", info_[key], kbinfo_[key], mute_freqs_[key]);
    kbinfo_[key] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.events_.size() < 1)
    throw runtime_error("ERROR no events for " + key + "\n");
  naive_seqs_[key] = result.events_[0].naive_seq_;  // NOTE it might be a bit wasteful to store these digitized? Or maybe it's better...
  if(result.boundary_error())
    errors_[key] = errors_[key] + ":boundary";
}

// ----------------------------------------------------------------------------------------
// add log prob for <name>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetLogProb(DPHandler &dph, string name, Sequences &seqs, KBounds &kbounds, double mean_mute_freq) {
  // NOTE that when this imporves the kbounds, that info doesn't get propagated to <kbinfo_>
  if(log_probs_.count(name))  // already did it
    return;
    
  Result result(kbounds);
  bool stop(false);
  // dph.PrintHMMS();
  do {
    result = dph.Run("forward", seqs, kbounds, mean_mute_freq);
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
Query Glomerator::GetMergedQuery(string name_a, string name_b) {
  Query qmerged;
  qmerged.name_ = JoinNames(name_a, name_b);
  qmerged.seqs_ = info_[name_a].Union(info_[name_b]);
  qmerged.kbounds_ = kbinfo_[name_a].LogicalOr(kbinfo_[name_b]);
  qmerged.only_genes_ = only_genes_[name_a];
  for(auto &g : only_genes_[name_b])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
    qmerged.only_genes_.push_back(g);
  qmerged.mean_mute_freq_ = (info_[name_a].n_seqs()*mute_freqs_[name_a] + info_[name_b].n_seqs()*mute_freqs_[name_b]) / float(qmerged.seqs_.n_seqs());  // simple weighted average
  qmerged.parents_ = pair<string, string>(name_a, name_b);
  return qmerged;
}

// ----------------------------------------------------------------------------------------
Query *Glomerator::ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen) {
  // first leave log space and normalize. NOTE instead of normalizing, we could just use rng::Uniform(0., total)
  vector<double> ratios;  // NOTE *not* a probability: it's the ratio of the probability together to the probability apart
  double total(0.0);
  // cout << "choosing" << endl;
  for(auto &pr : potential_merges) {
    double likelihood_ratio = exp(pr.first);
    // cout << "  " << pr.second.name_ << " " << likelihood_ratio << endl;
    total += likelihood_ratio;
    ratios.push_back(likelihood_ratio);
  }
  // cout << "   normalizing" << endl;
  for(auto &ratio : ratios)
    ratio /= total;
  // for(auto &ratio : ratios)
    // cout << "  " << ratio << endl;

  // then choose one at random according to the probs
  double drawpoint = rgen->Uniform(0., 1.);
  double sum(0.0);
  // pair<unsigned, unsigned> icls(9999, 9999);
  // double chosennetprob(0.0);
  // cout << "  choosing with " << drawpoint << endl;
  for(size_t im=0; im<ratios.size(); ++im) {
    double ratio(ratios[im]);
    sum += ratio;
    if(sum > drawpoint) {
      return &potential_merges[im].second;
    }
  }

  throw runtime_error("ERROR fell through in Glomerator::ChooseRandomMerge");
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinNames(string name1, string name2) {
  vector<string> names{name1, name2};
  sort(names.begin(), names.end());  // NOTE this doesn't sort *within* name1 or name2 when they're already comprised of several uids
  return names[0] + ":" + names[1];
}

// ----------------------------------------------------------------------------------------
// perform one merge step, i.e. find the two "nearest" clusters and merge 'em
void Glomerator::Merge(ClusterPath *path, smc::rng *rgen) {
  if(path->finished_) { // already finished this <path>, but we're still iterating 'cause some of the other paths aren't finished
    return;
  }
  double max_log_prob(-INFINITY);
  // pair<string, string> max_pair; // pair of clusters with largest log prob (i.e. the ratio of their prob together to their prob apart is largest)
  int imax(-1);

  vector<pair<double, Query> > potential_merges;

  int n_total_pairs(0), n_skipped_hamming(0);

  set<string> already_done;  // keeps track of which  a-b pairs we've already done
  for(auto &key_a : path->CurrentPartition()) {  // note that c++ maps are ordered
    for(auto &key_b : path->CurrentPartition()) {  // also, note that CurrentPartition() returns a reference
      if(key_a == key_b) continue;
      Query qmerged(GetMergedQuery(key_a, key_b));
      if(already_done.count(qmerged.name_))  // already did this pair
  	continue;
      else
  	already_done.insert(qmerged.name_);

      Sequences *a_seqs(&info_[key_a]), *b_seqs(&info_[key_b]);
      ++n_total_pairs;

      // NOTE it might help to also cache hamming fractions
      double hamming_fraction = float(NaiveHammingDistance(key_a, key_b)) / (*a_seqs)[0].size();  // hamming distance fcn will fail if the seqs aren't the same length
      if(hamming_fraction > args_->hamming_fraction_cutoff()) {
	// if(args_->debug()) printf("         hamming %.2f > %.2f for %s %s\n", hamming_fraction, args_->hamming_fraction_cutoff(), key_a.c_str(), key_b.c_str());
	++n_skipped_hamming;
	continue;
      }

      DPHandler dph(args_, gl_, hmms_, qmerged.only_genes_);  // NOTE it's an ok approximation to compare log probs for sequence sets that were run with different kbounds, but (I'm pretty sure) we do need to run them with the same set of genes. EDIT hm, well, maybe not. Anywa, it's ok for now

      // NOTE the error from using the single kbounds rather than the OR seems to be around a part in a thousand or less
      GetLogProb(dph, key_a, *a_seqs, kbinfo_[key_a], mute_freqs_[key_a]);
      GetLogProb(dph, key_b, *b_seqs, kbinfo_[key_b], mute_freqs_[key_b]);
      GetLogProb(dph, qmerged.name_, qmerged.seqs_, qmerged.kbounds_, qmerged.mean_mute_freq_);

      double bayes_factor(log_probs_[qmerged.name_] - log_probs_[key_a] - log_probs_[key_b]);  // REMINDER a, b not necessarily same order as names[0], names[1]
      if(args_->debug()) {
	printf("       %8.3f = ", bayes_factor);
	printf("%2s %8.2f", "", log_probs_[qmerged.name_]);
	printf(" - %8.2f - %8.2f", log_probs_[key_a], log_probs_[key_b]);
	printf("\n");
      }

      potential_merges.push_back(pair<double, Query>(bayes_factor, qmerged));

      if(bayes_factor > max_log_prob) {
	max_log_prob = bayes_factor;
	// max_pair = pair<string, string>(key_a, key_b);  // REMINDER not always same order as names[0], names[1]
	imax = potential_merges.size() - 1;
      }
    }
  }

  // if <info_> only has one cluster, if hamming is too large between all remaining clusters, or if remaining bayes factors are -INFINITY
  if(max_log_prob == -INFINITY) {
    path->finished_ = true;
    return;
  }

  Query *qmerged(nullptr);
  double chosen_logprob;
  if(args_->smc_particles() == 1) {
    qmerged = &potential_merges[imax].second;
    chosen_logprob = potential_merges[imax].first;
  } else {
    qmerged = ChooseRandomMerge(potential_merges, rgen);
    chosen_logprob = max_log_prob;
  }
  // cout << "chose: " << qmerged->name_ << endl;
  // Query qmerged(GetMergedQuery(max_pair.first, max_pair.second));  // NOTE I'm still kinda worried about ordering shenanigans here. I think the worst case is we calculate something twice that we don't need to, but I'm still nervous.
  
  // merge the two best clusters
  // string max_name_str = JoinNames(max_pair.first, max_pair.second);
  // if(info_.count(max_name_str) == 0) {  // if we're doing smc, this will happen once for each particle that wants to merge these two. NOTE you get a very, very strange seg fault at the Sequences::Union line above, I *think* inside the implicit copy constructor. Yes, I should just define my own copy constructor, but I can't work out how to navigate through the const jungle a.t.m.
  if(info_.count(qmerged->name_) == 0) {  // if we're doing smc, this will happen once for each particle that wants to merge these two. NOTE you get a very, very strange seg fault at the Sequences::Union line above, I *think* inside the implicit copy constructor. Yes, I should just define my own copy constructor, but I can't work out how to navigate through the const jungle a.t.m.
    // Query qmerged(GetMergedQuery(max_pair.first, max_pair.second));  // NOTE I'm still kinda worried about ordering shenanigans here. I think the worst case is we calculate something twice that we don't need to, but I'm still nervous.
    info_[qmerged->name_] = qmerged->seqs_;
    kbinfo_[qmerged->name_] = qmerged->kbounds_;
    mute_freqs_[qmerged->name_] = qmerged->mean_mute_freq_;
    only_genes_[qmerged->name_] = qmerged->only_genes_;

    GetNaiveSeq(qmerged->name_);
  }

  Partition new_partition(path->CurrentPartition());  // note: CurrentPartition() returns a reference
  // new_partition.erase(max_pair.first);
  // new_partition.erase(max_pair.second);
  // new_partition.insert(max_name_str);
  new_partition.erase(qmerged->parents_.first);
  new_partition.erase(qmerged->parents_.second);
  new_partition.insert(qmerged->name_);

  path->AddPartition(new_partition, LogProbOfPartition(new_partition));

  if(args_->debug()) {
    printf("          hamming skipped %d / %d\n", n_skipped_hamming, n_total_pairs);
    // printf("       merged %-8.2f %s and %s\n", max_log_prob, max_pair.first.c_str(), max_pair.second.c_str());
    printf("       merged %-8.2f %s and %s\n", chosen_logprob, qmerged->parents_.first.c_str(), qmerged->parents_.second.c_str());
    PrintPartition(new_partition, "current");
  }
}

}
