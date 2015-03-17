#include "glomerator.h"

namespace ham {

Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<Sequences> &qry_seq_list, Args *args, Track *track):  // reminder: <qry_seq_list> is a list of lists
  hmms_(hmms),
  gl_(gl),
  args_(args),
  max_log_prob_of_partition_(-INFINITY),
  finished_(false)
{
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

  ReadCachedLogProbs(args->incachefile(), track);
  ofs_.open(args_->outfile());
  ofs_ << "partition,score,errors" << endl;
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  // first add the initial partition
  vector<string> initial_partition(GetClusterList(info_));
  list_of_partitions_.push_back(pair<double, vector<string> >(LogProbOfPartition(initial_partition), initial_partition));
  if(args_->debug())
    PrintPartition(initial_partition, "initial");

  // then agglomerate
  if(args_->debug()) cout << "   glomerating" << endl;
  do {
    Merge();
  } while(!finished_);

  if(args_->debug()) {
    cout << "  -----------------" << endl;  
    PrintPartition(best_partition_, "best");
  }

  // TODO oh, damn, if the initial partition is better than any subsequent ones this breaks
  WritePartitions();

  // TODO really pass the cachefile inside args_? maybe do something cleaner
  WriteCachedLogProbs();

  ofs_.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::ReadCachedLogProbs(string fname, Track *track) {
  if(fname == "")
    return;

  ifstream ifs(fname);
  assert(ifs.is_open());
  string line;

  // check the header is right TODO should write a general csv reader
  getline(ifs, line);
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  vector<string> headstrs(SplitString(line, ","));
  // cout << "x" << headstrs[0] << "x" << headstrs[1] << "x" << headstrs[2] << "x" << endl;
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
      naive_seqs_[unique_ids] = naive_seq;  //Sequence(track, unique_ids, naive_seq_str);
  }
}

// ----------------------------------------------------------------------------------------
// return a vector consisting of the keys in <partinfo>
vector<string> Glomerator::GetClusterList(map<string, Sequences> &partinfo) {
  vector<string> clusters;
  for(auto &kv : partinfo)
    clusters.push_back(kv.first);
  return clusters;
}

// ----------------------------------------------------------------------------------------
void Glomerator::GetSoloLogProb(string key) {
  // TODO the only reason to have this separate from GetLogProb is the only_genes stuff... which needs to be thought about
  JobHolder jh(gl_, hmms_, "forward", only_genes_[key]);
  jh.SetDebug(args_->debug());
  jh.SetChunkCache(args_->chunk_cache());
  jh.SetNBestEvents(args_->n_best_events());
  string errors;
  GetLogProb(jh, key, info_[key], kbinfo_[key], errors);  // TODO make sure these errors will get stored somewhere
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(vector<string> &clusters) {
  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  for(auto &key : clusters) {
    if(log_probs_.count(key) == 0)
      GetSoloLogProb(key);  // should only happen for the initial partition
      // throw runtime_error("ERROR couldn't find key " + key + " in cached log probs\n");
    total_log_prob = AddWithMinusInfinities(total_log_prob, log_probs_[key]);
  }
  return total_log_prob;
}

// ----------------------------------------------------------------------------------------
void Glomerator::PrintPartition(vector<string> &clusters, string extrastr) {
  const char *extra_cstr(extrastr.c_str());  // dammit I shouldn't need this line
  printf("    %-8.2f %s partition\n", LogProbOfPartition(clusters), extra_cstr);
  for(auto &key : clusters)
    cout << "          " << key << endl;
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteCachedLogProbs() {
  ofstream log_prob_ofs(args_->outcachefile());
  assert(log_prob_ofs.is_open());
  log_prob_ofs << "unique_ids,score,naive-seq" << endl;
  for(auto &kv : log_probs_) {
    // if(naive_seqs_.count(kv.first) == 0) {
    //   cout << "LATE " << kv.first << endl;
    //   GetNaiveSeq(kv.first);
    //   // throw runtime_error("ERROR don't have naive seq for " + kv.first + "\n");
    // }
    log_prob_ofs << kv.first << "," << kv.second << "," << naive_seqs_[kv.first] << endl;  // NOTE if we didn't calculate the naive sequence, it'll print as null string, which is fine
  }
  log_prob_ofs.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WritePartitions() {
  for(auto &pr : list_of_partitions_) {
    for(size_t ic=0; ic<pr.second.size(); ++ic) {
      if(ic > 0)
	ofs_ << ";";
      ofs_ << pr.second[ic];
    }
    ofs_ << ",";
    ofs_ << pr.first << ",";
    ofs_ << "n/a\n";
  }
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
// ----------------------------------------------------------------------------------------
int Glomerator::MinimalHammingDistance(Sequences &seqs_a, Sequences &seqs_b) {
  // Minimal hamming distance between any sequence in <seqs_a> and any sequence in <seqs_b>
  // NOTE for now, we require sequences are the same length (otherwise we have to deal with alignming them which is what we would call a CAN OF WORMS.
  assert(seqs_a.n_seqs() > 0 && seqs_b.n_seqs() > 0);
  int min_distance(seqs_a[0].size());
  for(size_t is=0; is<seqs_a.n_seqs(); ++is) {
    for(size_t js=0; js<seqs_b.n_seqs(); ++js) {
      int distance = HammingDistance(seqs_a[is], seqs_b[js]);
      if(distance < min_distance)
	min_distance = distance;
    }
  }
  return min_distance;
}

// ----------------------------------------------------------------------------------------
// add log prob for <key>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetNaiveSeq(string key) {
  if(naive_seqs_.count(key))  // already did it
    return;

  string errors;  // TODO make sure these errors get collected as well
  JobHolder jh(gl_, hmms_, "viterbi", only_genes_[key]);
  jh.SetDebug(args_->debug());
  jh.SetChunkCache(args_->chunk_cache());
  jh.SetNBestEvents(args_->n_best_events());
  Result result(kbinfo_[key]);
  bool stop(false);
  do {
    result = jh.Run(info_[key], kbinfo_[key]);
    kbinfo_[key] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.events_.size() < 1)
    throw runtime_error("ERROR no events for" + key + "\n");
  // naive_seqs_[key] = Sequence(info_[key][0].track(), key, result.events_[0].naive_seq_);  // TODO it might be a bit wasteful to store these digitized? Or maybe it's better...
  naive_seqs_[key] = result.events_[0].naive_seq_;  // TODO it might be a bit wasteful to store these digitized? Or maybe it's better...
  if(result.boundary_error())
    errors = "boundary";
}

// ----------------------------------------------------------------------------------------
// add log prob for <name>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetLogProb(JobHolder &jh, string name, Sequences &seqs, KBounds &kbounds, string &errors) {
  // TODO note that when this imporves the kbounds, that info doesn't get propagated to <kbinfo_>
  if(log_probs_.count(name))  // already did it
    return;
    
  Result result(kbounds);
  bool stop(false);
  do {
    result = jh.Run(seqs, kbounds);
    kbounds = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  log_probs_[name] = result.total_score();
  if(result.boundary_error())
    errors = "boundary";
}

// ----------------------------------------------------------------------------------------
// perform one merge step, i.e. find the two "nearest" clusters and merge 'em
void Glomerator::Merge() {
  if(args_->naive_preclustering())
    cout << "merging with naive preclustering" << endl;
  else
    cout << "merging with minimal mature preclustering" << endl;
  double max_log_prob(-INFINITY);
  pair<string, string> max_pair; // pair of clusters with largest log prob (i.e. the ratio of their prob together to their prob apart is largest)
  vector<string> max_only_genes;
  KBounds max_kbounds;

  int n_total_pairs(0), n_skipped_hamming(0);

  set<string> already_done;  // keeps track of which  a-b pairs we've already done
  for(auto &kv_a : info_) {  // note that c++ maps are ordered
    for(auto &kv_b : info_) {
      if(kv_a.first == kv_b.first) continue;
      vector<string> names{kv_a.first, kv_b.first};
      sort(names.begin(), names.end());
      string bothnamestr(names[0] + ":" + names[1]);
      if(already_done.count(bothnamestr))  // already did this pair
	continue;
      else
	already_done.insert(bothnamestr);

      ++n_total_pairs;
      Sequences a_seqs(kv_a.second), b_seqs(kv_b.second);  // TODO cache hamming fraction as well

      if(args_->naive_preclustering()) {
	double hamming_fraction = float(NaiveHammingDistance(kv_a.first, kv_b.first)) / a_seqs[0].size();  // hamming distance fcn will fail if the seqs aren't the same length
	if(hamming_fraction > args_->hamming_fraction_cutoff()) {
	  ++n_skipped_hamming;
	  continue;
	}
      } else {
	double hamming_fraction = float(MinimalHammingDistance(a_seqs, b_seqs)) / a_seqs[0].size();  // minimal_hamming_distance() will fail if the seqs aren't all the same length
	if(hamming_fraction > args_->hamming_fraction_cutoff()) {  // if all sequences in a are too far away from all sequences in b
	  ++n_skipped_hamming;
	  continue;
	}
      }

      // TODO skip all this stuff if we already have all three of 'em cached
      Sequences ab_seqs(a_seqs.Union(b_seqs));
      vector<string> ab_only_genes = only_genes_[kv_a.first];
      for(auto &g : only_genes_[kv_b.first])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
	ab_only_genes.push_back(g);
      KBounds ab_kbounds = kbinfo_[kv_a.first].LogicalOr(kbinfo_[kv_b.first]);

      JobHolder jh(gl_, hmms_, "forward", ab_only_genes);  // NOTE it's an ok approximation to compare log probs for sequence sets that were run with different kbounds, but (I'm pretty sure) we do need to run them with the same set of genes. EDIT hm, well, maybe not. Anywa, it's ok for now
      // TODO make sure that using <ab_only_genes> doesn't introduce some bias
      jh.SetDebug(args_->debug());
      jh.SetChunkCache(args_->chunk_cache());
      jh.SetNBestEvents(args_->n_best_events());

      string errors;
      // NOTE the error from using the single kbounds rather than the OR seems to be around a part in a thousand or less
      GetLogProb(jh, kv_a.first, a_seqs, kbinfo_[kv_a.first], errors);
      GetLogProb(jh, kv_b.first, b_seqs, kbinfo_[kv_b.first], errors);
      GetLogProb(jh, bothnamestr, ab_seqs, ab_kbounds, errors);

      double bayes_factor(log_probs_[bothnamestr] - log_probs_[kv_a.first] - log_probs_[kv_b.first]);  // REMINDER a, b not necessarily same order as names[0], names[1]
      if(args_->debug()) {
	printf("       %8.3f = ", bayes_factor);
	printf("%2s %8.2f", "", log_probs_[bothnamestr]);
	printf(" - %8.2f - %8.2f", log_probs_[kv_a.first], log_probs_[kv_b.first]);
	printf("\n");
      }

      if(bayes_factor > max_log_prob) {
	max_log_prob = bayes_factor;
	max_pair = pair<string, string>(kv_a.first, kv_b.first);  // REMINDER not always same order as names[0], names[1]
	max_only_genes = ab_only_genes;
	max_kbounds = ab_kbounds;
      }
    }
  }

  // if <info_> only has one cluster, if hamming is too large between all remaining clusters/remaining bayes factors are -INFINITY
  if(max_log_prob == -INFINITY) {
    finished_ = true;
    return;
  }

  // then merge the two best clusters
  vector<string> max_names{max_pair.first, max_pair.second};
  sort(max_names.begin(), max_names.end());
  Sequences max_seqs(info_[max_names[0]].Union(info_[max_names[1]]));  // NOTE this will give the ordering {<first seqs>, <second seqs>}, which should be the same as in <max_name_str>. But I don't think it'd hurt anything if the sequences and names were in a different order
  string max_name_str(max_names[0] + ":" + max_names[1]);  // NOTE the names[i] are not sorted *within* themselves, but <names> itself is sorted. This is ok, because we will never again encounter these sequences separately
  info_[max_name_str] = max_seqs;
  kbinfo_[max_name_str] = max_kbounds;
  only_genes_[max_name_str] = max_only_genes;

  GetNaiveSeq(max_name_str);

  info_.erase(max_pair.first);
  info_.erase(max_pair.second);
  kbinfo_.erase(max_pair.first);
  kbinfo_.erase(max_pair.second);
  only_genes_.erase(max_pair.first);
  only_genes_.erase(max_pair.second);

  printf("          hamming skipped %d / %d = %.2f\n", n_skipped_hamming, n_total_pairs, float(n_skipped_hamming) / n_total_pairs);
  if(args_->debug()) {
    printf("          hamming skipped %d / %d\n", n_skipped_hamming, n_total_pairs);
    printf("       merged %-8.2f %s and %s\n", max_log_prob, max_pair.first.c_str(), max_pair.second.c_str());
  }

  vector<string> partition(GetClusterList(info_));
  if(args_->debug())
    PrintPartition(partition, "current");
  double total_log_prob = LogProbOfPartition(partition);

  list_of_partitions_.push_back(pair<double, vector<string> >(total_log_prob, partition));

  if(total_log_prob > max_log_prob_of_partition_) {
    best_partition_ = partition;
    max_log_prob_of_partition_ = total_log_prob;
  }

  if(max_log_prob_of_partition_ - total_log_prob > 1000.0) {  // stop if we've moved too far past the maximum
    cout << "    stopping after drop " << max_log_prob_of_partition_ << " --> " << total_log_prob << endl;
    finished_ = true;  // NOTE this will not play well with multiple maxima, but I'm pretty sure we shouldn't be getting those
  }
}

}
