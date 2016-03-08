#include "glomerator.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track) :
  track_(track),
  args_(args),
  vtb_dph_("viterbi", args_, gl, hmms),
  fwd_dph_("forward", args_, gl, hmms),
  i_initial_partition_(0),
  n_max_factor_(1.5),
  n_fwd_calculated_(0),
  n_vtb_calculated_(0),
  n_hfrac_calculated_(0),
  n_hamming_merged_(0),
  progress_file_(fopen((args_->outfile() + ".progress").c_str(), "w"))
{
  time(&last_status_write_time_);
  ReadCachedLogProbs();

  name_translations_[args_->biggest_naive_seq_cluster_to_calculate()] = map<string, string>();
  name_translations_[args_->biggest_logprob_cluster_to_calculate()] = map<string, string>();

  Partition tmp_partition;
  int last_ipath(0);
  assert(args_->integers_["path_index"].size() == qry_seq_list.size());
  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key = SeqNameStr(qry_seq_list[iqry], ":");

    int ipath(args_->integers_["path_index"][iqry]);
    if(last_ipath != ipath) {  // if we need to push back a new initial partition (i.e. this is a new path/particle) NOTE I'm assuming the the first one will have <path_index> zero
      initial_partitions_.push_back(tmp_partition);
      initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));
      initial_logweights_.push_back(args_->floats_["logweight"][iqry]);  // it's the same for each <iqry> with this <path_index>
      tmp_partition.clear();
      last_ipath = ipath;
    }
      
    tmp_partition.insert(key);

    if(seq_info_.count(key) > 0)  // already added this cluster to the maps (but it'll appear more than once if it's in multiple paths/particles), so we only need to add the cluster info to <initial_partitions_>
      continue;

    seq_info_[key] = qry_seq_list[iqry];
    only_genes_[key] = args_->str_lists_["only_genes"][iqry];

    KSet kmin(args_->integers_["k_v_min"][iqry], args_->integers_["k_d_min"][iqry]);
    KSet kmax(args_->integers_["k_v_max"][iqry], args_->integers_["k_d_max"][iqry]);
    KBounds kb(kmin, kmax);
    kbinfo_[key] = kb;

    vector<KBounds> kbvector(seq_info_[key].size(), kbinfo_[key]);

    mute_freqs_[key] = avgVector(args_->float_lists_["mute_freqs"][iqry]);
  }
  // add the last initial partition (i.e. for the last path/particle)
  assert(tmp_partition.size() > 0);
  initial_partitions_.push_back(tmp_partition);
  initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));
  initial_logweights_.push_back(args_->floats_["logweight"].back());

  // the *first* time, we only get one path/partition from partis, so we push back a bunch of copies
  if((int)initial_partitions_.size() == 1 && args_->smc_particles() > 1)  {
    for(int ip=1; ip<args_->smc_particles(); ++ip) {
      initial_partitions_.push_back(tmp_partition);
      initial_logprobs_.push_back(LogProbOfPartition(tmp_partition));
      initial_logweights_.push_back(args_->floats_["logweight"].back());
    }
  }

  if((int)initial_partitions_.size() != args_->smc_particles())
    throw runtime_error("wrong number of initial partitions " + to_string(initial_partitions_.size()) + " (should be " + to_string(args_->smc_particles()) + ")");

  if(args_->debug())
    for(auto &part : initial_partitions_)
      PrintPartition(part, "initial");
}

// ----------------------------------------------------------------------------------------
Glomerator::~Glomerator() {
  printf("         calculated   vtb %-4d   fwd %-4d   hamming merged %-4d   naive hfracs %d\n", n_vtb_calculated_, n_fwd_calculated_, n_hamming_merged_, n_hfrac_calculated_);
  if(args_->cachefile() != "")
    WriteCachedLogProbs();

  // // ----------------------------------------------------------------------------------------
  // ifstream smapfs("/proc/self/smaps");
  // string tmpline;
  // while(getline(smapfs, tmpline)) {
  //   cout << "MEM " << tmpline << endl;
  // }
  // // ----------------------------------------------------------------------------------------
  fclose(progress_file_);
  remove((args_->outfile() + ".progress").c_str());
}

// ----------------------------------------------------------------------------------------
void Glomerator::CacheNaiveSeqs() {  // they're written to file in the destructor, so we just need to calculate them here
  for(auto &kv : seq_info_)
    GetNaiveSeq(kv.first);
  ofs_.open(args_->outfile());  // a.t.m. I'm signalling that I finished ok by doing this
  ofs_.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  if(args_->debug()) cout << "   glomerating" << endl;

  if(args_->logprob_ratio_threshold() == -INFINITY)
    throw runtime_error("logprob ratio threshold not specified");

  assert((int)initial_partitions_.size() == 1);
  ClusterPath cp(initial_partitions_[0], LogProbOfPartition(initial_partitions_[0]), initial_logweights_[0]);
  do {
    Merge(&cp);
    // cout << ClusterSizeString(&cp) << endl;
    WriteStatus(&cp);
  } while(!cp.finished_);

  vector<ClusterPath> paths{cp};
  WritePartitions(paths);
  if(args_->annotationfile() != "")
    WriteAnnotations(paths);
}

// ----------------------------------------------------------------------------------------
Partition Glomerator::GetAnInitialPartition(int &initial_path_index, double &logweight) {
  initial_path_index = i_initial_partition_;
  assert(i_initial_partition_ < (int)initial_partitions_.size());
  logweight = initial_logweights_[i_initial_partition_];
  return initial_partitions_.at(i_initial_partition_++);
}

// ----------------------------------------------------------------------------------------
void Glomerator::ReadCachedLogProbs() {
  ifstream ifs(args_->cachefile());
  if(!ifs.is_open()) {  // this means we don't have any cached results to start with, but we'll write out what we have at the end of the run to this file
    cout << "        cachefile d.n.e." << endl;
    return;
  }
  string line;
  // check the header is right (no cached info)
  if(!getline(ifs, line)) {
    cout << "        empty cachefile" << endl;
    return;  // return for zero length file
  }
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  vector<string> headstrs(SplitString(line, ","));
  assert(headstrs[0].find("unique_ids") == 0);
  assert(headstrs[1].find("logprob") == 0);
  assert(headstrs[2].find("naive_seq") == 0);
  assert(headstrs[3].find("naive_hfrac") == 0);

  // NOTE there can be two lines with the same key (say if in one run we calculated the naive seq, and in a later run calculated the log prob)
  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> column_list = SplitString(line, ",");
    assert(column_list.size() == 5);
    string query(column_list[0]);

    string logprob_str(column_list[1]);
    if(logprob_str.size() > 0) {
      log_probs_[query] = stof(logprob_str);
      initial_log_probs_.insert(query);
    }

    string naive_seq(column_list[2]);

    string naive_hfrac_str(column_list[3]);
    if(naive_hfrac_str.size() > 0) {
      naive_hfracs_[query] = stof(naive_hfrac_str);
      initial_naive_hfracs_.insert(query);
    }

    if(naive_seq.size() > 0) {
      naive_seqs_[query] = naive_seq;
      initial_naive_seqs_.insert(query);
    }
  }
  cout << "        read " << log_probs_.size() << " cached logprobs and " << naive_seqs_.size() << " naive seqs" << endl;
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteCacheLine(ofstream &ofs, string query) {
  ofs << query << ",";
  if(log_probs_.count(query))
    ofs << log_probs_[query];
  ofs << ",";
  if(naive_seqs_.count(query))
    ofs << naive_seqs_[query];
  ofs << ",";
  if(args_->cache_naive_hfracs() && naive_hfracs_.count(query))
    ofs << naive_hfracs_[query];
  ofs << ",";
  if(errors_.count(query))
    ofs << errors_[query];
  ofs << endl;
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteCachedLogProbs() {
  ofstream log_prob_ofs(args_->cachefile());
  if(!log_prob_ofs.is_open())
    throw runtime_error("ERROR cache file (" + args_->cachefile() + ") d.n.e.\n");
  log_prob_ofs << "unique_ids,logprob,naive_seq,naive_hfrac,errors" << endl;

  log_prob_ofs << setprecision(20);

  set<string> keys_to_cache;
  for(auto &kv : log_probs_) {
    if(args_->only_cache_new_vals() && initial_log_probs_.count(kv.first))  // don't cache it if we had it in the initial cache file (this is just an optimization)
      continue;
    keys_to_cache.insert(kv.first);
  }
  for(auto &kv : naive_seqs_) {
    if(args_->only_cache_new_vals() && initial_naive_seqs_.count(kv.first))  // note that if we had an initial log prob, but not an initial naive seq, we *do* want to write it (if we calculated the naive seq)
      continue;
    keys_to_cache.insert(kv.first);
  }
  if(args_->cache_naive_hfracs()) {
    for(auto &kv : naive_hfracs_) {
      if(args_->only_cache_new_vals() && initial_naive_hfracs_.count(kv.first))
	continue;
      keys_to_cache.insert(kv.first);
    }
  }

  for(auto &key : keys_to_cache)
    WriteCacheLine(log_prob_ofs, key);

  log_prob_ofs.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WritePartitions(vector<ClusterPath> &paths) {
  ofs_.open(args_->outfile());
  ofs_ << setprecision(20);
  if(paths.size() > 0)
    ofs_ << "path_index,initial_path_index,";
  ofs_ << "partition,logprob,logweight" << endl;
  int ipath(0);
  for(auto &cp : paths) {
    for(unsigned ipart=0; ipart<cp.partitions().size(); ++ipart) {
      if(paths.size() > 0)
	ofs_ << ipath << "," << cp.initial_path_index_ << ",";
      int ic(0);
      for(auto &cluster : cp.partitions()[ipart]) {
	if(ic > 0)
	  ofs_ << ";";
	ofs_ << cluster;
	++ic;
      }
      ofs_ << "," << cp.logprobs()[ipart]
	   << "," << cp.logweights()[ipart] << endl;
    }
    ++ipath;
  }
  ofs_.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteAnnotations(vector<ClusterPath> &paths) {
  throw runtime_error("needs updating -- specifically need to make sure that when we replace one naive seq with another, we also replace other things like the event_ info ");
  ofstream annotation_ofs;
  annotation_ofs.open(args_->annotationfile());
  StreamHeader(annotation_ofs, "viterbi");

  assert(paths.size() == 1);  // would need to update this for smc
  int ipath(0);
  ClusterPath cp(paths[ipath]);
  for(unsigned ipart=0; ipart<cp.partitions().size(); ++ipart) {  // we don't work out which is the best partition until later (in the python), so darn, I guess I'll just write annotations for all the partitions
    for(auto &cluster : cp.partitions()[ipart]) {
      if(events_[cluster].genes_["d"] == "") {  // shouldn't happen any more, but it is a check that could fail at some point
	cout << "WTF " << cluster << " x" << events_[cluster].naive_seq_ << "x" << endl;
	assert(0);
      }
      vector<RecoEvent> event_list({events_[cluster]});
      StreamOutput(annotation_ofs, "viterbi", 1, event_list, seq_info_[cluster], 0., "");
    }
  }
  annotation_ofs.close();
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(Partition &partition, bool debug) {
  if(args_->no_fwd())  // if we're doing pure naive hamming glomeration, we don't want to calculate any forward probs
    return -INFINITY;

  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  if(debug)
    cout << "LogProbOfPartition: " << endl;
  for(auto &key : partition) {
    // assert(SameLength(seq_info_[key], true));
    double log_prob = GetLogProb(key);
    if(debug)
      cout << "  " << log_prob << "  " << key << endl;
    total_log_prob = AddWithMinusInfinities(total_log_prob, log_prob);
  }
  if(debug)
    cout << "  total: " << total_log_prob << endl;
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
void Glomerator::WriteStatus(ClusterPath *path) {
  time_t current_time;
  time(&current_time);
  if(difftime(current_time, last_status_write_time_) > 300) {  // write something every five minutes
    char buffer[200];
    strftime(buffer, 200, "%b %d %T", localtime(&current_time));  // %H:%M
    fprintf(progress_file_, "      %s    %4d clusters    fwd %-4d   vtb %-4d\n", buffer, (int)path->CurrentPartition().size(), n_fwd_calculated_, n_vtb_calculated_);

    fprintf(progress_file_, "%s\n", ClusterSizeString(path).c_str());

    fflush(progress_file_);
    last_status_write_time_ = current_time;
  }
}

// ----------------------------------------------------------------------------------------
string Glomerator::ParentalString(pair<string, string> *parents) {
  if(CountMembers(parents->first) > 5 || CountMembers(parents->second) > 5) {
    return to_string(CountMembers(parents->first)) + " and " + to_string(CountMembers(parents->second));
  } else {
    return parents->first + " and " + parents->second;
  }
}

// ----------------------------------------------------------------------------------------
// count the number of members in a cluster's colon-separated name string
int Glomerator::CountMembers(string namestr) {
  int n_colons = (int)count(namestr.begin(), namestr.end(), ':');
  return n_colons + 1;
}

// ----------------------------------------------------------------------------------------
string Glomerator::ClusterSizeString(ClusterPath *path) {
  vector<int> cluster_sizes;
  for(auto &cluster : path->CurrentPartition()) {
    cluster_sizes.push_back(CountMembers(cluster));
  }
  sort(cluster_sizes.begin(), cluster_sizes.end());
  reverse(cluster_sizes.begin(), cluster_sizes.end());
  string return_str("          clusters: ");
  for(size_t is=0; is<cluster_sizes.size(); ++is)
    return_str += " "  + to_string(cluster_sizes[is]);
  return return_str;
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinNames(string name1, string name2, string delimiter) {
  vector<string> names{name1, name2};
  sort(names.begin(), names.end());  // NOTE this doesn't sort *within* name1 or name2 when they're already comprised of several uids. In principle this will lead to unnecessary cache misses (if we later arrive at the same combination of sequences from a different starting point). In practice, this is very unlikely (unless we're dong smc) since we've already merged the constituents of name1 and name2 and we can't unmerge them.
  return names[0] + delimiter + names[1];
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinNameStrings(vector<Sequence> &strlist, string delimiter) {
  string return_str;
  for(size_t is=0; is<strlist.size(); ++is) {
    if(is > 0)
      return_str += delimiter;
    return_str += strlist[is].name();
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinSeqStrings(vector<Sequence> &strlist, string delimiter) {
  string return_str;
  for(size_t is=0; is<strlist.size(); ++is) {
    if(is > 0)
      return_str += delimiter;
    return_str += strlist[is].undigitized();
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
double Glomerator::CalculateHfrac(string &seq_a, string &seq_b) {
  ++n_hfrac_calculated_;
  if(seq_a.size() != seq_b.size())
    throw runtime_error("ERROR sequences different length in Glomerator::NaiveHfrac (" + seq_a + "," + seq_b + ")\n");
  int distance(0), len_excluding_ambigs(0);
  for(size_t ic=0; ic<seq_a.size(); ++ic) {
    uint8_t ch_a = track_->symbol_index(seq_a.substr(ic, 1));  // kind of hackey remnant left from when naive seqs were Sequence objects
    uint8_t ch_b = track_->symbol_index(seq_b.substr(ic, 1));
    if(ch_a == track_->ambiguous_index() || ch_b == track_->ambiguous_index())  // skip this position if either sequence has an ambiguous character (if not set, ambig-base should be the empty string)
      continue;
    ++len_excluding_ambigs;
    if(ch_a != ch_b)
      ++distance;
  }

  return distance / double(len_excluding_ambigs);
}

// ----------------------------------------------------------------------------------------
double Glomerator::NaiveHfrac(string key_a, string key_b) {
  string joint_key = JoinNames(key_a, key_b);  // NOTE since the cache is indexed by the joint key, this assumes we can arrive at this cluster via only one path. Which should be ok.
  if(naive_hfracs_.count(joint_key))  // if we've already calculated this distance
    return naive_hfracs_[joint_key];

  string &seq_a = GetNaiveSeq(key_a);
  string &seq_b = GetNaiveSeq(key_b);
  naive_hfracs_[joint_key] = CalculateHfrac(seq_a, seq_b);

  return naive_hfracs_[joint_key];
}

// ----------------------------------------------------------------------------------------
pair<string, vector<Sequence> > Glomerator::ChooseSubsetOfNames(string queries, int n_max) {
  vector<string> subqueries;
  vector<Sequence> subseqs;
  vector<string> namevector(SplitString(queries, ":"));

  srand(hash<string>{}(queries));  // make sure we get the same subset each time we pass in the same queries (well, if there's different thresholds for naive_seqs annd logprobs they'll each get their own [very correlated] subset)

  set<int> already_chosen;
  for(size_t iname=0; iname<unsigned(n_max); ++iname) {
    int ichosen(-1);
    while(ichosen < 0 || already_chosen.count(ichosen))
      ichosen = rand() % namevector.size();
    already_chosen.insert(ichosen);
    subqueries.push_back(namevector[ichosen]);
    subseqs.push_back(seq_info_[queries][ichosen]);
  }

  pair<string, vector<Sequence> > substuff(JoinStrings(subqueries), subseqs);
  return substuff;
}

// ----------------------------------------------------------------------------------------
string Glomerator::GetNameToCalculate(string actual_queries, int n_max) {
  // // NOTE we don't really need this, since we're setting the random seed. But it just seems so messy to go through the whole subset calculation every time, even though I profiled it and it's not a significant contributor
  if(name_translations_[n_max].count(actual_queries))  // make sure we always use the same translation if we already did one
    return name_translations_[n_max][actual_queries];

  // if cluster is more than half again larger than n_max, replace it with a cluster of this size
  if(CountMembers(actual_queries) > n_max_factor_ * n_max) {
    pair<string, vector<Sequence> > substuff = ChooseSubsetOfNames(actual_queries, n_max);
    string subqueries(substuff.first);
    seq_info_[subqueries] = substuff.second;
    if(args_->debug() > 0) {
      if(CountMembers(actual_queries) + CountMembers(subqueries) < 20)
	cout <<  "     replacing " << actual_queries << " --> " << subqueries << endl;
      else
	cout <<  "     replacing " << CountMembers(actual_queries) << " --> " << CountMembers(subqueries) << endl;
    }
    kbinfo_[subqueries] = kbinfo_[actual_queries];  // just use the entire/super cluster for this stuff. It's just overly conservative (as long as you keep the mute freqs the same)
    mute_freqs_[subqueries] = mute_freqs_[actual_queries];
    only_genes_[subqueries] = only_genes_[actual_queries];

    name_translations_[n_max][actual_queries] = subqueries;
    return subqueries;
  }

  // falling through to here means we want to just use <actual_queries>
  return actual_queries;
}

// ----------------------------------------------------------------------------------------
string &Glomerator::GetNaiveSeq(string queries, pair<string, string> *parents) {
  if(naive_seqs_.count(queries))
    return naive_seqs_[queries];

  // if we have naive seqs for both the parental clusters and they're the same, no reason to calculate this naive seq
  if(parents != nullptr && GetNaiveSeq(parents->first) == GetNaiveSeq(parents->second)) {  /// add an entry to <naive_seqs_>, filed under <queries>, which has the same naive sequence as <parentname>
    naive_seqs_[queries] = GetNaiveSeq(parents->first);  // copy the whole sequence object
    events_[queries] = events_[parents->first];  // also copy the event  // TODO this doesn't set all the event info correctly
    events_[queries].seq_name_ = queries;
    return naive_seqs_[queries];
  }

  string queries_to_calc = GetNameToCalculate(queries, args_->biggest_naive_seq_cluster_to_calculate());

  if(naive_seqs_.count(queries_to_calc) == 0)
    naive_seqs_[queries_to_calc] = CalculateNaiveSeq(queries_to_calc);

  if(queries_to_calc != queries) {
    // string full_nseq = CalculateNaiveSeq(queries);
    // string sub_nseq = CalculateNaiveSeq(queries_to_calc);
    // double hfrac = CalculateHfrac(full_nseq, sub_nseq);
    // printf("       use %3d instead of %3d (hamming %5.3f)\n", CountMembers(queries_to_calc), CountMembers(queries), hfrac);
    naive_seqs_[queries] = naive_seqs_[queries_to_calc];
  }

  return naive_seqs_[queries];
}

// // ----------------------------------------------------------------------------------------
// double Glomerator::NormFactor(string name) {
//     double nseqs(CountMembers(name));
//     return nseqs / (1. - 0.24 / pow(nseqs, 0.9));
// }

// ----------------------------------------------------------------------------------------
double Glomerator::GetLogProb(string queries) {
  if(log_probs_.count(queries))  // already did it
    return log_probs_[queries];

  double tmplp = CalculateLogProb(queries);  // NOTE this should be the *only* place (besides cache reading) that log_probs_ gets modified
  log_probs_[queries] = tmplp;  // tmp variable is just so we can assert that queries isn't already in log_probs_

  return log_probs_[queries];
}

// ----------------------------------------------------------------------------------------
double Glomerator::GetLogProbRatio(string key_a, string key_b) {
  // NOTE the error from using the single kbounds rather than the OR seems to be around a part in a thousand or less (it's only really important that the merged query has the OR)
  // NOTE also that the _a and _b results will be cached, but with their *individual* only_gene sets (rather than the OR)... but this seems to be ok.
  // NOTE if kbounds gets expanded in one of these three calls, we don't redo the others. Which is really ok, but could be checked again?
  string key_a_to_calc = GetNameToCalculate(key_a, args_->biggest_logprob_cluster_to_calculate());
  string key_b_to_calc = GetNameToCalculate(key_b, args_->biggest_logprob_cluster_to_calculate());

  // arg, this makes a whole new Query even if we're not replacing
  // TODO wait, don't I want to translate *this*, and *then* see about translating the denominator?
  Query qmerged_to_calc = GetMergedQuery(key_a_to_calc, key_b_to_calc);  // NOTE also enters the merged query's info into seq_info_, kbinfo_, mute_freqs_, and only_genes_

  double log_prob_a = GetLogProb(key_a_to_calc);
  double log_prob_b = GetLogProb(key_b_to_calc);
  double log_prob_ab = GetLogProb(qmerged_to_calc.name_);

  double lratio(log_prob_ab - log_prob_a - log_prob_b);
  if(args_->debug()) {
    printf("       %8.3f = ", lratio);
    printf("%2s %8.2f", "", log_prob_ab);
    printf(" - %8.2f - %8.2f", log_prob_a, log_prob_b);
    printf("\n");
    // if(key_a != key_a_to_calc || key_b != key_b_to_calc) {
    //   Query tmpq(GetMergedQuery(key_a, key_b));  // need this to insert the full merged query's info into seq_info_ and whatnot
    //   double full_lratio = CalculateLogProb(JoinNames(key_a, key_b)) - CalculateLogProb(key_a) - CalculateLogProb(key_b);
    //   double sub_lratio = CalculateLogProb(JoinNames(key_a_to_calc, key_b_to_calc)) - CalculateLogProb(key_a_to_calc) - CalculateLogProb(key_b_to_calc);
    //   double frac_diff = (sub_lratio - full_lratio) / full_lratio;
    //   printf("       use %8.2f instead of %8.2f (%5.3f)\n", sub_lratio, full_lratio, frac_diff);
    // }
  }

  return lratio;
}

// ----------------------------------------------------------------------------------------
string Glomerator::CalculateNaiveSeq(string queries) {
  // NOTE do *not* call this from anywhere except GetNaiveSeq()
  assert(naive_seqs_.count(queries) == 0);

  if(seq_info_.count(queries) == 0 || kbinfo_.count(queries) == 0 || only_genes_.count(queries) == 0 || mute_freqs_.count(queries) == 0)
    throw runtime_error("no info for " + queries);

  ++n_vtb_calculated_;

  Result result(kbinfo_[queries]);
  bool stop(false);
  do {
    // assert(SameLength(seq_info_[queries], true));
    result = vtb_dph_.Run(seq_info_[queries], kbinfo_[queries], only_genes_[queries], mute_freqs_[queries]);  // NOTE the sequences in <seq_info_[queries]> should already be the same length, since they've already been merged
    kbinfo_[queries] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.events_.size() < 1)
    throw runtime_error("no events for queries " + queries + "\n");
  events_[queries] = result.events_[0];  // NOTE keeping separate from naive_seqs_ (at least for now) because I only need the full event for the final partition (UPDATE or do I only use it to write annotations)
  // TODO get events_ on same footing as naive_seqs_
  if(result.boundary_error())
    errors_[queries] = errors_[queries] + ":boundary";

  return result.events_[0].naive_seq_;
}

// ----------------------------------------------------------------------------------------
double Glomerator::CalculateLogProb(string queries) {  // NOTE can modify kbinfo_
  // NOTE do *not* call this from anywhere except GetLogProb()
  assert(log_probs_.count(queries) == 0);

  if(seq_info_.count(queries) == 0 || kbinfo_.count(queries) == 0 || only_genes_.count(queries) == 0 || mute_freqs_.count(queries) == 0)
    throw runtime_error("no info for " + queries);
  
  ++n_fwd_calculated_;

  Result result(kbinfo_[queries]);
  bool stop(false);
  do {
    result = fwd_dph_.Run(seq_info_[queries], kbinfo_[queries], only_genes_[queries], mute_freqs_[queries]);  // NOTE <only_genes> isn't necessarily <only_genes_[queries]>, since for the denominator calculation we take the OR
    kbinfo_[queries] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.boundary_error() && !result.could_not_expand())  // could_not_expand means the max is at the edge of the sequence -- e.g. k_d min is 1
    errors_[queries] = errors_[queries] + ":boundary";

  return result.total_score();
}

// ----------------------------------------------------------------------------------------
vector<Sequence> Glomerator::MergeSeqVectors(string name_a, string name_b) {
  // first merge the two vectors
  vector<Sequence> merged_seqs;  // NOTE doesn't work if you preallocate and use std::vector::insert(). No, I have no *@#*($!#ing idea why
  for(size_t is=0; is<seq_info_[name_a].size(); ++is)
    merged_seqs.push_back(seq_info_[name_a][is]);
  for(size_t is=0; is<seq_info_[name_b].size(); ++is)
    merged_seqs.push_back(seq_info_[name_b][is]);

  // then make sure we don't have the same sequence twice
  set<string> all_names;
  for(size_t is=0; is<merged_seqs.size(); ++is) {
    string name(merged_seqs[is].name());
    if(all_names.count(name)) {
      if(args_->seed_unique_id() != "" && args_->seed_unique_id() == name) {
	// cout << "    found seed uid twice" << endl;
      } else {
	throw runtime_error("tried to add sequence with name " + name + " twice in Glomerator::MergeSeqVectors()");
      }
    } else {
      all_names.insert(name);
    }
  }

  return merged_seqs;
}

// ----------------------------------------------------------------------------------------
bool Glomerator::SameLength(vector<Sequence> &seqs, bool debug) {
  // are all the seqs in <seqs> the same length?
  bool same(true);
  int len(-1);
  for(auto &seq : seqs) {
    if(len < 0)
      len = (int)seq.size();
    if(len != (int)seq.size()) {
      same = false;
      break;
    }
  }
  if(debug) {
    if(same)
      cout << "same length: " << JoinNameStrings(seqs) << endl;
    else
      cout << "not same length! " << JoinNameStrings(seqs) << "\n" << JoinSeqStrings(seqs, "\n") << endl;
  }

  return same;
}  

// ----------------------------------------------------------------------------------------
Query Glomerator::GetMergedQuery(string name_a, string name_b) {
  Query qmerged;
  qmerged.name_ = JoinNames(name_a, name_b);  // sorts name_a and name_b, but *doesn't* sort within them
  qmerged.seqs_ = MergeSeqVectors(name_a, name_b);

  assert(SameLength(seq_info_[name_a]));  // all the seqs for name_a should already be the same length
  assert(SameLength(seq_info_[name_b]));  // ...same for name_b

  qmerged.kbounds_ = kbinfo_[name_a].LogicalOr(kbinfo_[name_b]);

  qmerged.only_genes_ = only_genes_[name_a];
  for(auto &g : only_genes_[name_b])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
    qmerged.only_genes_.push_back(g);
  qmerged.mean_mute_freq_ = (seq_info_[name_a].size()*mute_freqs_[name_a] + seq_info_[name_b].size()*mute_freqs_[name_b]) / double(qmerged.seqs_.size());  // simple weighted average (doesn't account for different sequence lengths)
  qmerged.parents_ = pair<string, string>(name_a, name_b);

  // NOTE now that I'm adding the merged query to the cache info here, I can maybe get rid of the qmerged entirely
  seq_info_[qmerged.name_] = qmerged.seqs_;
  kbinfo_[qmerged.name_] = qmerged.kbounds_;
  mute_freqs_[qmerged.name_] = qmerged.mean_mute_freq_;
  only_genes_[qmerged.name_] = qmerged.only_genes_;

  // // only do this if seed unique id is set
  // int biggest_cluster = args_->biggest_naive_seq_cluster_to_calculate();
  // string queries = name_a;
  // string queries_other = name_b;

  // int nseq = CountMembers(queries);
  // int nseq_other = CountMembers(queries_other);
  // int nmax = n_max_factor_ * biggest_cluster;
  // if(nseq > nmax && float(nseq) / nseq_other > 2. ) {  // if <nseq> is large, and if <nseq> more than twice the size of <nseq_other>, use the existing name translation (for which we should already have a logprob and a naive seq)
  //   name_translations_[nmax][JoinNames(queries, queries_other)] = name_translations_[nmax][queries];
  // }


  return qmerged;
}

// ----------------------------------------------------------------------------------------
pair<double, Query> *Glomerator::ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen) {
  // first leave log space and normalize. NOTE instead of normalizing, we could just use rng::Uniform(0., total)
  vector<double> ratios;  // NOTE *not* a probability: it's the ratio of the probability together to the probability apart
  double total(0.0);
  for(auto &pr : potential_merges) {
    double likelihood_ratio = exp(pr.first);
    total += likelihood_ratio;
    ratios.push_back(likelihood_ratio);
  }
  for(auto &ratio : ratios)
    ratio /= total;

  // then choose one at random according to the probs
  double drawpoint = rgen->Uniform(0., 1.);
  double sum(0.0);
  for(size_t im=0; im<ratios.size(); ++im) {
    sum += ratios[im];
    if(sum > drawpoint) {
      return &potential_merges[im];
    }
  }

  throw runtime_error("fell through in Glomerator::ChooseRandomMerge");
}

// ----------------------------------------------------------------------------------------
bool Glomerator::LikelihoodRatioTooSmall(double lratio, int candidate_cluster_size) {
  int ccs(candidate_cluster_size);  // shorthand
  bool lratio_too_small(false);
  double max_lratio(args_->logprob_ratio_threshold());
  if(ccs == 2 && lratio < max_lratio) {
    lratio_too_small = true;
  } else if(ccs == 3 && lratio < max_lratio - 2.) {  // this subtraction "scheme" is largely heuristic a.t.m.
    lratio_too_small = true;
  } else if(ccs == 4 && lratio < max_lratio - 3.) {
    lratio_too_small = true;
  } else if(ccs == 5 && lratio < max_lratio - 4.) {
    lratio_too_small = true;
  } else if(lratio < max_lratio - 5.) {  // just guessing on the 13... but I don't think the best threshold gets anywhere close to zero (like I had it before...)
    lratio_too_small = true;
  }

  return lratio_too_small;
}

// ----------------------------------------------------------------------------------------
Query Glomerator::ChooseMerge(ClusterPath *path, smc::rng *rgen, double *chosen_lratio) {
  double max_log_prob(-INFINITY), min_hamming_fraction(INFINITY);
  Query min_hamming_merge;
  int imax(-1);
  vector<pair<double, Query> > potential_merges;
  int n_total_pairs(0), n_skipped_hamming(0), n_small_lratios(0), n_inf_factors(0);
  for(Partition::iterator it_a = path->CurrentPartition().begin(); it_a != path->CurrentPartition().end(); ++it_a) {
    for(Partition::iterator it_b = it_a; ++it_b != path->CurrentPartition().end();) {
      string key_a(*it_a), key_b(*it_b);
      if(args_->seed_unique_id() != "") {
	if(key_a.find(args_->seed_unique_id()) == string::npos && key_b.find(args_->seed_unique_id()) == string::npos)  // if a seed unique id was set on the command line, and if it's not in either key (which are colon-joined strings of unique_ids)
	  continue;
      }

      ++n_total_pairs;

      double hfrac = NaiveHfrac(key_a, key_b);
      if(hfrac > args_->hamming_fraction_bound_hi()) {  // if naive hamming fraction too big, don't even consider merging the pair
	++n_skipped_hamming;
	continue;
      }

      if(args_->hamming_fraction_bound_lo() > 0.0 && hfrac < args_->hamming_fraction_bound_lo()) {  // if naive hamming is small enough, merge the pair without running hmm
	if(hfrac < min_hamming_fraction) {
	  min_hamming_fraction = hfrac;
	  min_hamming_merge = GetMergedQuery(key_a, key_b);  // NOTE also enters the merged query's info into seq_info_, kbinfo_, mute_freqs_, and only_genes_
	}
	continue;
      }
      if(min_hamming_fraction < INFINITY) {  // if we have any potential hamming merges, that means we'll do those before we do any hmm tomfoolery
	continue;
      }

      double lratio = GetLogProbRatio(key_a, key_b);

      // don't merge if lratio is small (less than zero, more or less)
      if(LikelihoodRatioTooSmall(lratio, CountMembers(key_a) + CountMembers(key_b))) {
	++n_small_lratios;
	continue;
      }

      potential_merges.push_back(pair<double, Query>(lratio, GetMergedQuery(key_a, key_b))); // NOTE also enters the merged query's info into seq_info_, kbinfo_, mute_freqs_, and only_genes_

      if(lratio == -INFINITY)
	++n_inf_factors;

      if(lratio > max_log_prob) {
	max_log_prob = lratio;
	imax = potential_merges.size() - 1;
      }
    }
  }

  if(min_hamming_fraction < INFINITY) {  // if lower hamming fraction was set, then merge the best (nearest) pair of clusters
    ++n_hamming_merged_;
    *chosen_lratio = -INFINITY;  // er... or something
    if(args_->debug())
      printf("           naive hamming merge %.3f\n", min_hamming_fraction);
    return min_hamming_merge;
  }

  if(args_->debug())
    printf("          hamming skipped %d / %d\n", n_skipped_hamming, n_total_pairs);

  // if <path->CurrentPartition()> only has one cluster, if hamming is too large between all remaining clusters, or if remaining likelihood ratios are -INFINITY
  if(max_log_prob == -INFINITY) {
    if(path->CurrentPartition().size() == 1)
      cout << "        stop with partition of size one" << endl;
    else if(n_skipped_hamming == n_total_pairs)
      cout << "        stop with all " << n_skipped_hamming << " / " << n_total_pairs << " hamming distances greater than " << args_->hamming_fraction_bound_hi() << endl;
    else if(n_inf_factors == n_total_pairs)
      cout << "        stop with all " << n_inf_factors << " / " << n_total_pairs << " likelihood ratios -inf" << endl;
    else
      cout << "        stop for some reason or other with -inf: " << n_inf_factors << "   ham skip: " << n_skipped_hamming << "   small lratios: " << n_small_lratios << "   total: " << n_total_pairs << endl;

    path->finished_ = true;
    return Query();
  }

  if(args_->smc_particles() == 1) {
    *chosen_lratio = potential_merges[imax].first;
    return potential_merges[imax].second;
  } else {
    pair<double, Query> *chosen_qmerge = ChooseRandomMerge(potential_merges, rgen);
    *chosen_lratio = chosen_qmerge->first;
    return chosen_qmerge->second;
  }
}

// ----------------------------------------------------------------------------------------
// perform one merge step, i.e. find the two "nearest" clusters and merge 'em (unless we're doing doing smc, in which case we choose a random merge accordingy to their respective nearnesses)
void Glomerator::Merge(ClusterPath *path, smc::rng *rgen) {
  if(path->finished_)  // already finished this <path>, but we're still iterating 'cause some of the other paths aren't finished
    return;
  double chosen_lratio;
  Query chosen_qmerge = ChooseMerge(path, rgen, &chosen_lratio);  // NOTE chosen_lratio is not set if we hamming merge
  if(path->finished_)
    return;

  assert(seq_info_.count(chosen_qmerge.name_) != 0);
  GetNaiveSeq(chosen_qmerge.name_, &chosen_qmerge.parents_);

  double last_partition_logprob(LogProbOfPartition(path->CurrentPartition()));
  Partition new_partition(path->CurrentPartition());  // note: CurrentPartition() returns a reference
  new_partition.erase(chosen_qmerge.parents_.first);
  new_partition.erase(chosen_qmerge.parents_.second);
  new_partition.insert(chosen_qmerge.name_);
  path->AddPartition(new_partition, LogProbOfPartition(new_partition), args_->max_logprob_drop());

  if(args_->debug()) {
    printf("       merged %-8.2f", chosen_lratio);
    double newdelta = LogProbOfPartition(new_partition) - last_partition_logprob;
    if(fabs(newdelta - chosen_lratio) > 1e-8)
      printf(" ( %-20.15f != %-20.15f)", chosen_lratio, LogProbOfPartition(new_partition) - last_partition_logprob);
    printf("   %s and %s\n", chosen_qmerge.parents_.first.c_str(), chosen_qmerge.parents_.second.c_str());
    string extrastr("current (logweight " + to_string(path->CurrentLogWeight()) + ")");
  }
}

// NOTE don't remove these (yet, at least)
// // ----------------------------------------------------------------------------------------
// ClusterPair Glomerator::GetClustersToMergeForNaiveSeqGlomerate(set<vector<string> > &clusters, int max_per_cluster, bool merge_whatever_you_got) {
//   double smallest_min_distance(9999);
//   ClusterPair clusters_to_merge;
//   int n_skipped(0);
//   for(set<vector<string> >::iterator clust_a = clusters.begin(); clust_a != clusters.end(); ++clust_a) {
//     for(set<vector<string> >::iterator clust_b = clust_a; ++clust_b != clusters.end();) {
//       if(!merge_whatever_you_got && clust_a->size() + clust_b->size() > (size_t)max_per_cluster) {  // merged cluster would be too big, so look for smaller (albeit further-apart) things to merge
// 	++n_skipped;
// 	continue;
//       }

//       double min_distance(9999);  // find the smallest hamming distance between any two sequences in the two clusters
//       for(auto &query_a : *clust_a) {
// 	for(auto &query_b : *clust_b) {
// 	  string key(query_a + "-" + query_b);
// 	  double hfrac;
//        throw runtime_error("needs updating -- should use naive_hamming_fractions_ or naive_hfracs_ like everybody else");
// 	  if(hamming_fractions_.count(key) == 0) {
// 	    hfrac = NaiveHfrac(query_a, query_b);
// 	    hamming_fractions_[key] = hfrac;
// 	    hamming_fractions_[query_b + "-" + query_a] = hfrac;  // also add the reverse-ordered key
// 	  } else {
// 	    hfrac = hamming_fractions_[key];
// 	  }
// 	  if(hfrac < min_distance)
// 	    min_distance = hfrac;
// 	}
//       }

//       if(min_distance < smallest_min_distance) {
// 	smallest_min_distance = min_distance;
// 	// vector<string> &ref_a(*clust_a), &ref_b(*clust_b);
// 	clusters_to_merge = ClusterPair(*clust_a, *clust_b);
//       }
//       // if(args_->debug() && n_skipped > 0)
//       // 	printf("      skipped: %d", n_skipped);
//     }
//   }
//   return clusters_to_merge;
// }

// // ----------------------------------------------------------------------------------------
// void Glomerator::PrintClusterSizes(set<vector<string> > &clusters) {
//   for(auto &cl : clusters)
//     cout << " " << cl.size();
//   cout << endl;
// }

// // ----------------------------------------------------------------------------------------
// ClusterPair Glomerator::GetSmallBigClusters(set<vector<string> > &clusters) { // return the samllest and biggest clusters
//   ClusterPair smallbig;
//   for(auto &clust : clusters) {
//     if(smallbig.first.size() == 0 || clust.size() < smallbig.first.size())
//       smallbig.first = clust;
//     if(smallbig.second.size() == 0 || clust.size() > smallbig.second.size())
//       smallbig.second = clust;
//   }
//   return smallbig;
// }

// // ----------------------------------------------------------------------------------------
// void Glomerator::NaiveSeqGlomerate(int n_clusters) {
//   clock_t run_start(clock());
//   // clusters = [[names,] for names in naive_seqs.keys()]
//   double seqs_per_cluster = double(seq_info_.size()) / n_clusters;
//   int max_per_cluster = ceil(seqs_per_cluster);
//   if(args_->debug())
//     printf("  making %d clusters (max %d per cluster)\n", n_clusters, max_per_cluster);

//   set<vector<string> > clusters;
//   for(auto &kv : seq_info_)
//     clusters.insert(vector<string>{kv.first});

//   if(args_->debug())
//     PrintClusterSizes(clusters);

//   bool merge_whatever_you_got(false);
//   while(clusters.size() > (size_t)n_clusters) {
//     // if(args_->debug())//   printf'    current ', ' '.join([str(len(cl)) for cl in clusters])
//     ClusterPair clusters_to_merge = GetClustersToMergeForNaiveSeqGlomerate(clusters, max_per_cluster, merge_whatever_you_got);
//     if(clusters_to_merge.first.size() == 0) {  // if we didn't find a suitable pair
//       // if debug://     print '    didn\'t find shiznitz'
//       merge_whatever_you_got = true;  // next time through, merge whatever's best regardless of size
//     } else {
//       // if debug:print '    merging', len(clusters_to_merge[0]), len(clusters_to_merge[1])
//       vector<string> new_cluster(clusters_to_merge.first);
//       new_cluster.insert(new_cluster.end(), clusters_to_merge.second.begin(), clusters_to_merge.second.end());  // add clusters from the second cluster to the first cluster
//       clusters.insert(new_cluster);
//       clusters.erase(clusters_to_merge.first);  // then erase the old clusters
//       clusters.erase(clusters_to_merge.second);
//       if(args_->debug())
// 	PrintClusterSizes(clusters);
//     }
//   }

//   size_t itries(0);
//   ClusterPair smallbig = GetSmallBigClusters(clusters);
//   while(float(smallbig.second.size()) / smallbig.first.size() > 1.1 && smallbig.second.size() - smallbig.first.size() > 3) {  // keep homogenizing while biggest cluster is more than 3/2 the size of the smallest (and while their sizes differ by more than 2)
//     if(args_->debug())
//       cout << "homogenizing" << endl;
//     int n_to_keep_in_biggest_cluster = ceil(double(smallbig.first.size() + smallbig.second.size()) / 2);
//     clusters.erase(smallbig.first);
//     clusters.erase(smallbig.second);
//     smallbig.first.insert(smallbig.first.end(), smallbig.second.begin() + n_to_keep_in_biggest_cluster, smallbig.second.end());
//     smallbig.second.erase(smallbig.second.begin() + n_to_keep_in_biggest_cluster, smallbig.second.end());
//     clusters.insert(smallbig.first);
//     clusters.insert(smallbig.second);
//     if(args_->debug())
//       PrintClusterSizes(clusters);
//     ++itries;
//     if(itries > clusters.size()) {
//       // if debug: print '  too many homogenization tries'
//       break;
//     }
//     smallbig = GetSmallBigClusters(clusters);
//   }

//   // cout << "  final bcrham divvy" << endl;
//   // int tmpic(0);
//   // for(auto &clust : clusters) {
//   //   cout << "      " << tmpic << endl;
//   //   for(auto &query : clust)
//   //     cout << "          " << query << endl;
//   //   ++tmpic;
//   // }

//   ofs_.open(args_->outfile());
//   ofs_ << "partition" << endl;
//   int ic(0);
//   for(auto &clust : clusters) {
//     if(ic > 0)
//       ofs_ << "|";
//     int iq(0);
//     for(auto &query : clust) {
//       if(iq > 0)
// 	ofs_ << ";";
//       ofs_ << query;
//       ++iq;
//     }
//     ++ic;
//   }
//   ofs_ << endl;
//   ofs_.close();

//   cout << "        naive hamming cluster time " << ((clock() - run_start) / (double)CLOCKS_PER_SEC) << "   to condense " << seq_info_.size() << " --> " << n_clusters <<  endl;
// }

}
