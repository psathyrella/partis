#include "glomerator.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track) :
  track_(track),
  args_(args),
  vtb_dph_("viterbi", args_, gl, hmms),
  fwd_dph_("forward", args_, gl, hmms),
  n_fwd_calculated_(0),
  n_vtb_calculated_(0),
  n_hfrac_calculated_(0),
  n_hfrac_merges_(0),
  n_lratio_merges_(0),
  asym_factor_(4.),
  current_partition_(nullptr),
  progress_file_(fopen((args_->outfile() + ".progress").c_str(), "w"))
{
  time(&last_status_write_time_);
  ReadCacheFile();

  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key = SeqNameStr(qry_seq_list[iqry], ":");

    initial_partition_.insert(key);
    seq_info_[key] = qry_seq_list[iqry];
    seed_missing_[key] = !InString(args_->seed_unique_id(), key);
    only_genes_[key] = args_->str_lists_["only_genes"][iqry];
    mute_freqs_[key] = avgVector(args_->float_lists_["mut_freqs"][iqry]);

    KSet kmin(args_->integers_["k_v_min"][iqry], args_->integers_["k_d_min"][iqry]);
    KSet kmax(args_->integers_["k_v_max"][iqry], args_->integers_["k_d_max"][iqry]);
    KBounds kb(kmin, kmax);
    kbinfo_[key] = kb;
  }

  current_partition_ = &initial_partition_;
}

// ----------------------------------------------------------------------------------------
Glomerator::~Glomerator() {
  cout << FinalString() << endl;
  if(args_->cachefile() != "")
    WriteCachedLogProbs();
  fclose(progress_file_);
  remove((args_->outfile() + ".progress").c_str());
}

// ----------------------------------------------------------------------------------------
void Glomerator::CacheNaiveSeqs() {  // they're written to file in the destructor, so we just need to calculate them here
  cout << "      caching all naive sequences" << endl;
  for(auto &kv : seq_info_)
    GetNaiveSeq(kv.first);
  ofs_.open(args_->outfile());  // a.t.m. I'm signalling that I finished ok by doing this
  ofs_.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  if(args_->debug()) {
    cout << "   hieragloming " << initial_partition_.size() << " clusters";
    if(args_->seed_unique_id() != "")
      cout << "  (" << GetSeededClusters(initial_partition_).size() << " seeded)";
    cout << endl;
  }

  if(args_->logprob_ratio_threshold() == -INFINITY)
    throw runtime_error("logprob ratio threshold not specified");

  ClusterPath cp(initial_partition_);
  do {
    Merge(&cp);
  } while(!cp.finished_);

  vector<ClusterPath> paths{cp};
  WritePartitions(paths);
  if(args_->annotationfile() != "")
    WriteAnnotations(paths);
}

// ----------------------------------------------------------------------------------------
Partition Glomerator::GetAnInitialPartition(int &initial_path_index, double &logweight) {
  assert(0);
  // initial_path_index = i_initial_partition_;
  // assert(i_initial_partition_ < (int)initial_partitions_.size());
  // logweight = initial_logweights_[i_initial_partition_];
  // return initial_partitions_.at(i_initial_partition_++);
  return Partition();
}

// ----------------------------------------------------------------------------------------
void Glomerator::ReadCacheFile() {
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
  assert(headstrs[4].find("errors") == 0);

  // NOTE there can be two lines with the same key (say if in one run we calculated the naive seq, and in a later run calculated the log prob)
  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> column_list = SplitString(line, ",");
    assert(column_list.size() == 5);
    string query(column_list[0]);
    string errors(column_list[4]);
    if(errors.find("no_path") != string::npos) {
      failed_queries_.insert(query);
      continue;
    }

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
  if(args_->debug())
    cout << "        writing partitions" << endl;
  ofs_.open(args_->outfile());
  ofs_ << setprecision(20);
  ofs_ << "partition,logprob" << endl;
  int ipath(0);
  for(auto &cp : paths) {
    unsigned istart(0);
    if((int)cp.partitions().size() > args_->n_partitions_to_write())
      istart = cp.partitions().size() - args_->n_partitions_to_write();
    for(unsigned ipart=istart; ipart<cp.partitions().size(); ++ipart) {
      if(args_->write_logprob_for_each_partition())  // only want to calculate this the last time through, i.e. when we're only one process
	cp.set_logprob(ipart, LogProbOfPartition(cp.partitions()[ipart]));
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
void Glomerator::WriteAnnotations(vector<ClusterPath> &paths) {
  cout << "      calculating and writing annotations" << endl;
  ofstream annotation_ofs;
  annotation_ofs.open(args_->annotationfile());
  StreamHeader(annotation_ofs, "viterbi");

  assert(paths.size() == 1);  // would need to update this for smc
  int ipath(0);
  ClusterPath cp(paths[ipath]);
  unsigned ipart(cp.partitions().size() - 1);  // just write the last (best) one. NOTE that we're no longer keeping track of the total log prob of the partition, since a) a bunch (most?) of the time we're merging with naive hfrac and b) we don't need it, anyway
  for(auto &cluster : cp.partitions()[ipart]) {
    if(args_->seed_unique_id() != "" && SeedMissing(cluster))
      continue;

    RecoEvent event;
    CalculateNaiveSeq(cluster, &event);  // calculate the viterbi path from scratch -- we're doing so much translation crap at this point it's just too hard to keep track of things otherwise

    if(event.genes_["d"] == "") {  // shouldn't happen any more, but it is a check that could fail at some point
      cout << "WTF " << cluster << " x" << event.naive_seq_ << "x" << endl;
      assert(0);
    }
    StreamViterbiOutput(annotation_ofs, event, seq_info_[cluster], "");
  }
  annotation_ofs.close();
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(Partition &partition, bool debug) {

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
  printf("    %-8.2f %s partition\n", -INFINITY/*LogProbOfPartition(partition)*/, extra_cstr);
  for(auto &key : partition)
    cout << "          " << key << endl;
}

// ----------------------------------------------------------------------------------------
string Glomerator::CacheSizeString() {
  char buffer[2000];
  sprintf(buffer, "      %8zu   %8zu   %8zu   %8zu   %8zu", log_probs_.size(), naive_hfracs_.size(), lratios_.size(), naive_seqs_.size(), errors_.size());
  return string(buffer);
}

// ----------------------------------------------------------------------------------------
string Glomerator::FinalString() {
    char buffer[2000];
    sprintf(buffer, "        calcd:   vtb %-4d  fwd %-4d  hfrac %-8d    merged:  hfrac %-4d lratio %-4d", n_vtb_calculated_, n_fwd_calculated_, n_hfrac_calculated_, n_hfrac_merges_, n_lratio_merges_);
    return string(buffer);
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteStatus() {

  // // ----------------------------------------------------------------------------------------
  // // print memory usage
  // ifstream smapfs("/proc/self/smaps");
  // string tmpline;
  // cout << "contents of /proc/self/smaps:" << endl;
  // while(getline(smapfs, tmpline)) {
  //   cout << "MEM " << tmpline << endl;
  // }
  // cout << endl;
  // // ----------------------------------------------------------------------------------------

  // cout << CacheSizeString() << endl;
  time_t current_time;
  time(&current_time);
  if(difftime(current_time, last_status_write_time_) > 300) {  // write something every five minutes
    char buffer[2000];
    strftime(buffer, 2000, "%b %d %T", localtime(&current_time));  // %H:%M
    // fprintf(progress_file_, "      %s    %4d clusters    calcd  fwd %-4d   vtb %-4d   hfrac %-8d    merged  hfrac %-4d\n", buffer, (int)path_->CurrentPartition().size(), n_fwd_calculated_, n_vtb_calculated_);
    fprintf(progress_file_, "      %s    %4d clusters", buffer, (int)current_partition_->size());
    fprintf(progress_file_, "   %s", FinalString().c_str());

    fprintf(progress_file_, "     %s\n", ClusterSizeString(current_partition_).c_str());

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
string Glomerator::ClusterSizeString(Partition *partition) {
  vector<int> cluster_sizes;
  for(auto &cluster : *partition) {
    cluster_sizes.push_back(CountMembers(cluster));
  }
  sort(cluster_sizes.begin(), cluster_sizes.end());
  reverse(cluster_sizes.begin(), cluster_sizes.end());
  string return_str("          clusters: ");
  int n_singletons(0);
  for(size_t is=0; is<cluster_sizes.size(); ++is) {
    if(cluster_sizes[is] == 1)
      ++n_singletons;
    else
      return_str += " "  + to_string(cluster_sizes[is]);
  }
  return return_str + "  (" + to_string(n_singletons) + " singletons)";
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
string Glomerator::PrintStr(string queries) {
  if(CountMembers(queries) < 10)
    return queries;
  else
    return "len(" + to_string(CountMembers(queries)) + ")";
}

// ----------------------------------------------------------------------------------------
bool Glomerator::SeedMissing(string queries, string delimiter) {
  return seed_missing_[queries];  // NOTE after refactoring the double loops, we probably don't really need to cache all these any more
  // set<string> queryset(SplitString(queries, delimiter));  // might be faster to look for :uid: and uid: and... hm, wait, that's kind of hard
  // return !InString(args_->seed_unique_id(), queries,  delimiter);
  //// oh, wait this (below) won't work without more checks
  // if(queries.find(delimiter + uid) != string::npos)
  //   return true;
  // else if(queries.find(uid + delimiter) != string::npos)
  //   return true;
  // else
  //   return false;
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
  double hfrac(INFINITY);
  if(failed_queries_.count(key_a) || failed_queries_.count(key_b))
    return hfrac;
  naive_hfracs_[joint_key] = CalculateHfrac(seq_a, seq_b);

  return naive_hfracs_[joint_key];
}

// ----------------------------------------------------------------------------------------
string Glomerator::ChooseSubsetOfNames(string queries, int n_max) {
  if(name_subsets_.count(queries))
    return name_subsets_[queries];

  assert(seq_info_.count(queries));
  vector<string> namevector(SplitString(queries, ":"));

  srand(hash<string>{}(queries));  // make sure we get the same subset each time we pass in the same queries (well, if there's different thresholds for naive_seqs annd logprobs they'll each get their own [very correlated] subset)

  // first decide which indices we'll choose
  set<int> ichosen;
  vector<int> ichosen_vec;  // don't really need both of these... but maybe it's faster
  set<string> chosen_strs;  // make sure we don't choose seed unique id more than once
  for(size_t iname=0; iname<unsigned(n_max); ++iname) {
    int ich(-1);
    int n_tries(0);
    while(ich < 0 || ichosen.count(ich) || chosen_strs.count(namevector[ich])) {
      ich = rand() % namevector.size();
      ++n_tries;
      if(n_tries > 1e6)
	throw runtime_error("too many tries in Glomerator::ChooseSubsetOfNames() -- maybe too many copies of the seed unique id?");
    }
    ichosen.insert(ich);
    chosen_strs.insert(namevector[ich]);
    ichosen_vec.push_back(ich);
  }

  // then sort 'em
  sort(ichosen_vec.begin(), ichosen_vec.end());

  // and finally make the new vectors
  vector<string> subqueryvec;
  vector<Sequence> subseqs;
  for(auto &ich : ichosen_vec) {
    subqueryvec.push_back(namevector[ich]);
    subseqs.push_back(seq_info_[queries][ich]);
  }

  string subqueries(JoinStrings(subqueryvec));

  seq_info_[subqueries] = subseqs;
  seed_missing_[subqueries] = !InString(args_->seed_unique_id(), subqueries);
  kbinfo_[subqueries] = kbinfo_[queries];  // just use the entire/super cluster for this stuff. It's just overly conservative (as long as you keep the mute freqs the same)
  mute_freqs_[subqueries] = mute_freqs_[queries];
  only_genes_[subqueries] = only_genes_[queries];

  if(args_->debug())
    cout << "                chose subset  " << queries << "  -->  " << subqueries << endl;

  name_subsets_[queries] = subqueries;
  return subqueries;
}

// ----------------------------------------------------------------------------------------
string Glomerator::GetNaiveSeqNameToCalculate(string actual_queries) {
  // NOTE we don't really need to cache the names like this, since we're setting the random seed when we choose a subset. But it just seems so messy to go through the whole subset calculation every time, even though I profiled it and it's not a significant contributor
  if(naive_seq_name_translations_.count(actual_queries))
    return naive_seq_name_translations_[actual_queries];

  // if cluster is less than half again larger than N, just use <actual_queries>
  if(CountMembers(actual_queries) < 1.5 * args_->biggest_naive_seq_cluster_to_calculate())  // if <<actual_queries>> is small return all of 'em
    return actual_queries;

  // but if it's bigger than this, replace it with a subset of size N
  string subqueries = ChooseSubsetOfNames(actual_queries, args_->biggest_naive_seq_cluster_to_calculate());
  if(args_->debug() > 0)
    cout << "                translate for naive seq  " << actual_queries << "  -->  " << subqueries << endl;

  naive_seq_name_translations_[actual_queries] = subqueries;
  return subqueries;
}

// ----------------------------------------------------------------------------------------
string Glomerator::GetLogProbNameToCalculate(string queries, int n_max) {
  string queries_to_calc(queries);
  if(logprob_asymetric_translations_.count(queries)) {
    if(args_->debug())
      cout << "             using asymetric translation  " << queries << "  -->  " << logprob_asymetric_translations_[queries] << endl;
    queries_to_calc = logprob_asymetric_translations_[queries];
  } 

  if(CountMembers(queries_to_calc) > n_max) {
    return ChooseSubsetOfNames(queries_to_calc, n_max);
  } else {
    return queries_to_calc;
  }
}

// ----------------------------------------------------------------------------------------
pair<string, string> Glomerator::GetLogProbPairOfNamesToCalculate(string actual_queries, pair<string, string> actual_parents) {
  // NOTE we don't really need to cache the names like this, since we're setting the random seed when we choose a subset. But it just seems so messy to go through the whole subset calculation every time, even though I profiled it and it's not a significant contributor
  if(logprob_name_translations_.count(actual_queries))
    return logprob_name_translations_[actual_queries];

  int n_max(args_->biggest_logprob_cluster_to_calculate());

  // if the merged cluster is small, just use the whole thing
  if(CountMembers(actual_queries) <= 2. * n_max)  // TODO don't hard code this
    return actual_parents;

  // replace either/both parents as necessary
  pair<string, string> queries_to_calc;
  queries_to_calc.first = GetLogProbNameToCalculate(actual_parents.first, n_max);  // note, no factor of 1.5, "since" this is more considering the lratio as a whole, and
  queries_to_calc.second = GetLogProbNameToCalculate(actual_parents.second, n_max);

  if(args_->debug())
    printf("                translate for lratio (%s)   %s  %s  -->  %s  %s\n", actual_queries.c_str(), actual_parents.first.c_str(), actual_parents.second.c_str(), queries_to_calc.first.c_str(), queries_to_calc.second.c_str());

  logprob_name_translations_[actual_queries] = queries_to_calc;
  return queries_to_calc;
}

// ----------------------------------------------------------------------------------------
bool Glomerator::FirstParentMuchBigger(string queries, string queries_other, int nmax) {
  int nseq(CountMembers(queries));
  int nseq_other(CountMembers(queries_other));
  if(nseq > nmax && float(nseq) / nseq_other > asym_factor_ ) {  // if <nseq> is large, and if <nseq> more than twice the size of <nseq_other>, use the existing name translation (for which we should already have a logprob and a naive seq)
    if(args_->debug()) {
      cout << "                asymetric  " << nseq << " " << nseq_other << "  use " << queries << "  instead of " << JoinNames(queries, queries_other) << endl;
      if(naive_seq_name_translations_.count(queries))
	cout << "                    naive seq translates to " << naive_seq_name_translations_[queries] << endl;
    }
    return true;
  }
  return false;
}

// ----------------------------------------------------------------------------------------
string Glomerator::FindNaiveSeqNameReplace(pair<string, string> *parents) {
  assert(parents != nullptr);

  // if both parents have the same naive sequence, just use the first one
  if(GetNaiveSeq(parents->first) == GetNaiveSeq(parents->second))
    return parents->first;

  // and if one of the parents is much bigger, just use that one
  int nmax = 1.5 * args_->biggest_naive_seq_cluster_to_calculate();  // TODO don't hard code the factor
  if(FirstParentMuchBigger(parents->first, parents->second, nmax))
    return parents->first;
  if(FirstParentMuchBigger(parents->second, parents->first, nmax))
    return parents->second;

  return string("");  // if we fall through, we don't want to replace the current query (but maybe we'll later decide to only use a subset of it)
}

// ----------------------------------------------------------------------------------------
string &Glomerator::GetNaiveSeq(string queries, pair<string, string> *parents) {
  if(naive_seqs_.count(queries))
    return naive_seqs_[queries];

  // see if we want to just straight up use the naive sequence from one of the parents
  if(parents != nullptr) {
    string name_with_which_to_replace = FindNaiveSeqNameReplace(parents);
    if(name_with_which_to_replace != "") {
      naive_seqs_[queries] = GetNaiveSeq(name_with_which_to_replace);  // copy the whole sequence object  TODO this doesn't follow/do the turtle thing
      return naive_seqs_[queries];
    }
  }

  // see if we want to calculate with only a subset of the queries
  string queries_to_calc = GetNaiveSeqNameToCalculate(queries);

  // actually calculate the viterbi path for whatever queries we've decided on
  if(naive_seqs_.count(queries_to_calc) == 0)
    naive_seqs_[queries_to_calc] = CalculateNaiveSeq(queries_to_calc);

  // if we did some translation, propagate the naive sequence back to the queries we were originally interested in
  if(queries_to_calc != queries)
    naive_seqs_[queries] = naive_seqs_[queries_to_calc];

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
  // NOTE we could avoid recalculating a lot of the denominators if we didn't randomly choose a subset, and instead looked to see what we already have (but then it would be a lot harder to have a representive sample...)

  string joint_name(JoinNames(key_a, key_b));

  if(lratios_.count(joint_name))  // NOTE as in other places, this assumes there's only *one* way to get to a given joint name (or at least that we'll get about the same answer each different way)
    return lratios_[joint_name];

  Query full_qmerged = GetMergedQuery(key_a, key_b);  // have to make this to get seq_info_ and whatnot filled up
  pair<string, string> parents_to_calc = GetLogProbPairOfNamesToCalculate(joint_name, full_qmerged.parents_);
  string key_a_to_calc = parents_to_calc.first;
  string key_b_to_calc = parents_to_calc.second;
  Query qmerged_to_calc = GetMergedQuery(key_a_to_calc, key_b_to_calc);  // NOTE also enters the merged query's info into seq_info_, kbinfo_, mute_freqs_, and only_genes_

  double log_prob_a = GetLogProb(key_a_to_calc);
  double log_prob_b = GetLogProb(key_b_to_calc);
  double log_prob_ab = GetLogProb(qmerged_to_calc.name_);

  double lratio(log_prob_ab - log_prob_a - log_prob_b);
  if(args_->debug()) {
    printf("             %8.3f =", lratio);
    printf(" %s - %s - %s", joint_name.c_str(), key_a.c_str(), key_b.c_str());
    if(qmerged_to_calc.name_ != joint_name || key_a_to_calc != key_a || key_b_to_calc != key_b)
      printf(" (calcd  %s - %s - %s)", qmerged_to_calc.name_.c_str(), key_a_to_calc.c_str(), key_b_to_calc.c_str());
    printf("\n");
  }

  return lratios_[joint_name] = lratio;
  return lratio;
}

// ----------------------------------------------------------------------------------------
string Glomerator::CalculateNaiveSeq(string queries, RecoEvent *event) {
  if(event == nullptr)  // if we're calling it with <event> set, then we know we're recalculating some things
    assert(naive_seqs_.count(queries) == 0);

  if(seq_info_.count(queries) == 0 || kbinfo_.count(queries) == 0 || only_genes_.count(queries) == 0 || mute_freqs_.count(queries) == 0)
    throw runtime_error("no info for " + queries);

  ++n_vtb_calculated_;

  Result result(kbinfo_[queries], args_->chain());
  bool stop(false);
  do {
    // assert(SameLength(seq_info_[queries], true));
    result = vtb_dph_.Run(seq_info_[queries], kbinfo_[queries], only_genes_[queries], mute_freqs_[queries]);  // NOTE the sequences in <seq_info_[queries]> should already be the same length, since they've already been merged
    kbinfo_[queries] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand() || result.no_path_;  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.no_path_) {
    AddFailedQuery(queries, "no_path");
    return "";
  }
  if(result.boundary_error())
    errors_[queries] = errors_[queries] + ":boundary";

  if(event != nullptr)
    *event = result.best_event();

  WriteStatus();
  return result.best_event().naive_seq_;
}

// ----------------------------------------------------------------------------------------
double Glomerator::CalculateLogProb(string queries) {  // NOTE can modify kbinfo_
  // NOTE do *not* call this from anywhere except GetLogProb()
  assert(log_probs_.count(queries) == 0);

  if(seq_info_.count(queries) == 0 || kbinfo_.count(queries) == 0 || only_genes_.count(queries) == 0 || mute_freqs_.count(queries) == 0)
    throw runtime_error("no info for " + queries);
  
  ++n_fwd_calculated_;

  Result result(kbinfo_[queries], args_->chain());
  bool stop(false);
  do {
    result = fwd_dph_.Run(seq_info_[queries], kbinfo_[queries], only_genes_[queries], mute_freqs_[queries]);  // NOTE <only_genes> isn't necessarily <only_genes_[queries]>, since for the denominator calculation we take the OR
    kbinfo_[queries] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand() || result.no_path_;  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "             expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.no_path_) {
    AddFailedQuery(queries, "no_path");
    return -INFINITY;
  }
  if(result.boundary_error() && !result.could_not_expand())  // could_not_expand means the max is at the edge of the sequence -- e.g. k_d min is 1
    errors_[queries] = errors_[queries] + ":boundary";

  WriteStatus();
  return result.total_score();
}

// ----------------------------------------------------------------------------------------
void Glomerator::AddFailedQuery(string queries, string error_str) {
    errors_[queries] = errors_[queries] + ":" + error_str;
    failed_queries_.insert(queries);
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
    // if(all_names.count(name)) {
    //   if(args_->seed_unique_id() == "" || !InString(args_->seed_unique_id(), name))  // if seed id isn't set, or if it is set but it isn't in <name>
    // 	throw runtime_error("tried to add sequence with name " + name + " twice in Glomerator::MergeSeqVectors()\n    when merging " + name_a + " and " + name_b);
    // } else {
      all_names.insert(name);
    // }
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
  seed_missing_[qmerged.name_] = !InString(args_->seed_unique_id(), qmerged.name_);
  kbinfo_[qmerged.name_] = qmerged.kbounds_;
  mute_freqs_[qmerged.name_] = qmerged.mean_mute_freq_;
  only_genes_[qmerged.name_] = qmerged.only_genes_;
  
  return qmerged;
}

// ----------------------------------------------------------------------------------------
pair<double, Query> *Glomerator::ChooseRandomMerge(vector<pair<double, Query> > &potential_merges) {
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
  double drawpoint = 0.;  assert(0); //rgen->Uniform(0., 1.);
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
Partition Glomerator::GetSeededClusters(Partition &partition) {
  Partition clusters;
  for(auto &queries : partition) {
    if(!SeedMissing(queries))
      clusters.insert(queries);
  }
  return clusters;
}

// ----------------------------------------------------------------------------------------
pair<double, Query> Glomerator::FindHfracMerge(ClusterPath *path) {
  double min_hamming_fraction(INFINITY);
  Query min_hamming_merge;

  Partition outer_clusters(path->CurrentPartition());  // for plain partitioning, outer loop is over everything in the current partition
  if(args_->seed_unique_id() != "")  // whereas if seed unique id is set, outer loop is only over those clusters that contain the seed
    outer_clusters = GetSeededClusters(path->CurrentPartition());
  for(Partition::iterator it_a = outer_clusters.begin(); it_a != outer_clusters.end(); ++it_a) {
    Partition::iterator it_b(it_a);  // for plain partitioning, inner loop starts at the position of the outer iterator + 1
    ++it_b;
    if(args_->seed_unique_id() != "")  // but if the seed's set, inner loop is over the *entire* current partition (including the seeded clusters). So we also need to skip key_a == key_b below.
      it_b = path->CurrentPartition().begin();

    for( ; true; ++it_b) {
      if(args_->seed_unique_id() == "") {
	if(it_b == outer_clusters.end())  // NOTE you can't replace this with path->CurrentPartition().end() because in this case <it_b> is looping over the *copy* of path->CurrentPartition()
	  break;
      } else {
	if(it_b == path->CurrentPartition().end())
	  break;
      }

      string key_a(*it_a), key_b(*it_b);
      if(key_a == key_b)  // otherwise we'd loop over the seeded ones twice
	continue;
      if(failed_queries_.count(key_a) || failed_queries_.count(key_b))
	continue;

      double hfrac = NaiveHfrac(key_a, key_b);
      if(hfrac > args_->hamming_fraction_bound_hi())  // if naive hamming fraction too big, don't even consider merging the pair
	continue;

      if(args_->hamming_fraction_bound_lo() <= 0.0 || hfrac >= args_->hamming_fraction_bound_lo())
	continue;

      if(hfrac < min_hamming_fraction) {
	  min_hamming_fraction = hfrac;
	  min_hamming_merge = GetMergedQuery(key_a, key_b);  // NOTE also enters the merged query's info into seq_info_, kbinfo_, mute_freqs_, and only_genes_
      }
    }
  }

  if(min_hamming_fraction != INFINITY) {
    ++n_hfrac_merges_;
    if(args_->debug())
      printf("          hfrac merge  %.3f   %s  %s\n", min_hamming_fraction, PrintStr(min_hamming_merge.parents_.first).c_str(), PrintStr(min_hamming_merge.parents_.second).c_str());
  }

  return pair<double, Query>(min_hamming_fraction, min_hamming_merge);
}

// ----------------------------------------------------------------------------------------
pair<double, Query> Glomerator::FindLRatioMerge(ClusterPath *path) {
  double max_lratio(-INFINITY);
  Query chosen_qmerge;
  int n_total_pairs(0), n_skipped_hamming(0), n_small_lratios(0);

  Partition outer_clusters(path->CurrentPartition());
  if(args_->seed_unique_id() != "")  // see comments in FindHfracMerge
    outer_clusters = GetSeededClusters(path->CurrentPartition());
  for(Partition::iterator it_a = outer_clusters.begin(); it_a != outer_clusters.end(); ++it_a) {
    Partition::iterator it_b(it_a);
    ++it_b;
    if(args_->seed_unique_id() != "")
      it_b = path->CurrentPartition().begin();

    for( ; true; ++it_b) {
      if(args_->seed_unique_id() == "") {
	if(it_b == outer_clusters.end())  // NOTE you can't replace this with path->CurrentPartition().end() because in this case <it_b> is looping over the *copy* of path->CurrentPartition()
	  break;
      } else {
	if(it_b == path->CurrentPartition().end())
	  break;
      }

      string key_a(*it_a), key_b(*it_b);
      if(key_a == key_b)  // otherwise we'd loop over the seeded ones twice
	continue;
      if(failed_queries_.count(key_a) || failed_queries_.count(key_b))
	continue;

      ++n_total_pairs;

      double hfrac = NaiveHfrac(key_a, key_b);
      if(hfrac > args_->hamming_fraction_bound_hi()) {  // if naive hamming fraction too big, don't even consider merging the pair
	++n_skipped_hamming;
	continue;
      }

      double lratio = GetLogProbRatio(key_a, key_b);

      // don't merge if lratio is small (less than zero, more or less)
      if(LikelihoodRatioTooSmall(lratio, CountMembers(key_a) + CountMembers(key_b))) {
	++n_small_lratios;
	continue;
      }

      if(lratio > max_lratio) {
	max_lratio = lratio;
	chosen_qmerge = GetMergedQuery(key_a, key_b);
      }
    }
  }

  if(max_lratio == -INFINITY) {  // if we didn't find any merges that we liked
    path->finished_ = true;
    cout << "        stop with:  big hfrac " << n_skipped_hamming << "   small lratio " << n_small_lratios << "   total " << n_total_pairs;
    if(failed_queries_.size() > 0)
      cout << "   (" << failed_queries_.size() << " failed queries)";
    cout << endl;
  } else {
    ++n_lratio_merges_;
    if(args_->debug())
      printf("          lratio merge  %.3f   %s  %s\n", max_lratio, PrintStr(chosen_qmerge.parents_.first).c_str(), PrintStr(chosen_qmerge.parents_.second).c_str());
  }

  return pair<double, Query>(max_lratio, chosen_qmerge);
}

// ----------------------------------------------------------------------------------------
void Glomerator::UpdateLogProbTranslationsForAsymetrics(Query &qmerge) {
  string queries("");

  // see if one of the parents is much bigger than the other
  int nmax = 1.5 * args_->biggest_logprob_cluster_to_calculate();  // TODO don't hard code the factor
  if(FirstParentMuchBigger(qmerge.parents_.first, qmerge.parents_.second, nmax))
    queries = qmerge.parents_.first;
  else if(FirstParentMuchBigger(qmerge.parents_.second, qmerge.parents_.first, nmax))
    queries = qmerge.parents_.second;

  if(queries != "") {

    // if the large parent itself was formed by an asymetric merge, keep following the chain of translations (note that since we do this every time, the chain can't get longer than 1 [erm, I think])
    string subqueries(queries);
    while(logprob_asymetric_translations_.count(subqueries)) {
      if(args_->debug())
	cout << "                  turtles " << subqueries << "  -->  " << logprob_asymetric_translations_[subqueries] << endl;
      subqueries = logprob_asymetric_translations_[subqueries];
    }

    // if we haven't added too many new sequences since we last calculated something, we can just reuse things  TODO don't hard code this factor
    if(float(CountMembers(queries)) / CountMembers(subqueries) < 2.)  {
      if(args_->debug())
	cout << "                logprob asymetric translation  " << qmerge.name_ << "  -->  " << subqueries << endl;
      logprob_asymetric_translations_[qmerge.name_] = subqueries;  // note that this just says *if* we need this logprob in the future, we should instead calculate this other one -- but we may never actually need it
    } else {
      if(args_->debug())
	cout << "                  ratio too big for asymetric " << CountMembers(queries) << " " << CountMembers(subqueries) << endl;
    }
  }
}

// ----------------------------------------------------------------------------------------
// perform one merge step, i.e. find the two "nearest" clusters and merge 'em (unless we're doing doing smc, in which case we choose a random merge accordingy to their respective nearnesses)
void Glomerator::Merge(ClusterPath *path) {
  if(path->finished_)  // already finished this <path>, but we're still iterating 'cause some of the other paths aren't finished
    return;
  pair<double, Query> qpair = FindHfracMerge(path);
  if(qpair.first == INFINITY)
    qpair = FindLRatioMerge(path);
  if(path->finished_)
    return;

  WriteStatus();
  Query chosen_qmerge = qpair.second;

  assert(seq_info_.count(chosen_qmerge.name_));
  GetNaiveSeq(chosen_qmerge.name_, &chosen_qmerge.parents_);  // this *needs* to happen here so it has the parental information
  UpdateLogProbTranslationsForAsymetrics(chosen_qmerge);

  Partition new_partition(path->CurrentPartition());
  new_partition.erase(chosen_qmerge.parents_.first);
  new_partition.erase(chosen_qmerge.parents_.second);
  new_partition.insert(chosen_qmerge.name_);
  path->AddPartition(new_partition, -INFINITY);
  current_partition_ = &path->CurrentPartition();

  if(args_->debug()) {
    printf("       merged   %s  %s\n", chosen_qmerge.parents_.first.c_str(), chosen_qmerge.parents_.second.c_str());
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
