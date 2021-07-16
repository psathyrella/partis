#include "glomerator.h"

namespace ham {

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track) :
  track_(track),
  args_(args),
  gl_(gl),
  hmms_(hmms),
  n_fwd_calculated_(0),
  n_vtb_calculated_(0),
  n_hfrac_calculated_(0),
  n_hfrac_merges_(0),
  n_lratio_merges_(0),
  asym_factor_(4.),
  force_merge_(false),
  current_partition_(nullptr),
  progress_file_(fopen((args_->outfile() + ".progress").c_str(), "w"))
{
  time(&last_status_write_time_);
  ReadCacheFile();

  for(size_t iqry = 0; iqry < qry_seq_list.size(); iqry++) {
    string key = SeqNameStr(qry_seq_list[iqry], ":");
    KSet kmin(args_->integers_["k_v_min"][iqry], args_->integers_["k_d_min"][iqry]);
    KSet kmax(args_->integers_["k_v_max"][iqry], args_->integers_["k_d_max"][iqry]);

    initial_partition_.insert(key);

    for(auto &seq_vec : qry_seq_list)
      for(auto &seq : seq_vec)
	single_seqs_[seq.name()] = seq;

    for(auto &uid : SplitString(key)) {
      single_seq_cachefo_[uid] = Query(uid,  // NOTE these are not necessarily the same as they would be (well, were) for the single seqs -- e.g. only_genes is now the OR for all the sequences
				       GetSeqs(uid),
				       !InString(args_->seed_unique_id(), uid),
				       args_->str_lists_["only_genes"][iqry],
				       KBounds(kmin, kmax),
				       args_->floats_["mut_freq"][iqry],
				       args_->integers_["cdr3_length"][iqry]);
    }

    cachefo_[key] = Query(key,
			  GetSeqs(key),
			  !InString(args_->seed_unique_id(), key),
			  args_->str_lists_["only_genes"][iqry],
			  KBounds(kmin, kmax),
			  args_->floats_["mut_freq"][iqry],
			  args_->integers_["cdr3_length"][iqry]);
  }

  current_partition_ = &initial_partition_;
}

// ----------------------------------------------------------------------------------------
Glomerator::~Glomerator() {
  cout << FinalString(true) << endl;
  WriteCacheFile();
  fclose(progress_file_);
  remove((args_->outfile() + ".progress").c_str());

// // ----------------------------------------------------------------------------------------
//   runps();
//   cout << "vtb dph" << endl;
//   vtb_dph_.Clear();
//   sleep(1);
//   runps();
//   cout << "fwd dph" << endl;
//   fwd_dph_.Clear();
//   sleep(1);
// // ----------------------------------------------------------------------------------------
}

// ----------------------------------------------------------------------------------------
void Glomerator::CacheNaiveSeqs() {  // they're written to file in the destructor, so we just need to calculate them here
  cout << "      caching all naive sequences" << endl;
  for(auto &kv : cachefo_)
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

  WritePartitions(cp);
  if(args_->annotationfile() != "")
    WriteAnnotations(cp);
}

// // ----------------------------------------------------------------------------------------
// Partition Glomerator::GetAnInitialPartition(double &logweight) {
//   assert(0);
//   return Partition();
// }

// ----------------------------------------------------------------------------------------
void Glomerator::ReadCacheFile() {
  if(args_->input_cachefname() == "") {
    cout << "        read-cache:  logprobs 0   naive-seqs 0" << endl;
    return;
  }

  ifstream ifs(args_->input_cachefname());
  if(!ifs.is_open())
    throw runtime_error("input cache file " + args_->input_cachefname() + " dne\n");

  string line;
  // check the header is right (no cached info)
  if(!getline(ifs, line)) {
    cout << "        empty cachefile" << endl;
    return;  // return for zero length file
  }
  line.erase(remove(line.begin(), line.end(), '\r'), line.end());
  vector<string> headstrs(SplitString(line, ","));
  assert(headstrs[0].find("unique_ids") == 0);  // these have to match the line in WriteCacheFile(), as well as partition_cachefile_headers in utils.py
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
    if(logprob_str.size() > 0) {  // NOTE <query> might already be in <log_probs_> (see above), but this won't replace it unless it's actually set in the file (we could also check that they're similar, but since we don't expect them to always be identical, that would be complicated)
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
  cout << "        read-cache:  logprobs " << log_probs_.size() << "   naive-seqs " << naive_seqs_.size() << endl;
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
void Glomerator::WriteCacheFile() {
  if(args_->output_cachefname() == "")
    return;

  ofstream log_prob_ofs(args_->output_cachefname());
  if(!log_prob_ofs.is_open())
    throw runtime_error("couldn't open output cache file " + args_->output_cachefname() + "\n");

  log_prob_ofs << "unique_ids,logprob,naive_seq,naive_hfrac,errors" << endl;  // these have to match the line in ReadCacheFile(), as well as partition_cachefile_headers in utils.py
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
void Glomerator::WritePartitions(ClusterPath &cp) {
  clock_t run_start(clock());
  if(args_->debug())
    cout << "        writing partitions" << endl;
  ofs_.open(args_->outfile());
  ofs_ << setprecision(20);
  ofs_ << "partition,logprob" << endl;
  unsigned istart(0);
  if((int)cp.partitions().size() > args_->n_partitions_to_write())
    istart = cp.partitions().size() - args_->n_partitions_to_write();
  for(unsigned ipart=istart; ipart<cp.partitions().size(); ++ipart) {
    if(args_->write_logprob_for_each_partition())  // only want to calculate this the last time through, i.e. when we're only one process NOTE this calculation can change the clustering (if we did an hfrac merge that logprob thinks we shouldn't have merged, when the python reads the partitions it'll notice this and choose the unmerged partition)
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
  ofs_.close();
  if(args_->write_logprob_for_each_partition())
    printf("        partition writing time (probably includes calculating a bunch of new logprobs) %.1f\n", ((clock() - run_start) / (double)CLOCKS_PER_SEC));
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteAnnotations(ClusterPath &cp) {
  cout << "DEPRECATED" << endl;  // for somewhat technical reasons -- it still basically works (see notes in partitiondriver.py) UPDATE I can't find any notes about this in partitiondriver.py, but I think the basic deal is I'm running a whole separate bcrham process (after I'm done partitioning) to get the annotations. I think one (perhaps the main?) reason was that translation is really complicated, but for the final annotations we typically want the real full cluster calculation
  clock_t run_start(clock());
  cout << "      calculating and writing annotations" << endl;
  ofstream annotation_ofs;
  annotation_ofs.open(args_->annotationfile());
  StreamHeader(annotation_ofs, "viterbi");

  // NOTE we're no longer calculating the logprob for *every* partition, but in Glomerator::WritePartitions() we *do* calculate them if we're told to (i.e. the last time through), and this can make it so the last partition isn't the most likely
  for(auto &cluster : cp.partitions()[cp.i_best()]) {
    if(args_->seed_unique_id() != "" && SeedMissing(cluster))
      continue;

    RecoEvent event;
    CalculateNaiveSeq(GetNaiveSeqNameToCalculate(cluster), &event);  // calculate the viterbi path from scratch to get the <event> set (should probably at some point start caching the events earlier)

    if(event.genes_["d"] == "") {  // shouldn't happen any more, but it is a check that could fail at some point
      cout << "WTF " << cluster << " x" << event.naive_seq_ << "x" << endl;
      assert(0);
    }
    StreamViterbiOutput(annotation_ofs, event, cachefo(cluster).seqs_, "");
  }
  annotation_ofs.close();
  printf("        annotation writing time (probably includes a bunch of new vtb calculations) %.1f\n", ((clock() - run_start) / (double)CLOCKS_PER_SEC));
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(Partition &partition, bool debug) {
  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  if(debug)
    cout << "LogProbOfPartition: " << endl;
  for(auto &key : partition) {
    double log_prob = GetLogProb(key);  // NOTE do *not* do any translation here -- we need the actual probability of the whole partition, to compare to the other partitions, so we need each and every sequence in each cluster (i.e. if you wanted to do translation, you'd have to coordinate the ignored sequences among the different partitions)
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
string Glomerator::FinalString(bool newline) {
    char buffer[2000];
    sprintf(buffer, "        calcd:   vtb %-4d  fwd %-4d  hfrac %-8d%s        merged:  hfrac %-4d lratio %-4d", n_vtb_calculated_, n_fwd_calculated_, n_hfrac_calculated_, newline ? "\n" : "", n_hfrac_merges_, n_lratio_merges_);
    return string(buffer);
}

// ----------------------------------------------------------------------------------------
string Glomerator::GetStatusStr(time_t current_time) {
  char timebuf[2000];
  strftime(timebuf, 2000, "%b %d %T", localtime(&current_time));
  stringstream ss;
  ss << "      " << timebuf;
  ss << "    " << setw(4) << current_partition_->size() << " clusters";
  ss << "    " << setw(9) << GetRss() << " / " << setw(1) << GetMemTot() << " kB = " << setw(6) << setprecision(3) << 100. * float(GetRss()) / GetMemTot() << " %";
  ss << "   " << FinalString(false);
  ss << "     " << ClusterSizeString(current_partition_).c_str();
  ss << endl;
  return ss.str();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WriteStatus() {

  // cout << CacheSizeString() << endl;
  time_t current_time;
  time(&current_time);
  if(difftime(current_time, last_status_write_time_) > 30) {  // write something every x seconds (if it crashes, partitiondriver prints the contents to stdout)
    string status_str(GetStatusStr(current_time));
    // cout << status_str;
    fprintf(progress_file_, "%s", status_str.c_str());
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
int Glomerator::CountMembers(string namestr, bool exclude_extra_seeds) {
  int n_colons = (int)count(namestr.begin(), namestr.end(), ':');
  int n_members(n_colons + 1);
  if(exclude_extra_seeds && args_->seed_unique_id() != "") {  // i'm adding this option long afterwards and remember nothing about this code, so there could be something already that could accomplish this
    vector<string> namevector(SplitString(namestr, ":"));
    int n_seeds = (int)count(namevector.begin(), namevector.end(), args_->seed_unique_id());
    if(n_seeds > 0)
      n_members -= n_seeds - 1;
  }
  return n_members;
}

// ----------------------------------------------------------------------------------------
// count the number of members in a cluster's colon-separated name string
unsigned Glomerator::LargestClusterSize(Partition &partition) {
  unsigned largest_cluster_size(0);
  for(auto &cluster : partition) {
    unsigned cluster_size = (unsigned)CountMembers(cluster);  // damnit, I wish I remembered why I made CountMembers() return a signed one
    if(largest_cluster_size == 0 or cluster_size > largest_cluster_size)
      largest_cluster_size = cluster_size;
  }
  return largest_cluster_size;
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
string Glomerator::JoinNameStrings(vector<Sequence*> &strlist, string delimiter) {
  string return_str;
  for(size_t is=0; is<strlist.size(); ++is) {
    if(is > 0)
      return_str += delimiter;
    return_str += strlist[is]->name();
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
string Glomerator::JoinSeqStrings(vector<Sequence*> &strlist, string delimiter) {
  string return_str;
  for(size_t is=0; is<strlist.size(); ++is) {
    if(is > 0)
      return_str += delimiter;
    return_str += strlist[is]->undigitized();
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
bool Glomerator::SeedMissing(string queries) {
  return cachefo(queries).seed_missing_;  // NOTE after refactoring the double loops, we probably don't really need to cache all these any more
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
    throw runtime_error("sequences different length in Glomerator::NaiveHfrac\n    " + to_string(seq_a.size()) + ": " + seq_a + "\n    " + to_string(seq_b.size()) + ": " + seq_b + "\n");
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
  if(seq_a.size() == 0 || seq_b.size() == 0 || failed_queries_.count(key_a) || failed_queries_.count(key_b))  // size 0 means there was no valid path (this *should* get picked up by looking in failed_queries_, but it's not because of translation, and I don't remember how this stuff works nearly well enough at this point to try and fix that)
    return hfrac;
  naive_hfracs_[joint_key] = CalculateHfrac(seq_a, seq_b);

  return naive_hfracs_[joint_key];
}

// ----------------------------------------------------------------------------------------
string Glomerator::ChooseSubsetOfNames(string queries, int n_max) {
  if(name_subsets_.count(queries))
    return name_subsets_[queries];

  // assert(seq_info_.count(queries) || tmp_cachefo_.count(queries));
  vector<string> namevector(SplitString(queries, ":"));
  set<string> nameset(namevector.begin(), namevector.end());
  if(nameset.size() < unsigned(n_max))  // with seed clustering you can get a lot of duplicate seqs (not just the seed seq), and I *think* this is maybe an ok way to avoid the runtime error below? But i'm adding this long afterwards so i'm not sure
    n_max = nameset.size();

  srand(hash<string>{}(queries));  // make sure we get the same subset each time we pass in the same queries (well, if there's different thresholds for naive_seqs annd logprobs they'll each get their own [very correlated] subset)

  // first decide which indices we'll choose
  set<int> ichosen;
  vector<int> ichosen_vec;  // don't really need both this and the set... but maybe it's faster
  set<string> chosen_strs;  // make sure we don't choose seed unique id more than once
  for(size_t iname=0; iname<unsigned(n_max); ++iname) {
    int ich(-1);
    int n_tries(0);
    while(ich < 0 || ichosen.count(ich) || chosen_strs.count(namevector[ich])) {
      ich = rand() % namevector.size();
      ++n_tries;
      if(n_tries > 1e6) {
	throw runtime_error("too many tries in Glomerator::ChooseSubsetOfNames() -- probably too many copies of some seqs during seed partitioning (choosing " + to_string(n_max) + " of " + to_string(namevector.size()) + " from " + queries);
      }
    }
    ichosen.insert(ich);
    chosen_strs.insert(namevector[ich]);
    ichosen_vec.push_back(ich);
  }

  // then sort 'em
  sort(ichosen_vec.begin(), ichosen_vec.end());

  Query &cacheref = cachefo(queries);

  // and finally make the new vectors
  vector<string> subqueryvec;
  for(auto &ich : ichosen_vec)
    subqueryvec.push_back(namevector[ich]);

  string subqueries(JoinStrings(subqueryvec));

  tmp_cachefo_[subqueries] = Query(subqueries,
				   GetSeqs(subqueries),
				   !InString(args_->seed_unique_id(), subqueries),
				   cacheref.only_genes_,
				   cacheref.kbounds_,
				   cacheref.mute_freq_,
				   cacheref.cdr3_length_);

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
  if(CountMembers(actual_queries, true) < 1.5 * args_->biggest_naive_seq_cluster_to_calculate())  // if <<actual_queries>> is small return all of 'em
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

  if(CountMembers(queries_to_calc, true) > n_max) {
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
  if(naive_seqs_.count(queries_to_calc) == 0) {
    string tmp_nseq = CalculateNaiveSeq(queries_to_calc);  // some compilers add <queries_to_calc> to <naive_seqs_> *before* calling CalculateNaiveSeq(), which causes that function's check to fail
    if(tmp_nseq.size() == 0) {
      cout << "  warning: zero length naive sequence for " << queries_to_calc << endl;
      return empty_string_;
    }
    naive_seqs_[queries_to_calc] = tmp_nseq;
  }

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
double Glomerator::GetLogProb(string queries) {  // NOTE this does *no* translation, so you better have done that already before you call it if you want it done
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

  Query full_qmerged = GetMergedQuery(key_a, key_b);
  pair<string, string> parents_to_calc = GetLogProbPairOfNamesToCalculate(joint_name, full_qmerged.parents_);
  string key_a_to_calc = parents_to_calc.first;
  string key_b_to_calc = parents_to_calc.second;
  Query qmerged_to_calc = GetMergedQuery(key_a_to_calc, key_b_to_calc);

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

  // if(seq_info_.count(queries) == 0 && tmp_cachefo_.count(queries) == 0)
  //   throw runtime_error("no info for " + queries);

  ++n_vtb_calculated_;

  DPHandler dph("viterbi", args_, gl_, hmms_);
  Query &cacheref = cachefo(queries);
  Result result = dph.Run(cacheref.seqs_, cacheref.kbounds_, cacheref.only_genes_, cacheref.mute_freq_);
  // if(FishyMultiSeqAnnotation(SplitString(queries).size(), result.best_event()))
  //   dph.HandleFishyAnnotations(result, cacheref.seqs_, cacheref.kbounds_, cacheref.only_genes_, cacheref.mute_freq_);
  if(result.no_path_) {
    AddFailedQuery(queries, "no_path");
    return "";
  }

  if(event != nullptr)
    *event = result.best_event();

  WriteStatus();
  return result.best_event().naive_seq_;
}

// ----------------------------------------------------------------------------------------
double Glomerator::CalculateLogProb(string queries) {  // NOTE can modify kbinfo_
  // NOTE do *not* call this from anywhere except GetLogProb()
  assert(log_probs_.count(queries) == 0);

  // if(seq_info_.count(queries) == 0 && tmp_cachefo_.count(queries) == 0)
  //   throw runtime_error("no info for " + queries);
  
  ++n_fwd_calculated_;

  DPHandler dph("forward", args_, gl_, hmms_);
  Query &cacheref = cachefo(queries);
  Result result = dph.Run(cacheref.seqs_, cacheref.kbounds_, cacheref.only_genes_, cacheref.mute_freq_);
  if(result.no_path_) {
    AddFailedQuery(queries, "no_path");
    return -INFINITY;
  }

  WriteStatus();
  return result.total_score();
}

// ----------------------------------------------------------------------------------------
void Glomerator::AddFailedQuery(string queries, string error_str) {
    errors_[queries] = errors_[queries] + ":" + error_str;
    failed_queries_.insert(queries);
}

// ----------------------------------------------------------------------------------------
Query &Glomerator::cachefo(string queries) {
  if(cachefo_.find(queries) != cachefo_.end())
    return cachefo_[queries];
  else if(tmp_cachefo_.find(queries) != tmp_cachefo_.end())
    return tmp_cachefo_[queries];
  else {  // if this is happening very frequently you've fucked up
    // throw runtime_error(queries + " not found in either cache\n");
    // cout << "hackadd to tmp cache " << queries << endl;
    set<string> only_gene_set;
    KBounds kbounds;
    double mute_freq_total(0.);
    size_t cdr3_length(0);
    vector<string> tmpvec(SplitString(queries));
    for(size_t is=0; is<tmpvec.size(); ++is) {
      Query &scache(single_seq_cachefo_[tmpvec[is]]);
      only_gene_set.insert(scache.only_genes_.begin(), scache.only_genes_.end());
      if(is==0)
	kbounds = scache.kbounds_;
      else
	kbounds = kbounds.LogicalOr(scache.kbounds_);
      mute_freq_total += scache.mute_freq_;
      if(is==0)
	cdr3_length = scache.cdr3_length_;
      if(cdr3_length != scache.cdr3_length_)
	throw runtime_error("cdr3 length mismatch " + to_string(cdr3_length) + " " + to_string(scache.cdr3_length_) + " for single query " + tmpvec[is] + " within " + queries);

      // cout << "    " << tmpvec[is] << "    " << kbounds.stringify() << "   " << scache.mute_freq_ << "   " << cdr3_length << "   ";
      // for(auto &g : only_gene_set)
      // 	cout << " " << g;
      // cout << endl;
    }
    // cout << "        final mute freq " << mute_freq_total / tmpvec.size() << endl;

    tmp_cachefo_[queries] = Query(queries,
				  GetSeqs(queries),
				  !InString(args_->seed_unique_id(), queries),
				  vector<string>(only_gene_set.begin(), only_gene_set.end()),
				  kbounds,
				  mute_freq_total / tmpvec.size(),
				  cdr3_length);
    return tmp_cachefo_[queries];
  }
}

// ----------------------------------------------------------------------------------------
bool Glomerator::SameLength(vector<Sequence*> &seqs, bool debug) {
  // are all the seqs in <seqs> the same length?
  bool same(true);
  int len(-1);
  for(auto &seq : seqs) {
    if(len < 0)
      len = (int)seq->size();
    if(len != (int)seq->size()) {
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
vector<Sequence*> Glomerator::GetSeqs(string query) {
  vector<string> queryvec(SplitString(query));
  vector<Sequence*> seqs(queryvec.size());
  for(size_t is=0; is<queryvec.size(); ++is) {
    if(single_seqs_.find(queryvec[is]) == single_seqs_.end())
      throw runtime_error("couldn't find query " + query + " in single seq vector");
    seqs[is] = &single_seqs_[queryvec[is]];
    if(seqs[is] == nullptr)
      throw runtime_error("null ptr to sequence for query " + query);
  }
  return seqs;
}

// ----------------------------------------------------------------------------------------
// when we're adding <query> to the permament cache in <cachefo_>, if it's been translated we also need it's subsets in <cachefo_>
void Glomerator::MoveSubsetsFromTmpCache(string query) {
  if(naive_seq_name_translations_.find(query) != naive_seq_name_translations_.end()) {
    string tquery(naive_seq_name_translations_[query]);
    // cout << "naive seq nt " << tquery << endl;
    CopyToPermanentCache(tquery, query);
  }

  if(logprob_name_translations_.find(query) != logprob_name_translations_.end()) {
    pair<string, string> tpair(logprob_name_translations_[query]);
    // cout << "logprob nt for: " << query << "    " << tpair.first << " " << tpair.second << endl;
    CopyToPermanentCache(tpair.first, query);
    CopyToPermanentCache(tpair.second, query);
  }

  if(logprob_asymetric_translations_.find(query) != logprob_asymetric_translations_.end()) {
    string tquery(logprob_asymetric_translations_[query]);
    // cout << "logprob asym t " << tquery << endl;
    CopyToPermanentCache(tquery, query);
  }
}

// ----------------------------------------------------------------------------------------
// Copy the entry for <translated_query> from <tmp_cachefo_> to <cachefo_>, unless it isn't there, in which case we reconstruct roughly what it should have been using <superquery> (the query for which <translated_query> is a translation).
// e.g. if <translated_query> is "is:hm" then <superquery> might be "az:fh:fi:is:fj:hm".
void Glomerator::CopyToPermanentCache(string translated_query, string superquery) {
  if(tmp_cachefo_.find(translated_query) != tmp_cachefo_.end()) {
    cachefo_[translated_query] = tmp_cachefo_[translated_query];
  } else {  // I think that if we don't have it even in the tmp cache, that we won't ever need the query info (I think it means to we already calculated everything for it) but it makes things more consistent and safer to make sure it's in the permanenet cache
    Query &supercache(cachefo(superquery));
    // cout << "scratchy! " << superquery << " --> " << translated_query << endl;
    cachefo_[translated_query] = Query(translated_query,
				       GetSeqs(translated_query),
				       !InString(args_->seed_unique_id(), translated_query),
				       supercache.only_genes_,
				       supercache.kbounds_,
				       supercache.mute_freq_,
				       supercache.cdr3_length_);
				       // supercache.parents_.first,
				       // supercache.parents_.second);
  }
}

// ----------------------------------------------------------------------------------------
Query &Glomerator::GetMergedQuery(string name_a, string name_b) {

  string joint_name = JoinNames(name_a, name_b);  // sorts name_a and name_b, but *doesn't* sort within them
  if(cachefo_.find(joint_name) != cachefo_.end())
    return cachefo_[joint_name];
  if(tmp_cachefo_.find(joint_name) != tmp_cachefo_.end())
    return tmp_cachefo_[joint_name];

  Query &ref_a = cachefo(name_a);
  Query &ref_b = cachefo(name_b);

  vector<string> joint_only_genes = ref_a.only_genes_;
  for(auto &g : ref_b.only_genes_) {
    if(find(joint_only_genes.begin(), joint_only_genes.end(), g) == joint_only_genes.end())
      joint_only_genes.push_back(g);
  }

  if(ref_a.cdr3_length_ != ref_b.cdr3_length_)
    throw runtime_error("cdr3 lengths different for " + name_a + " and " + name_b + " (" + to_string(ref_a.cdr3_length_) + " " + to_string(ref_b.cdr3_length_) + ")");

  // NOTE now that I'm adding the merged query to the cache info here, I can maybe get rid of the qmerged entirely UPDATE I have no idea if this is still relevant
  tmp_cachefo_[joint_name] = Query(joint_name,
				   GetSeqs(joint_name),
				   !InString(args_->seed_unique_id(), joint_name),
				   joint_only_genes,
				   ref_a.kbounds_.LogicalOr(ref_b.kbounds_),
				   (ref_a.seqs_.size()*ref_a.mute_freq_ + ref_b.seqs_.size()*ref_b.mute_freq_) / double(ref_a.seqs_.size() + ref_b.seqs_.size()),  // simple weighted average (doesn't account for different sequence lengths)
				   ref_a.cdr3_length_,
				   name_a,
				   name_b
				   );
  
  return tmp_cachefo_[joint_name];
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

      if(cachefo(key_a).cdr3_length_ != cachefo(key_b).cdr3_length_)
      	continue;

      double hfrac = NaiveHfrac(key_a, key_b);
      if(hfrac > args_->hamming_fraction_bound_hi())  // if naive hamming fraction too big, don't even consider merging the pair
	continue;

      if(args_->hamming_fraction_bound_lo() <= 0.0 || hfrac >= args_->hamming_fraction_bound_lo())
	continue;

      if(hfrac < min_hamming_fraction) {
	  min_hamming_fraction = hfrac;
	  min_hamming_merge = GetMergedQuery(key_a, key_b);
      }
    }
  }

  if(min_hamming_fraction != INFINITY) {  // (note that this is *plus* infinity, but in the lratio fcn it's -INFINITY)
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

      if(cachefo(key_a).cdr3_length_ != cachefo(key_b).cdr3_length_)
      	continue;

      double hfrac = NaiveHfrac(key_a, key_b);
      if(hfrac > args_->hamming_fraction_bound_hi())  // if naive hamming fraction too big, don't even consider merging the pair
	continue;

      double lratio = GetLogProbRatio(key_a, key_b);

      // don't merge if lratio is small (less than zero, more or less)
      if(!force_merge_ && LikelihoodRatioTooSmall(lratio, CountMembers(key_a) + CountMembers(key_b)))
	continue;

      if(lratio > max_lratio) {
	max_lratio = lratio;
	chosen_qmerge = GetMergedQuery(key_a, key_b);
      }
    }
  }

  if(max_lratio != -INFINITY) {  // if we found a merge that we liked (note that this is *minus* infinity, but in the hfrac fcn it's +INFINITY)
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
  pair<double, Query> qpair = FindHfracMerge(path);
  if(qpair.first == INFINITY)  // if there wasn't a good enough hfrac merge
    qpair = FindLRatioMerge(path);

  if(args_->max_cluster_size() > 0) {  // if we were told to stop if any clusters get too big
    for(auto &cluster : path->CurrentPartition()) {
      if(unsigned(CountMembers(cluster)) > args_->max_cluster_size()) {
	printf("    --max-cluster-size: stopping with a cluster of size %u (> %u)\n", unsigned(CountMembers(cluster)), args_->max_cluster_size());
	path->finished_ = true;
      }
    }
    if(path->finished_)
      return;
  }

  if(qpair.first == -INFINITY) {  // if there also wasn't a good lratio merge
    if(args_->n_final_clusters() == 0 && args_->min_largest_cluster_size() == 0) {  // default: stop when there's no good lratio merges, i.e. at (well, near) the maximum likelihood partition
      path->finished_ = true;
    } else if((args_->n_final_clusters() > 0 && path->CurrentPartition().size() > args_->n_final_clusters()) ||  // still have more clusters than we were asked for
	      (args_->min_largest_cluster_size() > 0 && LargestClusterSize(path->CurrentPartition()) < args_->min_largest_cluster_size())) {  // largest cluster is still too small
      if(force_merge_) {  // if we already set force merge on a previous iteration
	path->finished_ = true;
	if(args_->n_final_clusters() > 0)
	  printf("    couldn't merge beyond %lu clusters despite setting force merge\n", path->CurrentPartition().size());
	if(args_->min_largest_cluster_size() > 0)
	  printf("    couldn't merge past a biggest cluster of %u despite setting force merge\n", LargestClusterSize(path->CurrentPartition()));
      } else {
	force_merge_ = true;
	if(args_->n_final_clusters() > 0)
	  printf("    setting force merge (currently have %lu clusters, requested %u)\n", path->CurrentPartition().size(), args_->n_final_clusters());
	if(args_->min_largest_cluster_size() > 0)
	  printf("    setting force merge (current biggest cluster %u, requested %u)\n", LargestClusterSize(path->CurrentPartition()), args_->min_largest_cluster_size());
      }
    } else {  // we've gotten down to the requested number of clusters, so we can stop (shouldn't be possible to get here, since it'd require somehow missing setting the path to finished in the if clause below)
      path->finished_ = true;
    }
    return;
  }

  WriteStatus();
  Query chosen_qmerge = qpair.second;

  cachefo_[chosen_qmerge.name_] = chosen_qmerge;
  GetNaiveSeq(chosen_qmerge.name_, &chosen_qmerge.parents_);  // this *needs* to happen here so it has the parental information
  UpdateLogProbTranslationsForAsymetrics(chosen_qmerge);
  MoveSubsetsFromTmpCache(chosen_qmerge.name_);

  Partition new_partition(path->CurrentPartition());
  new_partition.erase(chosen_qmerge.parents_.first);
  new_partition.erase(chosen_qmerge.parents_.second);
  new_partition.insert(chosen_qmerge.name_);
  path->AddPartition(new_partition, -INFINITY, args_->n_partitions_to_write());
  current_partition_ = &path->CurrentPartition();

  if(args_->debug()) {
    printf("       merged   %s  %s\n", chosen_qmerge.parents_.first.c_str(), chosen_qmerge.parents_.second.c_str());
    cout << "          removing " << tmp_cachefo_.size() << " entries from tmp cache" << endl;
  }

  tmp_cachefo_.clear();  // NOTE I could simplify some other things if I only cleared the stuff from <tmp_cachefo_> that I thought I wouldn't later need.
  // naive_hfracs_.clear();  // so, this can reduce memory usage a *lot* for no cpu hit in (usually, I think) later steps, but in (usually, I think) earlier steps it can be prohibitively slower (I think, when early on you're doing a ton of hfrac merges)
  //  - bottom line, I think I'll need to buckle down and measure memory usage, decide how much I'm allowed to use, and only clear the caches that I need to
  //  - but also, it'd probably (maybe?) be ok to clear this cache after a logprob merge, since that would mean we're probably through with all the naive hfrac merges
  //  - or at least remove info for clusters we've merged out of existence

  if((args_->n_final_clusters() > 0 && path->CurrentPartition().size() <= args_->n_final_clusters()) ||
     (args_->min_largest_cluster_size() > 0 && LargestClusterSize(path->CurrentPartition()) >= args_->min_largest_cluster_size())) {  // largest cluster is still too small
    path->finished_ = true;
    if(args_->n_final_clusters() > 0)
      printf("    finished glomerating to %lu clusters (requested %u), force merge status %d\n", path->CurrentPartition().size(), args_->n_final_clusters(), force_merge_);
    if(args_->min_largest_cluster_size() > 0)
      printf("    finished glomerating to a biggest cluster of %u (requested %u), force merge status %d\n", LargestClusterSize(path->CurrentPartition()), args_->min_largest_cluster_size(), force_merge_);
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
