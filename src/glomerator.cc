#include "glomerator.h"

namespace ham {


// ----------------------------------------------------------------------------------------
string SeqStr(vector<Sequence> &seqs, string delimiter) {
  string seq_str;
  for(size_t iseq = 0; iseq < seqs.size(); ++iseq) {
    if(iseq > 0) seq_str += delimiter;
    seq_str += seqs[iseq].undigitized();
  }
  return seq_str;
}

// ----------------------------------------------------------------------------------------
string SeqNameStr(vector<Sequence> &seqs, string delimiter) {
  string name_str;
  for(size_t iseq = 0; iseq < seqs.size(); ++iseq) {
    if(iseq > 0) name_str += delimiter;
    name_str += seqs[iseq].name();
  }
  return name_str;
}

// ----------------------------------------------------------------------------------------
Glomerator::Glomerator(HMMHolder &hmms, GermLines &gl, vector<vector<Sequence> > &qry_seq_list, Args *args, Track *track) :
  track_(track),
  args_(args),
  dph_(args_, gl, hmms),
  i_initial_partition_(0),
  n_fwd_cached_(0),
  n_fwd_calculated_(0),
  n_vtb_cached_(0),
  n_vtb_calculated_(0)
{
  ReadCachedLogProbs();

  Partition tmp_partition;
  int last_ipath(0);
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
  cout << " cached vtb " << n_vtb_cached_ << "/" << (n_vtb_cached_ + n_vtb_calculated_)
       << "    fwd " << n_fwd_cached_ << "/" << (n_fwd_cached_ + n_fwd_calculated_) << endl;
  WriteCachedLogProbs();
}

// ----------------------------------------------------------------------------------------
void Glomerator::Cluster() {
  if(args_->debug()) cout << "   glomerating" << endl;

  assert((int)initial_partitions_.size() == 1);
  ClusterPath cp(initial_partitions_[0], LogProbOfPartition(initial_partitions_[0]), initial_logweights_[0]);
  do {
    Merge(&cp);
  } while(!cp.finished_);

  vector<ClusterPath> paths{cp};
  WritePartitions(paths);
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
  assert(headstrs[3].find("cyst_position") == 0);

  while(getline(ifs, line)) {
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    vector<string> column_list = SplitString(line, ",");
    assert(column_list.size() == 4 || column_list.size() == 5);  // if written by bcrham it'll have an error column, while if written by partitiondriver it won't (under normal circumstances we wouldn't see bcrham cachefiles right here, but it can be useful for testing)
    string unique_ids(column_list[0]);
    string logprob_str(column_list[1]);
    if(logprob_str.size() > 0)
      log_probs_[unique_ids] = stof(logprob_str);
    string naive_seq(column_list[2]);
    int cyst_position(atoi(column_list[3].c_str()));
    
    if(naive_seq.size() > 0)
      naive_seqs_[unique_ids] = Sequence(track_, unique_ids, naive_seq, cyst_position);
  }
  // cout << "      read " << log_probs_.size() << " cached results" << endl;
}

// ----------------------------------------------------------------------------------------
double Glomerator::LogProbOfPartition(Partition &partition) {
  // get log prob of entire partition given by the keys in <partinfo> using the individual log probs in <log_probs>
  double total_log_prob(0.0);
  for(auto &key : partition) {
    GetLogProb(key, seq_info_[key], kbinfo_[key], only_genes_[key], mute_freqs_[key]);  // immediately returns if we already have it
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
  log_prob_ofs << "unique_ids,score,naive-seq,cyst_position,errors" << endl;

  log_prob_ofs << setprecision(20);
  // first write everything for which we have log probs
  for(auto &kv : log_probs_) {
    string naive_seq, cyst_position_str;  // if we don't have the naive sequence, write empty strings for both
    if(naive_seqs_.count(kv.first)) {  // if we calculate the merged prob for two clusters, but don't end up merging them, then we will have the logprob but not the naive seq
      naive_seq = naive_seqs_[kv.first].undigitized();
      cyst_position_str = to_string(naive_seqs_[kv.first].cyst_position());
    }
    log_prob_ofs << kv.first << "," << kv.second << "," << naive_seq << "," << cyst_position_str << "," << errors_[kv.first] << endl;  // NOTE if there's no errors, it just prints the empty string, which is fine
  }

  // then write the queries for which we have naive seqs but not logprobs (if they were hamming-skipped )
  for(auto &kv : naive_seqs_) {
    if(log_probs_.count(kv.first))  // if it's in log_probs, we've already done it
      continue;
    assert(0);  // oh, wait, nevermind, I guess we *can't* actually have the naive seq if we don't have the log prob
    log_prob_ofs << kv.first << ",," << kv.second.undigitized() << "," << kv.second.cyst_position() << "," << errors_[kv.first] << endl;  // NOTE if there's no errors, it just prints the empty string, which is fine
  }

  log_prob_ofs.close();
}

// ----------------------------------------------------------------------------------------
void Glomerator::WritePartitions(vector<ClusterPath> &paths) {
  ofs_.open(args_->outfile());
  ofs_ << setprecision(20);
  ofs_ << "path_index,initial_path_index,partition,score,logweight" << endl;
  int ipath(0);
  for(auto &cp : paths) {
    for(unsigned ipart=0; ipart<cp.partitions().size(); ++ipart) {
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
double Glomerator::NaiveHammingFraction(string key_a, string key_b) {
  GetNaiveSeq(key_a);
  GetNaiveSeq(key_b);
  vector<Sequence> nseqs{naive_seqs_[key_a], naive_seqs_[key_b]};
  if(args_->truncate_seqs()) {
    vector<KBounds> kbvector{kbinfo_[key_a], kbinfo_[key_b]};  // don't need kbounds here, but we need to pass in something
    TruncateSeqs(nseqs, kbvector);
    printf("  truncate in NaiveHammingDistance: %d, %d --> %d, %d        %s %s\n",
	   int(naive_seqs_[key_a].size()), int(naive_seqs_[key_b].size()),
	   int(nseqs[0].size()), int(nseqs[1].size()),
	   key_a.c_str(), key_b.c_str());
  }
  return HammingDistance(nseqs[0], nseqs[1]) / double(nseqs[0].size());  // hamming distance fcn will fail if the seqs aren't the same length
}

// ----------------------------------------------------------------------------------------
// truncate sequences in <seqs> to the same length (on both ends of the conserved cysteine), and correspondingly modify <kbvector>
void Glomerator::TruncateSeqs(vector<Sequence> &seqs, vector<KBounds> &kbvector) {
  assert(seqs.size() == kbvector.size());  // one kbound for each sequence

  // first find min length to left and right of the cysteine position
  int min_left(-1), min_right(-1);
  for(size_t is=0; is<seqs.size(); ++is) {
    Sequence *seq(&seqs[is]);
    int cpos(seq->cyst_position());
    if(cpos < 0 || cpos >= (int)seq->size())
      throw runtime_error("cpos " + to_string(cpos) + " invalid for " + seq->name() + " (" + seq->undigitized() + ")");
    int dleft = cpos;  // NOTE <dright> includes <cpos>, i.e. dleft + dright = len(seq)
    int dright = seq->size() - cpos;
    printf("        %d %d  (%d, %d - %d)   %s\n", dleft, dright, (int)cpos, (int)seq->size(), (int)cpos, seq->name().c_str());
    if(min_left == -1 || dleft < min_left) {
      min_left = dleft;
    }
    if(min_right == -1 || dright < min_right) {
      min_right = dright;
    }
  }
  assert(min_left >= 0 && min_right >= 0);
  // printf("  min left %d right %d\n", min_left, min_right);

  // then truncate all the sequences to these lengths
  for(size_t is=0; is<seqs.size(); ++is) {
    Sequence *seq(&seqs[is]);
    int cpos(seq->cyst_position());
    int istart(cpos - min_left);
    int istop(cpos + min_right);
    int chopleft(istart);
    int chopright(seq->size() - istop);

    // string truncated_str(seq->undigitized().substr(istart, istop - istart));
    // Sequence trunc_seq(seq->track(), seq->name(), truncated_str, cpos - istart);
    // seqs[is] = trunc_seq;  // replace the old sequence with the truncated one
    seqs[is] = Sequence(*seq, istart, istop - istart);  // replace the old sequence with the truncated one (this sets cpos in the new sequence to cpos in the old sequence minus istart)
    assert(chopleft < (int)kbvector[is].vmin);  // kinda nonsensical if we start chopping off the entire v
    kbvector[is].vmin -= chopleft;
    kbvector[is].vmax -= chopleft;
    
    // printf("%s", seq->name());
    printf("      chop %d %d   %s\n", chopleft, chopright, seq->name().c_str());
    // printf("  before", self.sw_info[name]['k_v']['min'], self.sw_info[name]['k_v']['max'], self.sw_info[name]['v_5p_del'], self.sw_info[name]['j_3p_del'], self.sw_info[name]['seq']
    // self.sw_info[name]['seq'] = seq[istart : istop]
    // self.sw_info[name]['k_v']['min'] -= chopleft
    // self.sw_info[name]['k_v']['max'] -= chopleft
    // self.sw_info[name]['v_5p_del'] += chopleft
    // self.sw_info[name]['j_3p_del'] += chopright
    // print '   after', self.sw_info[name]['k_v']['min'], self.sw_info[name]['k_v']['max'], self.sw_info[name]['v_5p_del'], self.sw_info[name]['j_3p_del'], self.sw_info[name]['seq']
  }
}

// ----------------------------------------------------------------------------------------
// add log prob for <key>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetNaiveSeq(string key) {
  if(naive_seqs_.count(key)) {  // already did it (ok to cache naive seqs even when we're truncating, since each sequence, when part of a given group of sequence, always has the same length [it's different for forward because each key is compared in the likelihood ratio to many other keys, and each time its sequences can potentially have a different length])
    ++n_vtb_cached_;
    return;
  }
  ++n_vtb_calculated_;

  Result result(kbinfo_[key]);
  bool stop(false);
  do {
    result = dph_.Run("viterbi", seq_info_[key], kbinfo_[key], only_genes_[key], mute_freqs_[key]);
    kbinfo_[key] = result.better_kbounds();
    stop = !result.boundary_error() || result.could_not_expand();  // stop if the max is not on the boundary, or if the boundary's at zero or the sequence length
    if(args_->debug() && !stop)
      cout << "      expand and run again" << endl;  // note that subsequent runs are much faster than the first one because of chunk caching
  } while(!stop);

  if(result.events_.size() < 1)
    throw runtime_error("ERROR no events for " + key + "\n");

  naive_seqs_[key] = Sequence(track_, key, result.events_[0].naive_seq_, result.events_[0].cyst_position_);
  if(result.boundary_error())
    errors_[key] = errors_[key] + ":boundary";
}

// ----------------------------------------------------------------------------------------
// add log prob for <name>/<seqs> to <log_probs_> (if it isn't already there)
void Glomerator::GetLogProb(string name, vector<Sequence> &seqs, KBounds &kbounds, vector<string> &only_genes, double mean_mute_freq) {
  // NOTE that when this improves the kbounds, that info doesn't get propagated to <kbinfo_>
  if(log_probs_.count(name) && !args_->truncate_seqs()) {  // already did it (*no* caching if we're truncating sequences, since I don't want to make sequence truncation parameters part of <key> here [it's ok to cache viterbi, though -- see note in function above])
    ++n_fwd_cached_;
    return;
  }
  ++n_fwd_calculated_;
    
  Result result(kbounds);
  bool stop(false);
  do {
    result = dph_.Run("forward", seqs, kbounds, only_genes, mean_mute_freq);  // NOTE <only_genes> isn't necessarily <only_genes_[name]>, since for the denominator calculation we take the OR
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
vector<Sequence> Glomerator::MergeSeqVectors(string name_a, string name_b) {
  // NOTE does *not* truncate anything

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
    if(all_names.count(name))
      throw runtime_error("tried to add sequence with name " + name + " twice in Glomerator::MergeSeqVectors()");
    else
      all_names.insert(name);
  }

  return merged_seqs;
}

// ----------------------------------------------------------------------------------------
bool Glomerator::SameLength(vector<Sequence> &seqs) {
  // are all the seqs in <seqs> the same length?
  int len(-1);
  for(auto &seq : seqs) {
    if(len < 0)
      len = (int)seq.size();
    if(len != (int)seq.size())
      return false;
  }
  return true;
}  

// ----------------------------------------------------------------------------------------
Query Glomerator::GetMergedQuery(string name_a, string name_b) {
  // NOTE truncates sequences!

  Query qmerged;
  qmerged.name_ = JoinNames(name_a, name_b);
  qmerged.seqs_ = MergeSeqVectors(name_a, name_b);  // doesn't truncate anything

  // truncation
  assert(SameLength(seq_info_[name_a]));  // all the seqs for name_a should already be the same length
  assert(SameLength(seq_info_[name_b]));  // ...same for name_b
  vector<KBounds> kbvector;  // make vector of kbounds, with first chunk corresponding to <name_a>, and the second to <name_b>. This shenaniganery is necessary so we can take the OR of the two kbounds *after* truncation
  for(size_t is=0; is<qmerged.seqs_.size(); ++is) {
    if(is < seq_info_[name_a].size())  // the first part of the vector is for sequences from name_a
      kbvector.push_back(kbinfo_[name_a]);
    else
      kbvector.push_back(kbinfo_[name_b]);  // ... and second part is for those from name_b
  }

  if(args_->truncate_seqs()) {
    TruncateSeqs(qmerged.seqs_, kbvector);
    printf("  truncate in GetMergedQuery: %d, %d --> %d, %d (%d-%d, %d-%d --> %d-%d, %d-%d)      %s %s\n",
	   int(seq_info_[name_a][0].size()), int(seq_info_[name_b][0].size()), int(qmerged.seqs_[0].size()), int(qmerged.seqs_.back().size()),
	   int(kbinfo_[name_a].vmin), int(kbinfo_[name_a].vmax),
	   int(kbinfo_[name_b].vmin), int(kbinfo_[name_b].vmax),
	   int(kbvector[0].vmin), int(kbvector[0].vmax),
	   int(kbvector.back().vmin), int(kbvector.back().vmax),
	   name_a.c_str(), name_b.c_str());
  }
  qmerged.kbounds_ = kbvector[0].LogicalOr(kbvector.back());  // OR of the kbounds for name_a and name_b, after they've been properly truncated

  qmerged.only_genes_ = only_genes_[name_a];
  for(auto &g : only_genes_[name_b])  // NOTE this will add duplicates (that's no big deal, though) OPTIMIZATION
    qmerged.only_genes_.push_back(g);
  // TODO mute freqs aren't really right any more if we truncated things above
  qmerged.mean_mute_freq_ = (seq_info_[name_a].size()*mute_freqs_[name_a] + seq_info_[name_b].size()*mute_freqs_[name_b]) / float(qmerged.seqs_.size());  // simple weighted average (doesn't account for different sequence lengths)
  qmerged.parents_ = pair<string, string>(name_a, name_b);
  return qmerged;
}

// ----------------------------------------------------------------------------------------
Query *Glomerator::ChooseRandomMerge(vector<pair<double, Query> > &potential_merges, smc::rng *rgen) {
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
      return &potential_merges[im].second;
    }
  }

  throw runtime_error("fell through in Glomerator::ChooseRandomMerge");
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
  if(path->finished_)  // already finished this <path>, but we're still iterating 'cause some of the other paths aren't finished
    return;

  double max_log_prob(-INFINITY);
  int imax(-1);

  vector<pair<double, Query> > potential_merges;

  int n_total_pairs(0), n_skipped_hamming(0);

  set<string> already_done;  // keeps track of which  a-b pairs we've already done
  for(auto &key_a : path->CurrentPartition()) {  // note that c++ maps are ordered
    for(auto &key_b : path->CurrentPartition()) {  // also, note that CurrentPartition() returns a reference
      if(key_a == key_b) continue;
      if(already_done.count(JoinNames(key_a, key_b)))  // already did this pair (this is kind of a weird way of looping only over unique pairs... oh, well. Wish I had python itertools here)
  	continue;
      else
  	already_done.insert(JoinNames(key_a, key_b));

      ++n_total_pairs;

      // NOTE it might help to also cache hamming fractions
      if(NaiveHammingFraction(key_a, key_b) > args_->hamming_fraction_cutoff()) {  // truncates naive sequences before passing to actual hamming distance function
	++n_skipped_hamming;
	continue;
      }

      Query qmerged(GetMergedQuery(key_a, key_b));  // truncates <key_a> and <key_b> sequences to the same length

      // NOTE need to get the a_seqs and b_seqs from <qmerged> (instead of <seq_info_>) so they're properly truncated
      vector<Sequence> a_seqs, b_seqs;
      for(size_t is=0; is<qmerged.seqs_.size(); ++is) {
	if(is < seq_info_[key_a].size())  // the first part of the vector is for sequences from key_a
	  a_seqs.push_back(qmerged.seqs_[is]);
	else
	  b_seqs.push_back(qmerged.seqs_[is]);  // ... and second part is for those from key_b
      }
	
      // NOTE the error from using the single kbounds rather than the OR seems to be around a part in a thousand or less
      // NOTE also that the _a and _b results will be cached (unless we're truncating), but with their *individual* only_gene sets (rather than the OR)... but this seems to be ok.

      GetLogProb(key_a, a_seqs, qmerged.kbounds_, qmerged.only_genes_, mute_freqs_[key_a]);  // it is ABSOLUTELY CRUCIAL that caching is turned off here if we're truncating, since a_seqs is different for each choice of b_seqs
      GetLogProb(key_b, b_seqs, qmerged.kbounds_, qmerged.only_genes_, mute_freqs_[key_b]);
      GetLogProb(qmerged.name_, qmerged.seqs_, qmerged.kbounds_, qmerged.only_genes_, qmerged.mean_mute_freq_);

      // it's a likelihood ratio, not a bayes factor
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
	imax = potential_merges.size() - 1;
      }
    }
  }

  // if <seq_info_> only has one cluster, if hamming is too large between all remaining clusters, or if remaining bayes factors are -INFINITY
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
    chosen_logprob = log_probs_[qmerged->name_];
  }
  
  // add query info for the chosen merge, unless we already did it (e.g. if another particle already did this merge)
  if(seq_info_.count(qmerged->name_) == 0) {  // if we're doing smc, this will happen once for each particle that wants to merge these two. NOTE you get a very, very strange seg fault at the Sequences::Union line above, I *think* inside the implicit copy constructor. Yes, I should just define my own copy constructor, but I can't work out how to navigate through the const jungle a.t.m.
    seq_info_[qmerged->name_] = qmerged->seqs_;
    kbinfo_[qmerged->name_] = qmerged->kbounds_;
    mute_freqs_[qmerged->name_] = qmerged->mean_mute_freq_;
    only_genes_[qmerged->name_] = qmerged->only_genes_;

    GetNaiveSeq(qmerged->name_);  // TODO hm, wait, actually maybe I can leave naive seq caching on even if we're truncating sequences?
  }

  Partition new_partition(path->CurrentPartition());  // note: CurrentPartition() returns a reference
  new_partition.erase(qmerged->parents_.first);
  new_partition.erase(qmerged->parents_.second);
  new_partition.insert(qmerged->name_);
  path->AddPartition(new_partition, LogProbOfPartition(new_partition));

  if(args_->debug()) {
    // cout << "    path " << path->initial_path_index_ << endl;
    printf("          hamming skipped %d / %d\n", n_skipped_hamming, n_total_pairs);
    // printf("       merged %-8.2f %s and %s\n", max_log_prob, max_pair.first.c_str(), max_pair.second.c_str());
    printf("       merged %-8.2f %s and %s\n", chosen_logprob, qmerged->parents_.first.c_str(), qmerged->parents_.second.c_str());
    string extrastr("current (logweight " + to_string(path->CurrentLogWeight()) + ")");
    PrintPartition(new_partition, extrastr);
  }
}

}
