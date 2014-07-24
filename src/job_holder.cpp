#include "job_holder.h"

// ----------------------------------------------------------------------------------------
model *HMMHolder::Get(string gene) {
  if (hmms_.find(gene) == hmms_.end()) {  // if we don't already have it, read it from disk
    hmms_[gene] = new model;
    hmms_[gene]->import(hmm_dir_ + "/" + gl_.SanitizeName(gene) + ".hmm");
    hmms_[gene]->n_seqs_per_track_ = n_seqs_per_track_;
  }
  return hmms_[gene];
}  

// ----------------------------------------------------------------------------------------
HMMHolder::~HMMHolder() {
  for (auto &entry: hmms_)
    delete entry.second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JobHolder::JobHolder(size_t n_seqs_per_track, string algorithm, sequences *seqs, HMMHolder *hmms, string only_gene_str):
  n_seqs_per_track_(n_seqs_per_track),
  algorithm_(algorithm),
  seqs_(seqs),
  hmms_(hmms),
  debug_(0),
  n_best_events_(5)
{
  if (only_gene_str.size() > 0) {
    for (auto &region: gl_.regions_)
      only_genes_[region] = set<string>();
    while (true) {
      size_t i_next_colon(only_gene_str.find(":"));
      string gene = only_gene_str.substr(0,i_next_colon);  // get the next gene name
      only_genes_[gl_.GetRegion(gene)].insert(gene);
      only_gene_str = only_gene_str.substr(i_next_colon+1);  // then excise it from only_gene_str
      if (i_next_colon == string::npos)
	break;
    }
    for (auto &region: gl_.regions_)
      assert(only_genes_[region].size() > 0);
  }
}

// ----------------------------------------------------------------------------------------
JobHolder::~JobHolder() {
  for (auto &gene_map: trellisi_) {
    string gene(gene_map.first);
    for (auto &query_str_map: gene_map.second) {
      delete query_str_map.second;
      delete paths_[gene][query_str_map.first];
    }
  }
}

// ----------------------------------------------------------------------------------------
sequences JobHolder::GetSubSeqs(KSet kset, string region) {
  // get subsequences for one region
  size_t k_v(kset.first),k_d(kset.second);
  if (region=="v")
    return seqs_->getSubSequences(0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  else if (region=="d")
    return seqs_->getSubSequences(k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  else if (region=="j")
    return seqs_->getSubSequences(k_v + k_d, seqs_->size() - k_v - k_d);  // j region runs from k_v + k_d to end
  else
    assert(0);
}

// ----------------------------------------------------------------------------------------
map<string,sequences> JobHolder::GetSubSeqs(KSet kset) {
  // get subsequences for all regions
  map<string,sequences> subseqs;
  for (auto &region: gl_.regions_)
    subseqs[region] = GetSubSeqs(kset, region);
  return subseqs;
}

// ----------------------------------------------------------------------------------------
void JobHolder::Run(size_t k_v_start, size_t n_k_v, size_t k_d_start, size_t n_k_d) {
  // loop over k_v k_d space
  double best_score(-INFINITY),total_score_(-INFINITY);
  KSet best_kset(make_pair(0,0));
  if (debug_) cout << "    k_v: " << k_v_start << "-" << k_v_start+n_k_v-1 << "  k_d: " << k_d_start << "-" << k_d_start+n_k_d-1 << endl;
  for (size_t k_v=k_v_start; k_v<k_v_start+n_k_v; ++k_v) {
    for (size_t k_d=k_d_start; k_d<k_d_start+n_k_d; ++k_d) {
      KSet kset(k_v,k_d);
      RunKSet(kset);
      total_score_ = AddInLogSpace(total_scores_[kset], total_score_);  // sum up the probabilities for each kset, log P_tot = log \sum_i P_k_i
      if (debug_==2 and algorithm_=="forward") printf("            %9.2f (%.1e)  tot: %7.2f\n", total_scores_[kset], exp(total_scores_[kset]), total_score_);
      if (best_scores_[kset] > best_score) {
      	best_score = best_scores_[kset];
      	best_kset = kset;
      }
    }
  }

  // return if no valid path
  if (best_kset.first == 0) {
    if (debug_) cout << "  nothing sensical" << endl;
    return;
  }

  // sort vector of events by score, stream info to stderr, and print the top n_best_events_
  if (algorithm_ == "viterbi") {
    sort(events_.begin(), events_.end());
    reverse(events_.begin(), events_.end());
    if (debug_==2) {
      for (size_t ievt=0; ievt<n_best_events_; ++ievt) {
	events_[ievt].Print(gl_, 0, 0, false, false, "          ");  // man, I wish I had keyword args
	if (n_seqs_per_track_ == 2)
	  events_[ievt].Print(gl_, 0, 0, true, true, "          ");
      }
    }
  }
  StreamOutput(total_score_);  // NOTE this must happen after sorting in viterbi

  // warn if we're on the k_v k_d space boundary
  if (best_kset.first == k_v_start ||
      best_kset.first == k_v_start+n_k_v-1 ||
      best_kset.second == k_d_start ||
      best_kset.second == k_d_start+n_k_d-1) {
    cout << "    WARNING maximum at boundary"
	 << "  k_v: " << best_kset.first << "(" << k_v_start << "-" << k_v_start+n_k_v-1 << ")"
	 << "  k_d: " << best_kset.second << "(" << k_d_start << "-" << k_d_start+n_k_d-1 << ")" << endl;
  }

  // print debug info
  if (debug_) {
    if (algorithm_=="viterbi")
      cout << "    best kset: " << setw(4) << best_kset.first << setw(4) << best_kset.second << setw(12) << best_score << endl;
    else
      cout << "    sum over ksets: " << total_score_ << endl;
  }
}

// ----------------------------------------------------------------------------------------
void JobHolder::FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score) {
  // if we already did this gene for this query sequence. NOTE if we have one gene for this <query_strs>, we're always gonna have the rest of 'em too
  bool already_cached = trellisi_.find(gene) != trellisi_.end() && trellisi_[gene].find(query_strs) != trellisi_[gene].end();
  if (already_cached) {
    if (algorithm_=="viterbi") {
      *score = paths_[gene][query_strs] ? paths_[gene][query_strs]->getScore() : -INFINITY;  // it's set to nullptr if no valid path through hmm
      if (debug_ == 2) PrintPath(query_strs, gene, *score, "(cached)");
    } else if (algorithm_=="forward") {
      *score = trellisi_[gene][query_strs]->getForwardProbability();
    } else {
      assert(0);
    }
  } else {
    // initialize trellis and path
    if (trellisi_.find(gene) == trellisi_.end()) {
      trellisi_[gene] = map<StrPair,trellis*>();
      paths_[gene] = map<StrPair,traceback_path*>();
    }
    trellisi_[gene][query_strs] = new trellis(hmms_->Get(gene), query_seqs);

    if (algorithm_=="viterbi") {
      trellisi_[gene][query_strs]->viterbi();
      if (trellisi_[gene][query_strs]->ending_viterbi_score == -INFINITY) {  // no valid path through hmm. TODO fix this in a more general way
	*score = -INFINITY;
	paths_[gene][query_strs] = nullptr;
	if (debug_ == 2) cout << "                    " << gene << " " << *score << endl;
      } else {
	paths_[gene][query_strs] = new traceback_path(hmms_->Get(gene));
	trellisi_[gene][query_strs]->traceback(*paths_[gene][query_strs]);
	*score = paths_[gene][query_strs]->getScore();
	if (debug_ == 2) PrintPath(query_strs, gene, *score);
      }
      assert(fabs(*score) > 1e-200);
      assert(*score == -INFINITY || paths_[gene][query_strs]->size() > 0);
    } else if (algorithm_=="forward") {
      trellisi_[gene][query_strs]->forward();
      paths_[gene][query_strs] = nullptr;  // avoids violating the assumption that paths_ and trellisi_ have the same entries
      *score = trellisi_[gene][query_strs]->getForwardProbability();
    } else {
      assert(0);
    }
  }
}

// ----------------------------------------------------------------------------------------
void JobHolder::PrintPath(StrPair query_strs, string gene, double score, string extra_str) {  // NOTE query_str is seq1xseq2 for pair hmm
  if (score == -INFINITY) {
    cout << "                    " << gene << " " << score << endl;
    return;
  }
  vector<string> path_labels;
  paths_[gene][query_strs]->label(path_labels);
  if (path_labels.size() == 0) {
    if (debug_) cout << "                     " << gene << " has no valid path" << endl;
    return; // TODO fix this upstream. well, it isn't *broken*, but, you know, whatever
  }
  assert(path_labels.size() > 0);  // this will happen if the ending viterbi prob is 0, i.e. if there's no valid path through the hmm (probably the sequence or hmm lengths are screwed up)
  assert(path_labels.size() == query_strs.first.size());
  size_t insert_length = GetInsertLength(path_labels);
  size_t left_erosion_length = GetErosionLength("left", path_labels, gene);
  size_t right_erosion_length = GetErosionLength("right", path_labels, gene);

  string germline(gl_.seqs_[gene]);
  string modified_seq = germline.substr(left_erosion_length, germline.size() - right_erosion_length - left_erosion_length);
  for (size_t i=0; i<insert_length; ++i)
    modified_seq = "i" + modified_seq;
  assert(modified_seq.size() == query_strs.first.size());
  assert(germline.size() + insert_length - left_erosion_length - right_erosion_length == query_strs.first.size());
  TermColors tc;
  cout
    << "                    "
    << (left_erosion_length>0 ? ".." : "  ") << tc.ColorMutants("red", query_strs.first, modified_seq, query_strs.second) << (right_erosion_length>0 ? ".." : "  ")
    << "  " << extra_str
    << setw(12) << paths_[gene][query_strs]->getScore()
    << setw(25) << gene
    << endl;
}

// ----------------------------------------------------------------------------------------
void JobHolder::PushBackRecoEvent(KSet kset, map<string,string> &best_genes, double score) {
  RecoEvent event(FillRecoEvent(kset, best_genes, score));
  events_.push_back(event);
}

// ----------------------------------------------------------------------------------------
RecoEvent JobHolder::FillRecoEvent(KSet kset, map<string,string> &best_genes, double score) {
  RecoEvent event;
  StrPair seqs;
  for (auto &region: gl_.regions_) {
    StrPair query_strs(GetQueryStrs(kset, region));
    assert(best_genes.find(region) != best_genes.end());
    string gene(best_genes[region]);
    vector<string> path_labels;
    paths_[gene][query_strs]->label(path_labels);
    if (path_labels.size() == 0) {
      if (debug_) cout << "                     " << gene << " has no valid path" << endl;
      event.SetScore(-INFINITY);
      return event; // TODO fix this upstream. well, it isn't *broken*, but, you know, whatever
    }
    assert(path_labels.size() > 0);
    assert(path_labels.size() == query_strs.first.size());
    event.SetGene(region, gene);

    // set right-hand deletions
    event.SetDeletion(region + "_3p", GetErosionLength("right", path_labels, gene));
    // and left-hand deletions
    event.SetDeletion(region + "_5p", GetErosionLength("left", path_labels, gene));
    if (region=="d")
      event.SetInsertion("vd", query_strs.first.substr(0, GetInsertLength(path_labels)));
    if (region=="j")
      event.SetInsertion("dj", query_strs.first.substr(0, GetInsertLength(path_labels)));
    seqs.first += query_strs.first;
    seqs.second += query_strs.second;
  }

  // if (seqs_->size() > 0)  // erm, this segfaults a.t.m. I must be forgetting something somewhere else
  //   assert((*seqs_)[0].name_ == (*seqs_)[1].name_);  // er, another somewhat neurotic consistency check

  event.SetSeq((*seqs_)[0].name_, seqs.first);
  if (n_seqs_per_track_ == 2) {
    // assert((*seqs_)[0].name_ == (*seqs_)[1].name_);  don't recall at this point precisely why it was that I wanted this here
    event.SetSecondSeq((*seqs_)[1].name_, seqs.second);
  }
  event.SetScore(score);
}

// ----------------------------------------------------------------------------------------
void JobHolder::StreamOutput(double score) {
  if (algorithm_ == "viterbi") {
    cerr << "unique_id,v_gene,d_gene,j_gene,vd_insertion,dj_insertion,v_3p_del,d_5p_del,d_3p_del,j_5p_del,score,seq,second_seq" << endl;
    for (size_t ievt=0; ievt<n_best_events_; ++ievt) {
      RecoEvent *event = &events_[ievt];
      cerr
	<< event->seq_name_
	<< "," << event->genes_["v"]
	<< "," << event->genes_["d"]
	<< "," << event->genes_["j"]
	<< "," << event->insertions_["vd"]
	<< "," << event->insertions_["dj"]
	<< "," << event->deletions_["v_3p"]
	<< "," << event->deletions_["d_5p"]
	<< "," << event->deletions_["d_3p"]
	<< "," << event->deletions_["j_5p"]
	<< "," << event->score_
	<< "," << event->seq_
	<< ",";
      if (n_seqs_per_track_ == 2)
	cerr << event->second_seq_;
      cerr << endl;
    }
  } else {
    // cerr << total_score_ << endl;  // TODO seriously, wtf? this doesn't work, i.e. total_score_ is set before and after the call the StreamOutput, but is magically unset inside here
    cerr << score << endl;
  }
}
// ----------------------------------------------------------------------------------------
StrPair JobHolder::GetQueryStrs(KSet kset, string region) {
  sequences query_seqs(GetSubSeqs(kset,region));
  StrPair query_strs;
  query_strs.first = query_seqs[0].undigitize();
  if (query_seqs.size() == 2) {  // the sequences class should already ensure that both seqs are the same length
    assert(n_seqs_per_track_ == 2);
    query_strs.second = query_seqs[1].undigitize();
  }
  return query_strs;  
}
    
        
// ----------------------------------------------------------------------------------------
// add two numbers, treating -INFINITY as zero, i.e. calculates log a*b = log a + log b, i.e. a *and* b
double JobHolder::AddWithMinusInfinities(double first, double second) {
  if (first == -INFINITY || second == -INFINITY)
    return -INFINITY;
  else
    return first + second;
}

// ----------------------------------------------------------------------------------------
void JobHolder::RunKSet(KSet kset) {
  map<string,sequences> subseqs(GetSubSeqs(kset));
  best_scores_[kset] = -INFINITY;
  total_scores_[kset] = -INFINITY;  // total log prob of this kset, i.e. log(P_v * P_d * P_j), where e.g. P_v = \sum_i P(v_i k_v)
  map<string,double> regional_best_scores;  // the best score for each region
  map<string,double> regional_total_scores;  // the total score for each region, i.e. log P_v
  if (debug_ == 2)
    cout << "            " << kset.first << " " << kset.second << " -------------------" << endl;
  for (auto &region : gl_.regions_) {
    StrPair query_strs(GetQueryStrs(kset, region));

    TermColors tc;
    if (debug_ == 2) {
      if (algorithm_=="viterbi") {
	cout << "              " << region << " query " << tc.ColorMutants("purple", query_strs.second, query_strs.first) << endl;
	if (n_seqs_per_track_==2)
	  cout << "              " << region << " query " << tc.ColorMutants("purple", query_strs.first, query_strs.second) << endl;
      } else {
	cout << "              " << region << endl;
      }
    }

    regional_best_scores[region] = -INFINITY;
    regional_total_scores[region] = -INFINITY;
    size_t igene(0),n_short_v(0),n_long_erosions(0);
    for (auto &gene : gl_.names_[region]) {
      if (only_genes_[region].size()>0 and only_genes_[region].find(gene)==only_genes_[region].end())
	continue;
      igene++;
      
      if (region=="v" && query_strs.first.size()>gl_.seqs_[gene].size()) {  // query sequence too long for this v version to make any sense (ds and js have inserts so this doesn't affect them)
	if (debug_ == 2) cout << "                     " << gene << " too short" << endl;
	n_short_v++;
	continue;
      }
      if (query_strs.first.size() < gl_.seqs_[gene].size() - 10)  // entry into the left side of the v hmm is a little hacky, and is governed by a gaussian with width 5 (hmmwriter::fuzz_around_v_left_edge)
	n_long_erosions++;
	
      double gene_score(-INFINITY);
      FillTrellis(&subseqs[region], query_strs, gene, &gene_score);
      double gene_choice_score = log(hmms_->Get(gene)->overall_gene_prob_);  // TODO think through this again, and make sure it's correct for forward score, as well. I mean, I *think* it's right, but I could stand to go over it again
      gene_score = AddWithMinusInfinities(gene_score, gene_choice_score);
      regional_total_scores[region] = AddInLogSpace(gene_score, regional_total_scores[region]);  // (log a, log b) --> log a+b, i.e. here we are summing probabilities in log space, i.e. a *or* b
      if (debug_ == 2 && algorithm_=="forward") printf("                %9.2f (%.1e)  tot: %7.2f  %s\n", gene_score, exp(gene_score), regional_total_scores[region], tc.ColorGene(gene).c_str());
      if (gene_score > regional_best_scores[region]) {
	regional_best_scores[region] = gene_score;
	best_genes_[kset][region] = gene;
      }
    }

    // return if we didn't find a valid path for this region
    if (best_genes_[kset].find(region) == best_genes_[kset].end()) {
      assert(n_long_erosions == 0);  // adding assert because if it happens, that means my v_right_length was screwed up
      if (debug_ == 2) {
	cout << "                  found no gene for " << region << " so skip"
	     << " (" << n_short_v << "/" << igene << " v germlines too short, " << n_long_erosions << "/" << igene << " would require more than 10 erosions)" << endl;
      }
      return;
    }
    
  }

  best_scores_[kset] = AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]));  // i.e. best_prob = v_prob * d_prob * j_prob, v *and* d *and* j
  total_scores_[kset] = AddWithMinusInfinities(regional_total_scores["v"], AddWithMinusInfinities(regional_total_scores["d"], regional_total_scores["j"]));
  if (algorithm_=="viterbi")
    PushBackRecoEvent(kset, best_genes_[kset], best_scores_[kset]);
  // if (debug_ == 2 && algorithm_=="forward") {
  //   cout << "          ";
  //   for (auto &region: gl_.regions_) printf("%5.2f", regional_total_scores[region]);
  //   for (auto &region: gl_.regions_) printf("  %s", tc.ColorGene(gene).c_str());
  //   cout << endl;
  // }
}

// ----------------------------------------------------------------------------------------
size_t JobHolder::GetInsertLength(vector<string> labels) {
  size_t n_inserts(0);
  for (auto &label: labels) {
    if (label=="i")
      n_inserts++;
  }
  return n_inserts;
}

// ----------------------------------------------------------------------------------------
size_t JobHolder::GetErosionLength(string side, vector<string> labels, string gene_name) {
  string germline(gl_.seqs_[gene_name]);
  // first find the index in <labels> up to which we eroded
  size_t istate(0);  // index (in path) of first non-eroded state
  if (side=="left") { // to get left erosion length we look at the first non-insert state in the path
    if (labels[labels.size()-1] == "i") return floor(float(germline.size()) / 2);  // eroded entire sequence -- can't say how much was left and how much was right, so just divide by two
                                                                                   // TODO here (*and* below) do something better than floor/ceil to divide it up.
    for (size_t il=0; il<labels.size(); ++il) {  // loop over each state from left to right
      if (labels[il] == "i") {  // skip any insert states on the left
	continue;
      } else {  // found the leftmost non-insert state -- that's the one we want
	istate = il;
	break;
      }
    }
  } else if (side=="right") {  // and for the righthand one we need the last non-insert state
    if (labels[labels.size()-1] == "i") return ceil(float(germline.size()) / 2);  // eroded entire sequence -- can't say how much was left and how much was right, so just divide by two
    istate = labels.size()-1;
  } else {
    assert(0);
  }

  // then find the state number (in the hmm's state numbering scheme) of the state found at that index in the viterbi path
  assert(labels[istate].find("IGH") == 0);  // start of state label should be IGH[VDJ]
  string state_index_str = labels[istate].substr(labels[istate].find_last_of("_") + 1);
  size_t state_index = atoi(state_index_str.c_str());

  size_t length(0);
  if (side=="left") {
    length = state_index;
  } else if (side=="right") {
    size_t germline_length = gl_.seqs_[gene_name].size();
    length = germline_length - state_index - 1;
  } else {
    assert(0);
  }

  return length;
}
