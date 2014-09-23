#include "job_holder.h"

// ----------------------------------------------------------------------------------------
model *HMMHolder::Get(string gene) {
  if (hmms_.find(gene) == hmms_.end()) {  // if we don't already have it, read it from disk
    hmms_[gene] = new model;
    hmms_[gene]->parse(hmm_dir_ + "/" + gl_.SanitizeName(gene) + ".hmm");
  }
  return hmms_[gene];
}  

// ----------------------------------------------------------------------------------------
HMMHolder::~HMMHolder() {
  for (auto &entry: hmms_)
    delete entry.second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JobHolder::JobHolder(GermLines &gl, HMMHolder &hmms, string algorithm, string only_gene_str):
  gl_(gl),
  hmms_(hmms),
  algorithm_(algorithm),
  // total_score_(-INFINITY),
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
  Clear();
}

// ----------------------------------------------------------------------------------------
void JobHolder::Clear() {
  for (auto &gene_map: trellisi_) {
    string gene(gene_map.first);
    for (auto &query_str_map: gene_map.second) {
      delete query_str_map.second;
      if (paths_[gene][query_str_map.first])  // set to nullptr if no valid path
	delete paths_[gene][query_str_map.first];
    }
  }
  trellisi_.clear();
  paths_.clear();
  all_scores_.clear();
  best_per_gene_scores_.clear();
  errors_ = "";
}

// ----------------------------------------------------------------------------------------
sequences JobHolder::GetSubSeqs(sequences &seqs, KSet kset, string region) {
  // get subsequences for one region
  size_t k_v(kset.first),k_d(kset.second);
  if (region=="v")
    return seqs.getSubSequences(0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  else if (region=="d")
    return seqs.getSubSequences(k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  else if (region=="j")
    return seqs.getSubSequences(k_v + k_d, seqs.GetSequenceLength() - k_v - k_d);  // j region runs from k_v + k_d to end
  else
    assert(0);
}

// ----------------------------------------------------------------------------------------
map<string,sequences> JobHolder::GetSubSeqs(sequences &seqs, KSet kset) {
  // get subsequences for all regions
  map<string,sequences> subseqs;
  for (auto &region: gl_.regions_)
    subseqs[region] = GetSubSeqs(seqs, kset, region);
  return subseqs;
}

// ----------------------------------------------------------------------------------------
Result JobHolder::Run(sequence &seq, size_t k_v_min, size_t k_v_max, size_t k_d_min, size_t k_d_max) {
  sequences seqs;
  sequence *newseq = new sequence(seq);  // seriously wtf does <sequences> need to own its sequences?
  seqs.addSeq(newseq);
  return Run(seqs, k_v_min, k_v_max, k_d_min, k_d_max);
}

// ----------------------------------------------------------------------------------------
Result JobHolder::Run(sequences &seqs, size_t k_v_min, size_t k_v_max, size_t k_d_min, size_t k_d_max) {
  assert(k_v_max>k_v_min && k_d_max>k_d_min);
  assert(k_v_min > 0 && k_d_min > 0);
  assert(seqs.n_seqs() == 1 || seqs.n_seqs() == 2);
  Clear();
  assert(trellisi_.size()==0 && paths_.size()==0 && all_scores_.size()==0);
  map<KSet,double> best_scores;  // best score for each kset (summed over regions)
  map<KSet,double> total_scores;  // total score for each kset (summed over regions)
  map<KSet,map<string,string> > best_genes;  // map from a kset to its corresponding triplet of best genes

  Result result;

  // loop over k_v k_d space
  double best_score(-INFINITY);
  KSet best_kset(make_pair(0,0));
  double *total_score = &result.total_score_;  // total score for all ksets
  for (size_t k_v=k_v_min; k_v<k_v_max; ++k_v) {
    for (size_t k_d=k_d_min; k_d<k_d_max; ++k_d) {
      if (k_v + k_d >= seqs.GetSequenceLength()) {
	cout << "      skipping " << k_v << " + " << k_d << " = " << k_v + k_d << " >= " << seqs.GetSequenceLength() << endl;
	continue;
      }
      KSet kset(k_v,k_d);
      RunKSet(seqs, kset, &best_scores, &total_scores, &best_genes);
      *total_score = AddInLogSpace(total_scores[kset], *total_score);  // sum up the probabilities for each kset, log P_tot = log \sum_i P_k_i
      if (debug_==2 && algorithm_=="forward") printf("            %9.2f (%.1e)  tot: %7.2f\n", total_scores[kset], exp(total_scores[kset]), *total_score);
      if (best_scores[kset] > best_score) {
      	best_score = best_scores[kset];
      	best_kset = kset;
      }
      if (algorithm_=="viterbi" && best_scores[kset] != -INFINITY)
	PushBackRecoEvent(seqs, kset, best_genes[kset], best_scores[kset], &result.events_);
    }
  }

  // return if no valid path
  if (best_kset.first == 0) {
    if (debug_) cout << "  nothing sensical" << endl;
    return result;
  }

  // sort vector of events by score, stream info to stderr, and print the top n_best_events_
  if (algorithm_ == "viterbi") {
    sort(result.events_.begin(), result.events_.end());
    reverse(result.events_.begin(), result.events_.end());
    if (debug_==2) {
      assert(n_best_events_ <= result.events_.size());
      for (size_t ievt=0; ievt<n_best_events_; ++ievt) {
	result.events_[ievt].Print(gl_, 0, 0, false, false, "          ");  // man, I wish I had keyword args
	if (seqs.n_seqs() == 2)
	  result.events_[ievt].Print(gl_, 0, 0, true, true, "          ");
      }
    }
  }
  // StreamOutput(total_score_);  // NOTE this must happen after sorting in viterbi

  // warn if we're on the k_v k_d space boundary
  if (k_v_max-k_v_min > 1 && k_d_max-k_d_min > 2) {  // for the stripped-down hmm we kinda expect the max to be on the boundary
    if (best_kset.first == k_v_min)
      result.boundary_error_ = "k_v_min";
    else if (best_kset.first == k_v_max-1)
      result.boundary_error_ = "k_v_max";
    else if (best_kset.second == k_d_min)
      result.boundary_error_ = "k_d_min";
    else if (best_kset.second == k_d_max-1)
      result.boundary_error_ = "k_d_max";
    else
      result.boundary_error_ = "";

    if (result.boundary_error_ != "") {
      errors_ += ":boundary";
      if (debug_) {
	cout << "              WARNING maximum at boundary for "
	     << seqs[0].name();
	if (seqs.n_seqs() == 2)
	  cout << " " << seqs[1].name();
	cout
	  << "  k_v: " << best_kset.first << "(" << k_v_min << "-" << k_v_max-1 << ")"
	  << "  k_d: " << best_kset.second << "(" << k_d_min << "-" << k_d_max-1 << ")" << " " << result.boundary_error_ << endl;
	// cout << "ok, I changed my mind, ERROR!" << endl;
      }
      // assert(0);  TODO make sure it doesn't happen
    }
  }
  
  // print debug info
  if (debug_) {
    cout << "    " << setw(22) << seqs[0].name() << " " << setw(22) << (seqs.n_seqs()==2 ? seqs[1].name() : "") << "   " << k_v_min << "-" << k_v_max-1 << "   " << k_d_min << "-" << k_d_max-1;  // exclusive...
    if (algorithm_=="viterbi")
      cout << "    best kset: " << setw(4) << best_kset.first << setw(4) << best_kset.second << setw(12) << best_score << endl;
    else
      cout << "        " << *total_score << endl;
  }

  return result;
}

// ----------------------------------------------------------------------------------------
void JobHolder::FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score) {
  *score = -INFINITY;
  // initialize trellis and path
  if (trellisi_.find(gene) == trellisi_.end()) {
    trellisi_[gene] = map<StrPair,trellis*>();
    paths_[gene] = map<StrPair,traceback_path*>();
  }
  trellisi_[gene][query_strs] = new trellis(hmms_.Get(gene), query_seqs);
  trellis *trell(trellisi_[gene][query_strs]);

  if (algorithm_=="viterbi") {
    trell->viterbi();
    *score = trell->ending_viterbi_score;  // NOTE still need to add the gene choice prob to this score (it's done in RunKSet)
    if (trell->ending_viterbi_score == -INFINITY) {  // no valid path through hmm. TODO fix this in a more general way
      paths_[gene][query_strs] = nullptr;
      if (debug_ == 2) cout << "                    " << gene << " " << *score << endl;
    } else {
      paths_[gene][query_strs] = new traceback_path(hmms_.Get(gene));
      trell->traceback(*paths_[gene][query_strs]);
      assert(trell->ending_viterbi_score == paths_[gene][query_strs]->getScore());  // TODO remove this assertion
      if (debug_ == 2) PrintPath(query_strs, gene, *score);
    }
    assert(fabs(*score) > 1e-200);
    assert(*score == -INFINITY || paths_[gene][query_strs]->size() > 0);
  } else if (algorithm_=="forward") {
    trell->forward();
    paths_[gene][query_strs] = nullptr;  // avoids violating the assumption that paths_ and trellisi_ have the same entries
    *score = trell->getForwardProbability();  // NOTE still need to add the gene choice prob to this score
  } else {
    assert(0);
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
void JobHolder::PushBackRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score, vector<RecoEvent> *events) {
  RecoEvent event(FillRecoEvent(seqs, kset, best_genes, score));
  events->push_back(event);
}

// ----------------------------------------------------------------------------------------
RecoEvent JobHolder::FillRecoEvent(sequences &seqs, KSet kset, map<string,string> &best_genes, double score) {
  RecoEvent event;
  StrPair seq_strs;
  for (auto &region: gl_.regions_) {
    StrPair query_strs(GetQueryStrs(seqs, kset, region));
    if (best_genes.find(region) == best_genes.end()) {
      cout << seqs.stringifyWOHeader() << endl;
    }
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
    seq_strs.first += query_strs.first;
    seq_strs.second += query_strs.second;
  }

  // if (seqs_->size() > 0)  // erm, this segfaults a.t.m. I must be forgetting something somewhere else
  //   assert((*seqs_)[0].name() == (*seqs_)[1].name());  // er, another somewhat neurotic consistency check

  event.SetSeq(seqs[0].name(), seq_strs.first);
  if (seqs.n_seqs() == 2) {
    // assert((*seqs_)[0].name() == (*seqs_)[1].name());  don't recall at this point precisely why it was that I wanted this here
    event.SetSecondSeq(seqs[1].name(), seq_strs.second);
  }
  event.SetScore(score);
  return event;
}

// ----------------------------------------------------------------------------------------
StrPair JobHolder::GetQueryStrs(sequences &seqs, KSet kset, string region) {
  sequences query_seqs(GetSubSeqs(seqs, kset,region));
  StrPair query_strs;
  query_strs.first = query_seqs[0].undigitize();
  if (query_seqs.n_seqs() == 2) {  // the sequences class should already ensure that both seqs are the same length
    assert(seqs.n_seqs() == 2);
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
 void JobHolder::RunKSet(sequences &seqs, KSet kset, map<KSet,double> *best_scores, map<KSet,double> *total_scores, map<KSet,map<string,string> > *best_genes) {
  map<string,sequences> subseqs(GetSubSeqs(seqs, kset));
  (*best_scores)[kset] = -INFINITY;
  (*total_scores)[kset] = -INFINITY;  // total log prob of this kset, i.e. log(P_v * P_d * P_j), where e.g. P_v = \sum_i P(v_i k_v)
  (*best_genes)[kset] = map<string,string>();
  map<string,double> regional_best_scores;  // the best score for each region
  map<string,double> regional_total_scores;  // the total score for each region, i.e. log P_v
  if (debug_ == 2)
    cout << "            " << kset.first << " " << kset.second << " -------------------" << endl;
  for (auto &region : gl_.regions_) {
    StrPair query_strs(GetQueryStrs(seqs, kset, region));

    TermColors tc;
    if (debug_ == 2) {
      if (algorithm_=="viterbi") {
	cout << "              " << region << " query " << tc.ColorMutants("purple", query_strs.second, query_strs.first) << endl;
	if (seqs.n_seqs() == 2)
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
	
      double *gene_score(&all_scores_[gene][query_strs]);
      bool already_cached = trellisi_.find(gene) != trellisi_.end() && trellisi_[gene].find(query_strs) != trellisi_[gene].end();
      if (already_cached) {
	if (debug_==2 && algorithm_=="viterbi")
	  PrintPath(query_strs, gene, *gene_score, "(cached)");
      } else {
	FillTrellis(&subseqs[region], query_strs, gene, gene_score);  // set *gene_score to uncorrected score
	double gene_choice_score = log(hmms_.Get(gene)->overall_gene_prob());  // TODO think through this again, and make sure it's correct for forward score, as well. I mean, I *think* it's right, but I could stand to go over it again
	*gene_score = AddWithMinusInfinities(*gene_score, gene_choice_score);  // then correct it for gene choice probs
      }

      // set regional total scores
      regional_total_scores[region] = AddInLogSpace(*gene_score, regional_total_scores[region]);  // (log a, log b) --> log a+b, i.e. here we are summing probabilities in log space, i.e. a *or* b
      if (debug_ == 2 && algorithm_=="forward")
	printf("                %9.2f (%.1e)  tot: %7.2f  %s\n", *gene_score, exp(*gene_score), regional_total_scores[region], tc.ColorGene(gene).c_str());

      // set best regional scores
      if (*gene_score > regional_best_scores[region]) {
	regional_best_scores[region] = *gene_score;
	(*best_genes)[kset][region] = gene;
      }

      // set best per-gene scores
      if (best_per_gene_scores_.find(gene) == best_per_gene_scores_.end())
	best_per_gene_scores_[gene] = -INFINITY;
      if (*gene_score > best_per_gene_scores_[gene])
	best_per_gene_scores_[gene] = *gene_score;
    }

    // return if we didn't find a valid path for this region
    if ((*best_genes)[kset].find(region) == (*best_genes)[kset].end()) {
      // TODO put this assert back in. Well, I think. At least figure out a better scheme for dealing with v_right_length-related stuff
      // if(n_long_erosions != 0) {
      // 	cout << seqs[0].name() << endl;
      // 	if (seqs.n_seqs() > 1)
      // 	  cout << seqs[1].name() << endl;
      // }
      // assert(n_long_erosions == 0);  // adding assert because if it happens, that means my v_right_length was screwed up
      if (debug_ == 2) {
	cout << "                  found no gene for " << region << " so skip"
	     << " (" << n_short_v << "/" << igene << " v germlines too short, " << n_long_erosions << "/" << igene << " would require more than 10 erosions)" << endl;
      }
      return;
    }
  }

  (*best_scores)[kset] = AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]));  // i.e. best_prob = v_prob * d_prob * j_prob, v *and* d *and* j
  (*total_scores)[kset] = AddWithMinusInfinities(regional_total_scores["v"], AddWithMinusInfinities(regional_total_scores["d"], regional_total_scores["j"]));
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

// ----------------------------------------------------------------------------------------
void JobHolder::WriteBestGeneProbs(ofstream &ofs, string query_name) {
  ofs << query_name << ",";
  stringstream ss;
  for (auto &gene_map: best_per_gene_scores_)
    ss << gene_map.first << ":" << gene_map.second << ";";
  ofs << ss.str().substr(0, ss.str().size()-1) << endl;  // remove the last semicolon
}
