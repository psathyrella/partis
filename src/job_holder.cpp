#include "job_holder.h"

// ----------------------------------------------------------------------------------------
model *HMMHolder::Get(string gene) {
  if (hmms_.find(gene) == hmms_.end()) {  // if we don't already have it, read it from disk
    hmms_[gene] = new model;
    hmms_[gene]->import(hmm_dir_ + "/" + gl_.GetRegion(gene) + "/" + gl_.SanitizeName(gene) + ".hmm");
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
JobHolder::JobHolder(size_t n_seqs_per_track, string algorithm, sequences *seqs, HMMHolder *hmms, bool debug, string only_genes):
  n_seqs_per_track_(n_seqs_per_track),
  algorithm_(algorithm),
  seqs_(seqs),
  hmms_(hmms),
  debug_(debug)
{
  while (true) {
    size_t i_next_colon(only_genes.find(":"));
    if (i_next_colon == string::npos)
      break;
    only_genes_.insert(only_genes.substr(0,i_next_colon));  // get the next gene name
    only_genes = only_genes.substr(i_next_colon+1);  // then excise it from only_genes
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
  map<string,sequences> subseqs;
  for (auto &region: gl_.regions_)
    subseqs[region] = GetSubSeqs(kset, region);
  return subseqs;
}

// ----------------------------------------------------------------------------------------
void JobHolder::Run(size_t k_v_start, size_t n_k_v, size_t k_d_start, size_t n_k_d) {
  double best_score(-INFINITY),total_score(-INFINITY);
  KSet best_kset;
  if (debug_)
    cout
      << "  k_v: " << k_v_start << "-" << k_v_start+n_k_v-1
      << "  k_d: " << k_d_start << "-" << k_d_start+n_k_d-1
      << endl;
  for (size_t k_v=k_v_start; k_v<k_v_start+n_k_v; ++k_v) {
    for (size_t k_d=k_d_start; k_d<k_d_start+n_k_d; ++k_d) {
      // k_v = 296; k_d = 17;
      KSet kset(k_v,k_d);
      RunKSet(kset);
      total_score = AddInLogSpace(total_scores_[kset], total_score);  // sum up the probabilities for each kset, log P_tot = log \sum_i P_k_i
      if (best_scores_[kset] > best_score) {
      	best_score = best_scores_[kset];
      	best_kset = kset;
      }
    }
  }
  if (debug_)
    if (algorithm_=="viterbi")
      cout << "  best: " << setw(4) << best_kset.first << setw(4) << best_kset.second << setw(12) << best_score << endl;
    else
      cout << "  total: " << total_score << endl;
  if (algorithm_ == "viterbi")
    FillRecoEvent(best_kset, best_genes_[best_kset]);
  else
    cout << "score " << best_score << endl;
  if (best_kset.first == k_v_start ||
      best_kset.first == k_v_start+n_k_v-1 ||
      best_kset.second == k_d_start ||
      best_kset.second == k_d_start+n_k_d-1)
    cout << "\nWARNING you're at the boundary of the specified k_v k_d space" << endl
	 << "  k_v: " << best_kset.first << "(" << k_v_start << "-" << k_v_start+n_k_v-1 << ")" << endl
	 << "  k_d: " << best_kset.second << "(" << k_d_start << "-" << k_d_start+n_k_d-1 << ")" << endl;
}

// ----------------------------------------------------------------------------------------
void JobHolder::FillTrellis(sequences *query_seqs, StrPair query_strs, string gene, double *score) {
  if (trellisi_.find(gene) != trellisi_.end() &&
      trellisi_[gene].find(query_strs) != trellisi_[gene].end()) {  // if we already did this gene for this query sequence. NOTE if we have one gene for this <query_strs>, we're always gonna have the rest of 'em too
    // um, don't do anything a.t.m.
    if (algorithm_=="viterbi")
      *score = paths_[gene][query_strs]->getScore();
    else if (algorithm_=="forward")
      *score = trellisi_[gene][query_strs]->getForwardProbability();
    else
      assert(0);
  } else {
    if (trellisi_.find(gene) == trellisi_.end()) {
      trellisi_[gene] = map<StrPair,trellis*>();
      paths_[gene] = map<StrPair,traceback_path*>();
    }
    trellisi_[gene][query_strs] = new trellis(hmms_->Get(gene), query_seqs);
    if (algorithm_=="viterbi") {
      trellisi_[gene][query_strs]->viterbi();
      paths_[gene][query_strs] = new traceback_path(hmms_->Get(gene));
      trellisi_[gene][query_strs]->traceback(*paths_[gene][query_strs]);
      *score = paths_[gene][query_strs]->getScore();
      if (trellisi_[gene][query_strs]->ending_viterbi_score == -INFINITY) *score = -INFINITY;  // no valid path through hmm. TODO fix this in a more general way
      if (debug_)
	PrintPath(query_strs, gene);
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
void JobHolder::PrintPath(StrPair query_strs, string gene) {  // NOTE query_str is seq1xseq2 for pair hmm
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
    modified_seq = modified_seq + "i";
  assert(modified_seq.size() == query_strs.first.size());
  assert(germline.size() + insert_length - left_erosion_length - right_erosion_length == query_strs.first.size());
  TermColors tc;
  cout
    << "                    " << (left_erosion_length>0 ? ".." : "  ") << tc.ColorMutants("red", query_strs.first, modified_seq, query_strs.second) << (right_erosion_length>0 ? ".." : "  ")
    << setw(12) << paths_[gene][query_strs]->getScore()
    << setw(25) << gene
    << endl;
}

// ----------------------------------------------------------------------------------------
void JobHolder::FillRecoEvent(KSet kset, map<string,string> &best_genes) {
  event_.Clear();
  StrPair seqs;
  for (auto &region: gl_.regions_) {
    StrPair query_strs(GetQueryStrs(kset, region));
    assert(best_genes.find(region) != best_genes.end());
    string gene(best_genes[region]);
    vector<string> path_labels;
    paths_[gene][query_strs]->label(path_labels);
    if (path_labels.size() == 0) {
      if (debug_) cout << "                     " << gene << " has no valid path" << endl;
      return; // TODO fix this upstream. well, it isn't *broken*, but, you know, whatever
    }
    assert(path_labels.size() > 0);
    assert(path_labels.size() == query_strs.first.size());
    event_.SetGene(region, gene);

    // if (region=="v" || region=="d" || region=="j")  // set right-hand deletions
    event_.SetDeletion(region + "_3p", GetErosionLength("right", path_labels, gene));
    // if (region=="v" || region=="d" || region=="j")  // and left-hand deletions
    event_.SetDeletion(region + "_5p", GetErosionLength("left", path_labels, gene));
    if (region=="v")
      event_.SetInsertion("vd", query_strs.first.substr(query_strs.first.size() - GetInsertLength(path_labels)));
    if (region=="d")
      event_.SetInsertion("dj", query_strs.first.substr(query_strs.first.size() - GetInsertLength(path_labels)));
    seqs.first += query_strs.first;
    seqs.second += query_strs.second;
  }
  event_.SetSeq(seqs.first);
  event_.Print(gl_);
  if (n_seqs_per_track_ == 2) {
    event_.SetSeq(seqs.second);
    event_.Print(gl_, 0, 0, true);
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
// add two numbers, treating -INFINITY as zero, i.e. calculates log a*b = log a + log b
double JobHolder::AddWithMinusInfinities(double first, double second) {
  if (first == -INFINITY)
    return second;
  else if(second == -INFINITY)
    return first;
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
  if (debug_)
    cout << "            " << kset.first << " " << kset.second << " -------------------" << endl;
  for (auto &region : gl_.regions_) {
    StrPair query_strs(GetQueryStrs(kset, region));
    TermColors tc;
    if (debug_) {
      cout << "              " << region << " query " << tc.ColorMutants("purple", query_strs.second, query_strs.first) << endl;
      if (n_seqs_per_track_==2)
	cout << "              " << region << " query " << tc.ColorMutants("purple", query_strs.first, query_strs.second) << endl;
    }

    regional_best_scores[region] = -INFINITY;
    regional_total_scores[region] = -INFINITY;
    size_t igene(0),n_short_j(0);
    for (auto &gene : gl_.names_[region]) {
      if (only_genes_.size()>0 and only_genes_.find(gene)==only_genes_.end())
	continue;
      igene++;
      if (region=="j" && query_strs.first.size()>gl_.seqs_[gene].size()) {  // query sequence too long for this j version to make any sense
	if (debug_) cout << "                     " << gene << " too short" << endl;
	n_short_j++;
	continue;
      }
      double gene_score(-INFINITY);
      FillTrellis(&subseqs[region], query_strs, gene, &gene_score);
      regional_total_scores[region] = AddInLogSpace(gene_score, regional_total_scores[region]);  // (log a, log b) --> log a+b, i.e. here we are summing probabilities in log space
      if (gene_score > regional_best_scores[region]) {
	regional_best_scores[region] = gene_score;
	best_genes_[kset][region] = gene;
      }
    }
    if (best_genes_[kset].find(region) == best_genes_[kset].end()) {
      if (debug_)
	cout << "                  found no gene for " << region
	     << " (" << n_short_j << " / " << igene << " j versions were too short for the query sequence)" << endl;
      return;
    }
  }
  if (debug_ && algorithm_=="viterbi")
    FillRecoEvent(kset, best_genes_[kset]);
  best_scores_[kset] = AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]));  // i.e. best_prob = v_prob * d_prob * j_prob
  total_scores_[kset] = AddWithMinusInfinities(regional_total_scores["v"], AddWithMinusInfinities(regional_total_scores["d"], regional_total_scores["j"]));
  if (debug_)
    cout
      << "        "
      << " " << best_genes_[kset]["v"]
      << " " << best_genes_[kset]["d"]
      << " " << best_genes_[kset]["j"]
      << "  --> "
      << AddWithMinusInfinities(regional_best_scores["v"], AddWithMinusInfinities(regional_best_scores["d"], regional_best_scores["j"]))
      << endl;
}

// // ----------------------------------------------------------------------------------------
// size_t JobHolder::GetInsertion(vector<string> labels, string gene_name) {
//   // find index (in hmm state numbering scheme) of last non-insertion state
//   size_t i_last_germline(0);
//   for (auto &state: labels) {
//     if (state!="i") { // if we're still among germline states, set i_last_germline to the current state
//       i_last_germline = atoi(state.substr(state.find_last_of("_")+1));
//     } else { // otherwise break -- i_last_germline should be set to the last germline state
//       break;
//     }
//   }
//   assert(i_last_germline>0);
//   string germline(gl_.seqs[gene]);
//   assert(i_last_germline<germline.size());
  
  
// }
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
  // first find the index in <labels> up to which we eroded
  size_t istate;
  if (side=="left") { // to get left erosion length we look at the first state in the path
    if (labels[0] == "i") return 0;  // eroded entire sequence
    istate = 0;
  } else if (side=="right") {  // and for the righthand one we need the last non-insert state
    for (size_t il=labels.size()-1; il>=0; --il) {  // loop over each state from the right
      if (labels[il] == "i") {  // skip any insert states on the right
	continue;
      } else {  // found the rightmost non-insert state -- that's the one we want
	istate = il;
	break;
      }
    }
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
