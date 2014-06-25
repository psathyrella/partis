#include "job_holder.h"

// ----------------------------------------------------------------------------------------
model *HMMHolder::Get(string gene) {
  if (hmms_.find(gene) == hmms_.end()) {  // if we don't already have it, read it from disk
    hmms_[gene] = new model;
    hmms_[gene]->import(hmm_dir_ + "/" + gl_.GetRegion(gene) + "/" + gl_.SanitizeName(gene) + ".hmm");
  }
  return hmms_[gene];
}  

// ----------------------------------------------------------------------------------------
HMMHolder::~HMMHolder() {
  for (auto &entry: hmms_)
    delete entry.second;
}

// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
JobHolder::JobHolder(string hmmtype, string algorithm, sequences *seqs, HMMHolder *hmms, size_t n_max_versions):
  n_max_versions_(n_max_versions),
  algorithm_(algorithm),
  seqs_(seqs),
  hmms_(hmms),
  debug_(false)
{
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
map<string,sequences> JobHolder::GetSubSeqs(KSet kset) {
  size_t k_v(kset.first),k_d(kset.second);
  map<string,sequences> subseqs;
  subseqs["v"] = seqs_->getSubSequences(0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  subseqs["d"] = seqs_->getSubSequences(k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  subseqs["j"] = seqs_->getSubSequences(k_v + k_d, seqs_->size() - k_v - k_d);  // j region runs from k_v + k_d to end
  return subseqs;
}

// ----------------------------------------------------------------------------------------
void JobHolder::Run(size_t k_v_start, size_t n_k_v, size_t k_d_start, size_t n_k_d) {
  double best_score(-INFINITY);
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
      if (scores_[kset] > best_score) {
      	best_score = scores_[kset];
      	best_kset = kset;
      }
    }
  }
  if (debug_)
    cout << "  best: " << setw(4) << best_kset.first << setw(4) << best_kset.second << setw(12) << best_score << endl;
  map<string,string> query_strs(GetQueryStrs(best_kset));
  FillRecoEvent(best_genes_[best_kset], query_strs);
  if (best_kset.first == k_v_start ||
      best_kset.first == k_v_start+n_k_v-1 ||
      best_kset.second == k_d_start ||
      best_kset.second == k_d_start+n_k_d-1)
    cout << "\nWARNING you're at the boundary of the specified k_v k_d space" << endl
	 << "  k_v: " << best_kset.first << "(" << k_v_start << "-" << k_v_start+n_k_v-1 << ")" << endl
	 << "  k_d: " << best_kset.second << "(" << k_d_start << "-" << k_d_start+n_k_d-1 << ")" << endl;
 }

// ----------------------------------------------------------------------------------------
void JobHolder::FillTrellis(sequences *query_seqs, string region, string gene) {
  assert(algorithm_=="viterbi");  // de-implemented the others for the moment
  string query_str((*query_seqs)[0].undigitize());
  if (trellisi_.find(gene) != trellisi_.end() &&
      trellisi_[gene].find(query_str) != trellisi_[gene].end()) {  // if we already did this gene for this query sequence. NOTE if we have one gene for this query_str, we're always gonna have the rest of 'em too
  } else {
    if (trellisi_.find(gene) == trellisi_.end()) {
      trellisi_[gene] = map<string,trellis*>();
      paths_[gene] = map<string,traceback_path*>();
    }
    trellisi_[gene][query_str] = new trellis(hmms_->Get(gene), query_seqs);
    trellisi_[gene][query_str]->viterbi();
    paths_[gene][query_str] = new traceback_path(hmms_->Get(gene));
    trellisi_[gene][query_str]->traceback(*paths_[gene][query_str]);
    if (debug_)
      PrintPath(query_str, gene);
  }
  // } else if (algorithm_=="forward") {
  // 	trellisi_[kset][gene]->forward();
  // 	best_region_score = trellisi_[kset][gene]->getForwardProbability();
  // }
}

// ----------------------------------------------------------------------------------------
void JobHolder::PrintPath(string query_str, string gene) {
  vector<string> path_labels;
  paths_[gene][query_str]->label(path_labels);
  assert(path_labels.size() > 0);
  assert(path_labels.size() == query_str.size());
  size_t insert_length = GetInsertLength(path_labels);
  size_t left_erosion_length = GetErosionLength("left", path_labels, gene);
  size_t right_erosion_length = GetErosionLength("right", path_labels, gene);

  string germline(gl_.seqs_[gene]);
  string modified_seq = germline.substr(left_erosion_length, germline.size() - right_erosion_length - left_erosion_length);
  for (size_t i=0; i<insert_length; ++i)
    modified_seq = modified_seq + "i";
  assert(modified_seq.size() == query_str.size());
  assert(germline.size() + insert_length - left_erosion_length - right_erosion_length == query_str.size());
  TermColors tc;
  cout
    << "                    " << (left_erosion_length>0 ? ".." : "  ") << tc.redify_mutants(query_str, modified_seq) << (right_erosion_length>0 ? ".." : "  ")
    << setw(12) << paths_[gene][query_str]->getScore()
    << setw(25) << gene
    << endl;
}

// ----------------------------------------------------------------------------------------
void JobHolder::FillRecoEvent(map<string,string> &best_genes, map<string,string> &query_strs) {
  event_.Clear();
  for (auto &region: gl_.regions_) {
    string query_str(query_strs[region]);
    assert(best_genes.find(region) != best_genes.end());
    string gene(best_genes[region]);
    vector<string> path_labels;
    paths_[gene][query_str]->label(path_labels);
    assert(path_labels.size() > 0);
    assert(path_labels.size() == query_str.size());
    event_.SetGene(region, gene);

    if (region=="v" || region=="d" || region=="j")  // set right-hand deletions
      event_.SetDeletion(region + "_3p", GetErosionLength("right", path_labels, gene));
    if (region=="v" || region=="d" || region=="j")  // and left-hand deletions
      event_.SetDeletion(region + "_5p", GetErosionLength("left", path_labels, gene));
    if (region=="v")
      event_.SetInsertion("vd", query_str.substr(query_str.size() - GetInsertLength(path_labels)));
    if (region=="d")
      event_.SetInsertion("dj", query_str.substr(query_str.size() - GetInsertLength(path_labels)));
  }
  event_.SetSeq(query_strs["v"] + query_strs["d"] + query_strs["j"]);
  event_.Print(gl_);
}

// ----------------------------------------------------------------------------------------
map<string,string> JobHolder::GetQueryStrs(KSet kset) {
  map<string,sequences> subseqs = GetSubSeqs(kset);
  map<string,string> query_strs;
  for (auto &region: gl_.regions_) {
    sequences *query_seqs(&subseqs[region]);
    query_strs[region] = ((*query_seqs)[0].undigitize());
  }
  return query_strs;  
}
// ----------------------------------------------------------------------------------------
void JobHolder::RunKSet(KSet kset) {
  map<string,sequences> subseqs = GetSubSeqs(kset);
  map<string,string> query_strs(GetQueryStrs(kset));
  map<string,double> best_scores;
  for (auto &region : gl_.regions_) {
    sequences *query_seqs(&subseqs[region]);
    assert(query_seqs->size() == 1);  // testing stuff a.t.m. -- no pair hmms
    if (debug_)
      cout << "              " << region << " query " << query_strs[region] << endl;

    best_scores[region] = -INFINITY;
    size_t igene(0),n_short_j(0);
    for (auto &gene : gl_.names_[region]) {
      if (igene >= n_max_versions_)
	break;
      igene++;
      if (region=="j" && query_strs[region].size()>gl_.seqs_[gene].size()) {  // query sequence too long for this j version to make any sense
	if (debug_) cout << "                     " << gene << " too short" << endl;
	n_short_j++;
	continue;
      }
      FillTrellis(query_seqs, region, gene);
      if (paths_[gene][query_strs[region]]->getScore() > best_scores[region]) {
	best_scores[region] = paths_[gene][query_strs[region]]->getScore();
	best_genes_[kset][region] = gene;
      }
    }
    if (best_genes_[kset].find(region) == best_genes_[kset].end()) {
      if (debug_)
	cout << "        found no gene for " << region
	     << " (" << n_short_j << " / " << min(n_max_versions_, gl_.names_["j"].size()) << " j versions were too short for the query sequence)" << endl;
      scores_[kset] = -INFINITY;
      return;
    }
  }
  if (debug_)
    FillRecoEvent(best_genes_[kset], query_strs);
  scores_[kset] = best_scores["v"] + best_scores["d"] + best_scores["j"];  // set total score for this kset
  if (debug_)
    cout
      << "        "
      << " " << best_genes_[kset]["v"]
      << " " << best_genes_[kset]["d"]
      << " " << best_genes_[kset]["j"]
      << "  --> "
      << best_scores["v"] + best_scores["d"] + best_scores["j"]
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
