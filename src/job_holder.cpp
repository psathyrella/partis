#include "job_holder.h"

// ----------------------------------------------------------------------------------------
JobHolder::JobHolder(string hmmtype, string hmm_dir, string seqfname, size_t n_max_versions):
  hmm_dir_(hmm_dir),
  n_max_versions_(n_max_versions)
{
  vector<string> characters{"A","C","G","T"};
  track_ = track("NUKES", hmmtype == "single" ? 1 : 2, characters);

  // read sequences
  ifstream ifss(seqfname);
  assert(ifss.is_open());
  for (size_t iseq=0; iseq<track_.n_seqs; ++iseq) {  // for pair hmm, we push back two seqs for each track
    sequence *sq = new(nothrow) sequence(false);
    assert(sq);
    sq->getFasta(ifss, &track_);
    seqs_.addSeq(sq);
  }
  ifss.close();
}

// ----------------------------------------------------------------------------------------
map<string,sequences> JobHolder::GetSubSeqs(size_t k_v, size_t k_d) {
  map<string,sequences> subseqs;
  subseqs["v"] = seqs_.getSubSequences(0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  subseqs["d"] = seqs_.getSubSequences(k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  subseqs["j"] = seqs_.getSubSequences(k_v + k_d, seqs_.size() - k_v - k_d);  // j region runs from k_v + k_d to end
  return subseqs;
}

// ----------------------------------------------------------------------------------------
void JobHolder::Run(size_t k_v, size_t k_d, string algorithm) {
  map<string,sequences> subseqs = GetSubSeqs(k_v, k_d);
  double total_score(-INFINITY);
  for (auto &region : gl_.regions_) {
    double best_region_score(-INFINITY);
    string best_gene;
    size_t igene(0);
    sequences *query_seqs(&subseqs[region]);
    assert(query_seqs->size() == 1);  // testing stuff a.t.m. -- no pair hmms
    string query_str((*query_seqs)[0].undigitize());
    cout << "                query " << query_str << endl;
    for (auto &gene : gl_.names_[region]) {
      if (igene >= n_max_versions_)
	break;
      igene++;
      assert(gl_.seqs_.find(gene) != gl_.seqs_.end());
      string germline(gl_.seqs_[gene]);
      model hmm;
      hmm.import(hmm_dir_ + "/" + region + "/" + gl_.SanitizeName(gene) + ".hmm");
      if (region=="j" && query_str.size()>germline.size())  // query sequence too long for this j version to make any sense
	continue;
      trellis trell(&hmm, query_seqs);
      if (algorithm=="viterbi") {
	// run viterbi
	trell.viterbi();
	// trace back the path
	traceback_path path(&hmm);
	trell.traceback(path);
	vector<string> path_labels;
	path.label(path_labels);
	assert(path_labels.size() > 0);
	if (path.getScore() > best_region_score) {
	  best_region_score = path.getScore();
	  best_gene = gene;
	}
	size_t insert_length = GetInsertLength(path_labels);
	size_t left_erosion_length = GetErosionLength("left", path_labels, gene);
	size_t right_erosion_length = GetErosionLength("right", path_labels, gene);
	// assert(query_str

	assert(query_str.size() == path_labels.size());

	string modified_seq = germline.substr(left_erosion_length, germline.size() - right_erosion_length - left_erosion_length);
	// for (size_t i=0; i<left_erosion_length; ++i)
	//   modified_seq = "i" + modified_seq;
	for (size_t i=0; i<insert_length; ++i)
	  modified_seq = modified_seq + "i";
	assert(modified_seq.size() == query_str.size());
	assert(germline.size() + insert_length - left_erosion_length - right_erosion_length == query_str.size());
	TermColors tc;
	cout
	  << "                      " << tc.redify_mutants(query_str, modified_seq)
	  << setw(12) << path.getScore()
	  << setw(25) << gene
	  << endl;
      } else if (algorithm=="forward") {
	trell.forward();
	best_region_score = trell.getForwardProbability();
      }
    }
    cout << "  " << setw(12) << best_region_score << setw(29) << best_gene << endl;
    if (total_score == -INFINITY)
      total_score = best_region_score;
    else
      total_score += best_region_score;
  }
  cout << "  overall " << total_score << "  for k set " << k_v << " " << k_d << endl;
}

// ----------------------------------------------------------------------------------------
size_t JobHolder::GetInsertLength(vector<string> labels) {
  size_t n_inserts(0);
  for (auto &label : labels) {
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
