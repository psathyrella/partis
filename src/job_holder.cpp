#include "job_holder.h"

// ----------------------------------------------------------------------------------------
JobHolder::JobHolder(string hmmtype, string hmm_dir, string seqfname):
  finished_(false),
  hmm_dir_(hmm_dir)
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
// chop up query sequences into subseqs in preparation for running them through each hmm
void JobHolder::InitJobs(size_t k_v, size_t k_d) {
  subseqs_ = GetSubSeqs(k_v, k_d);
  i_current_region_ = gl_.regions_.begin();
  i_current_gene_ = gl_.names_[*i_current_region_].begin();
}

// ----------------------------------------------------------------------------------------
// get the next job to process, i.e. the next germline hmm through which to run the query sequence
void JobHolder::GetNextHMM() {
  assert(subseqs_.size() > 0);  // make sure InitJobs was run first

  current_seqs_ = &subseqs_[*i_current_region_];
  current_hmm_ = model();
  current_hmm_.import(hmm_dir_ + "/" + (*i_current_region_) + "/" + sanitize_name(*i_current_gene_) + ".hmm");
    
  cout
    << setw(12) << (*i_current_region_)
    << setw(40) << (*i_current_gene_);

  // set iterators to the *next* ones, or else set finished_ to true if we're done
  i_current_gene_++;
  if (i_current_gene_ == gl_.names_[*i_current_region_].end()) {
     i_current_region_++;
     if (i_current_region_ == gl_.regions_.end()) {
       finished_ = true;
     } else {
       i_current_gene_ = gl_.names_[*i_current_region_].begin();
     }
  }
}

// ----------------------------------------------------------------------------------------
// replace * with _star_ and / with _slash_
string JobHolder::sanitize_name(string gene_name) {
  size_t istar(gene_name.find("*"));
  if (istar != string::npos)
    gene_name.replace(istar, 1, "_star_");

  size_t islash(gene_name.find("/"));
  if (islash != string::npos)
    gene_name.replace(islash, 1, "_slash_");

  return gene_name;
}
