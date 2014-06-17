#include "job_holder.h"

JobHolder::JobHolder(string input_dir, vector<string> regions, string seqfname):regions_(regions) {
  for (auto &region : regions_) {
    hmms_[region] = model();
    hmms_[region].import(input_dir + "/" + region + ".hmm");
    for (size_t itrk=0; itrk<hmms_[region].track_size(); ++itrk) {  // push back a sequence (or more than one) from the input file for each track which was defined in the model
      track *trk = hmms_[region].getTrack(itrk);
      seqs_[region] = sequences();  // init with size zero
      // NOTE this loads the *same* sequences for each hmm -- but since the seqs require a track on import, it's simpler a.t.m. to do it this way
      for (size_t iseq=0; iseq<trk->get_n_seqs(); ++iseq) {  // for eg pair hmm, we push back two seqs for each track
	sequence *sq = new(std::nothrow) sequence(false);
	assert(sq);
	sq->getFasta(seqfname, trk);
	seqs_[region].addSeq(sq);
      }
    }
  }
}

// ----------------------------------------------------------------------------------------
map<string,sequences> JobHolder::GetSubSeqs(size_t k_v, size_t k_d) {
  map<string,sequences> subseqs;
  subseqs["v"] = seqs_["v"].getSubSequences(0, k_v);  // v region (plus vd insert) runs from zero up to k_v
  subseqs["d"] = seqs_["d"].getSubSequences(k_v, k_d);  // d region (plus dj insert) runs from k_v up to k_v + k_d
  subseqs["j"] = seqs_["j"].getSubSequences(k_v + k_d, seqs_["j"].size() - k_v - k_d);  // j region runs from k_v + k_d to end
  // for (auto &region : regions_) {
  //   std::cout << region << "  --------" << std::endl;
  //   std::cout << subseqs[region].stringify() << std::endl;
  // }
  return subseqs;
}
