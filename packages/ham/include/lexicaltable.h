#ifndef HAM_LEXICALTABLE_H
#define HAM_LEXICALTABLE_H

#include <string>
#include <vector>
#include <cassert>
#include "track.h"
#include "sequences.h"

using namespace std;
namespace ham {

class LexicalTable {
public:
  LexicalTable();
  void Init();
  ~LexicalTable();

  void AddTrack(Track* trk, int order) { tracks.push_back(trk); }
  void AddColumn(vector<double> logprobs);

  inline double LogProb(size_t letter) { assert(letter < (*log_probs_)[0].size()); return (*log_probs_)[0][letter]; }
  inline double LogProb(size_t letter1, size_t letter2) { assert(letter1 < log_probs_->size()); assert(letter2 < (*log_probs_)[letter1].size()); return (*log_probs_)[letter1][letter2]; }
  inline double LogProb(Sequences &seqs, size_t pos) { return (*log_probs_)[seqs[0][pos]][seqs[1][pos]]; }  // NOTE <seqs> *must* have length of two, and init() *must* have been called. I could check this, but I'm prematurely optimising. Good thing I'm not NASA, eh?
  inline double LogProb(Sequence &seq, size_t pos) { return (*log_probs_)[0][seq[pos]]; }  // see NOTE above
  inline Track* track(size_t iter) { return tracks.at(iter); }
  inline size_t n_tracks(){ return tracks.size(); }
  inline uint8_t alphabet_size(size_t i){ return tracks[i]->alphabet_size(); }
private:
  vector<Track*> tracks;  // tracks which are used by emissions in this table
  // log_probs_ scheme:
  //   first index: seq 1 emission
  //   second index: seq 2 emission
  //   NOTE the two sequences cannot (in current implementation) be distinguishable
  vector<vector<double> >* log_probs_;
};

}
#endif
