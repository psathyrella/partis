#ifndef HAM_LEXICALTABLE_H
#define HAM_LEXICALTABLE_H

#include <string>
#include <vector>
#include <math.h>
#include <cassert>
#include "track.h"
#include "sequences.h"

using namespace std;
namespace ham {

class LexicalTable {
public:
  LexicalTable();
  void Init();
  void ReplaceLogProbs(vector<vector<double> > new_log_probs);  // replace values in <log_probs_> with thsoe in <new_log_probs>
  void UnReplaceLogProbs();  // revert to the original log probs
  ~LexicalTable();

  void AddTrack(Track* trk, int order) { tracks.push_back(trk); }
  void AddColumn(vector<double> logprobs);
  vector<vector<double> > log_probs() { return log_probs_; }  // NOTE returns a *copy*

  double LogProb(size_t letter) { assert(letter < log_probs_[0].size()); return log_probs_[0][letter]; }
  double LogProb(Sequence *seq, size_t pos);
  Track* track(size_t iter) { return tracks.at(iter); }
  size_t n_tracks() { return tracks.size(); }
  uint8_t alphabet_size(size_t i) { return tracks[i]->alphabet_size(); }

private:
  vector<Track*> tracks;  // tracks which are used by emissions in this table
  vector<vector<double> > log_probs_;
  vector<vector<double> > original_log_probs_;  // i.e. before we rescaled the mute freq
};

}
#endif
