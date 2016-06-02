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
  void Init(Track *track);
  void ReplaceLogProbs(vector<double> new_log_probs);  // replace values in <log_probs_> with thsoe in <new_log_probs>
  void UnReplaceLogProbs();  // revert to the original log probs
  ~LexicalTable();

  void SetLogProbs(vector<double> logprobs) { log_probs_ = logprobs; }
  vector<double> log_probs() { return log_probs_; }  // NOTE returns a *copy*

  inline double LogProb(uint8_t index) { return log_probs_[index]; }
  double LogProb(Sequence *seq, size_t pos);
  Track *track() { return track_; }

private:
  Track *track_;
  vector<double> log_probs_;
  vector<double> original_log_probs_;  // i.e. before we rescaled the mute freq
};

}
#endif
