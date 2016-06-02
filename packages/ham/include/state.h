#ifndef HAM_STATE_H
#define HAM_STATE_H

#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include <stdlib.h>
#include <bitset>

#include "text.h"
#include "emission.h"
#include "transitions.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {
class Transition;

class State {
public:
  State();
  void Parse(YAML::Node node, vector<string> state_names, Track *track);
  void RescaleOverallMuteFreq(double factor);  // Rescale emissions by the ratio <factor>
  void UnRescaleOverallMuteFreq();  // undo the above
  ~State();

  inline string name() { return name_; }
  inline string abbreviation() { return name_.substr(0, 1); }
  inline size_t index() { return index_; }  // index of this state in the HMM model
  inline vector<Transition*> *transitions() { return transitions_; }
  inline bitset<STATE_MAX> *to_states() { return &to_states_; }
  inline bitset<STATE_MAX> *from_states() { return &from_states_; }
  inline vector<size_t> *from_state_indices() { return &from_state_indices_; }
  inline Transition *transition(size_t iter) { return (*transitions_)[iter]; }
  inline Transition *trans_to_end() { return trans_to_end_; }

  double EmissionLogprob(uint8_t ch);
  double EmissionLogprob(Sequences *seqs, size_t pos);
  inline double transition_logprob(size_t to_state) { return (*transitions_)[to_state]->log_prob(); }
  double end_transition_logprob();

  // property-setters for use in model::finalize()
  inline void AddToState(State *st) { to_states_[st->index()] = 1; }  // set bit in <to_states_> corresponding to <st>
  inline void AddFromState(State *st) { from_states_[st->index()] = 1; }  // set bit in <from_states_> corresponding to <st>
  inline void SetIndex(size_t val) { index_ = val; }
  void ReorderTransitions(map<string, State*>& state_indices);

  void SetFromStateIndices();

  void Print();
private:
  string name_, germline_nuc_;
  double ambiguous_emission_logprob_;
  string ambiguous_char_;
  vector<Transition*> *transitions_;
  Transition *trans_to_end_;
  Emission emission_;

  // hmm model-level information (assigned in model::finalize)
  size_t index_;  // position of this state in the vector model::states_ (set in model::finalize)
  bitset<STATE_MAX> to_states_;
  bitset<STATE_MAX> from_states_;
  vector<size_t> from_state_indices_;  // same information as <from_states_>, but hopefully faster to iterate over
};

}
#endif
