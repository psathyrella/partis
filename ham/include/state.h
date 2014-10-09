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
  void Parse(YAML::Node node, vector<string> state_names, Tracks trks);
  ~State();

  inline string &name() { return name_; }
  inline string &label() { return label_; }
  inline size_t index() { return index_; }  // index of this state in the HMM model
  inline vector<Transition*> *transitions() { return transitions_; }
  inline bitset<STATE_MAX> *to_states() { return &to_states_; }
  inline bitset<STATE_MAX> *from_states() { return &from_states_; }
  inline Transition *transition(size_t iter) { return (*transitions_)[iter]; }
  inline Transition *end_trans() { return end_trans_; }
              
  double emission_logprob(Sequences &seqs, size_t pos);
  inline double transition_logprob(size_t to_state) { return (*transitions_)[to_state]->log_prob(); }
  double end_transition_logprob();
              
  // property-setters for use in model::finalize()
  inline void AddToState(State *st) { to_states_[st->index()] = 1; }  // set bit in <to_states_> corresponding to <st>
  inline void AddFromState(State *st) { from_states_[st->index()] = 1; }  // set bit in <from_states_> corresponding to <st>
  inline void SetIndex(size_t val){ index_ = val; }
  void ReorderTransitions(map<string,State*>& state_indices);

  void Print();
private:
  string name_;       /* State name */
  string label_;      /* State feature path label */
              
  vector<Transition*>* transitions_;
  Transition* end_trans_;
  Emission emission_;
  Emission pair_emission_;

  // hmm model-level information (assigned in model::finalize)
  size_t index_;  // index of state in HMM (set in model:finalize)
  bitset<STATE_MAX> to_states_;
  bitset<STATE_MAX> from_states_;
};

}
#endif
