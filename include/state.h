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
  State(string&,stringList&,Tracks&); //!Create state from string
  void Parse(YAML::Node node, vector<string> state_names, Tracks trks);
  ~State();

  // accessors
  inline string& name(){return name_;}
  inline string& label(){return label_;}
  inline size_t index(){return index_;}  // index of this state in the HMM model
  inline vector<Transition*>* transitions(){return transitions_;}
  inline bitset<STATE_MAX>* to_states(){return &to_states_;}
  inline bitset<STATE_MAX>* from_states(){return &from_states_;}
  inline Transition* getTrans(size_t iter) { return (*transitions_)[iter]; }
  inline Transition* end_trans(){return end_trans_;}
              
  inline double emission_score(Sequences& seqs, size_t pos) {  // seqs better have 1 or 2 seqs, 'cause I ain't checkin it here
    if (seqs.n_seqs() == 2)
      return pair_emission_.score(seqs, pos);
    else
      return emission_.score(seqs[0], pos);
  }
  inline double transition_score(size_t to_state) { return (*transitions_)[to_state]->score(); }
  double getEndTrans();
              
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
              
  bool _parseHeader(string&);
  bool _parseTransition(string&,stringList&, Tracks&);
  bool _parseEmissions(string&,stringList&, Tracks&);
};

}
#endif
