#ifndef HAM_MODEL_H
#define HAM_MODEL_H

#include <fstream>
#include "state.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {

class Model {
public:
  Model();
  void Parse(string);
  void AddState(State*);
  void Finalize();

  inline string &name() { return name_; }
  inline Track *track(size_t it) { return tracks_[it]; }
  inline size_t n_states() { return states_.size(); }
  inline State *state(string name) { assert(states_by_name_.count(name)); return states_by_name_[name]; }
  inline State *state(size_t ist) {
    if(ist >= states_.size())
      throw runtime_error("ERROR state index too large in model::state(): " + to_string(ist));
    return states_[ist];
  }
  inline bitset<STATE_MAX> *initial_to_states() { return initial_->to_states(); }  //!Get vector of states that the initial state transitions to
  inline State *init_state() { return initial_; }  //!Get pointer to the initial state
  inline double overall_prob() { return overall_prob_; }

private:
  string name_;
  double overall_prob_;  // overall (say, prior) probability of this hmm
  Tracks tracks_;
  vector<State*> states_; //!  All the states contained in the model
  map<string, State*> states_by_name_; //Ptr to state stored by State name;
  State* initial_;
  State* ending_;
  bool finalized_;

  void FinalizeState(State *st);
  void CheckTopology();
  void AddToStateIndices(State* st, vector<uint16_t>& visited); // that's 'to-state', as in, 'here we push back the to-state indices onto <visited>'
};

}
#endif
