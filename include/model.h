#ifndef HAM_HMM_H
#define HAM_HMM_H

#include "state.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {

class model {
public:
  model();
  void parse(string);
  void add_state(State*);
  inline void set_init(State* st) { initial_=st; }
  inline void set_end(State* st) { ending_=st; }
  void finalize();
		
  inline string& name() { return name_; }
  inline size_t n_states() { return states_.size(); }
  inline State* state(string name) { assert(states_by_name_.count(name)); return states_by_name_[name]; }
  inline State* state(size_t ist) { assert(ist < states_.size()); return states_[ist]; }
  inline bitset<STATE_MAX>* initial_to_states() { return initial_->to_states(); }  //!Get vector of states that the initial state transitions to
  inline State* init_state() { return initial_; }  //!Get pointer to the initial state
  inline State* end_state() { return ending_; }  //!Get pointer to the ending state
  inline double overall_prob() { return overall_prob_; }

private:
  string name_;
  double overall_prob_;  // overall (say, prior) probability of this hmm
  Tracks tracks_;
  vector<State*> states_; //!  All the states contained in the model
  map<string,State*> states_by_name_; //Ptr to state stored by State name;
  State* initial_; //!Initial state q0
  State* ending_;	//!Ending state
  bool finalized_;
		
  void _addStateToFromTransition(State*);
  bool checkTopology();
  void _checkTopology(State* st, vector<uint16_t>& visited); //!Checks to see that all states are connected and there
};

}
#endif
