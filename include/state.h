#ifndef HAM_STATE_H
#define HAM_STATE_H

#include <string>
#include <vector>
#include <set>
#include <stdint.h>
#include <stdlib.h>
#include <bitset>

#include "text.h"
#include "emm.h"
#include "transitions.h"
#include "yaml-cpp/yaml.h"

using namespace std;
namespace ham {
class Transition;
  
class State {
public:
  State();
  State(string&,stringList&,Tracks&); //!Create state from string
  void parse(YAML::Node node, vector<string> state_names, Tracks trks);
  ~State();

  // accessors
  inline size_t index(){return index_;}  // index of this state in the HMM model
  inline string& name(){return name_;}
  inline string& label(){return label_;}
  inline vector<Transition*>* getTransitions(){return transi;}
  inline bitset<STATE_MAX>* to_states(){return &to_states_;}
  inline bitset<STATE_MAX>* from_states(){return &from_states_;}
  inline Transition* getTrans(size_t iter) { return (*transi)[iter]; }
  inline Transition* getEnding(){return endi;}
              
  inline double emission_score(Sequences& seqs, size_t pos) {  // seqs better have 1 or 2 seqs, 'cause I ain't checkin it here
    if (seqs.n_seqs() == 2)
      return pair_emission_.score(seqs, pos);
    else
      return emission_.score(seqs[0], pos);
  }
  inline double transition_score(size_t to_state) { return (*transi)[to_state]->score(); }
  double getEndTrans();
              
  void print();
              
  inline void setEndingTransition(Transition* trans){endi=trans;};
  inline void setName(string& txt){name_=txt;};
  inline void setLabel(string& txt){label_=txt;};
  inline void addToState(State* st){to_states_[st->index()]=1;};
  inline void addFromState(State* st){from_states_[st->index()]=1;};
  inline void setIter(size_t val){index_=val;};
              
  void _finalizeTransitions(map<string,State*>& state_index);

private:
  string name_;       /* State name */
  string label_;      /* State feature path label */
              
  vector<Transition*>* transi;
  Transition* endi;
  emm emission_;
  emm pair_emission_;

  // hmm model-level information (assigned in model::finalize)
  size_t index_;  //index of state in HMM (set in model:finalize)
  bitset<STATE_MAX> to_states_;
  bitset<STATE_MAX> from_states_;
              
  bool _parseHeader(string&);
  bool _parseTransition(string&,stringList&, Tracks&);
  bool _parseEmissions(string&,stringList&, Tracks&);
};
}
#endif
