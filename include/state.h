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
class transition;
  
class State {
  friend class model;
public:
  State();
  State(string&,stringList&,Tracks&); //!Create state from string
  void parse(YAML::Node node, vector<string> state_names, Tracks trks);
  ~State();

  // accessors
  inline size_t getIterator(){return stateIterator;};
  inline string& getName(){return name;};
  inline string& getLabel(){return label;};
  inline vector<transition*>* getTransitions(){return transi;};
  inline bitset<STATE_MAX>* getTo(){return &to;};
  inline bitset<STATE_MAX>* getFrom(){return &from;};
  inline transition* getTrans(size_t iter) { return (*transi)[iter]; }
  inline transition* getEnding(){return endi;};
              
  inline double emission_score(Sequences& seqs, size_t pos) {  // seqs better have 1 or 2 seqs, 'cause I ain't checkin it here
    if (seqs.n_seqs() == 2)
      return pair_emission_.score(seqs, pos);
    else
      return emission_.score(seqs[0], pos);
  }
  inline double transition_score(size_t to_state) { return (*transi)[to_state]->score(); }
  double getEndTrans();
              
  void print();
              
  inline void setEndingTransition(transition* trans){endi=trans;};
  inline void setName(string& txt){name=txt;};
  inline void setLabel(string& txt){label=txt;};
  inline void addToState(State* st){to[st->getIterator()]=1;};
  inline void addFromState(State* st){from[st->getIterator()]=1;};
  inline void setIter(size_t val){stateIterator=val;};
              
  void _finalizeTransitions(map<string,State*>& state_index);

private:
  string name;       /* State name */
  string label;      /* State feature path label */
              
  vector<transition*>* transi;
  transition* endi;
  emm emission_;
  emm pair_emission_;

  //Linking State Information (These are assigned at model finalization)
  size_t stateIterator;  //index of state in HMM (set in model:finalize)
  bitset<STATE_MAX> to;
  bitset<STATE_MAX> from;
              
  bool _parseHeader(string&);
  bool _parseTransition(string&,stringList&, Tracks&);
  bool _parseEmissions(string&,stringList&, Tracks&);
};
}
#endif
