#ifndef STATE_H
#define STATE_H

#include <string>
#include <vector>
#include <set>
#include "text.h"
#include "emm.h"
#include "transitions.h"
#include <stdint.h>
#include <stdlib.h>
#include <bitset>

#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

//Define size of bitset.  If not set in makefile then will use default 1024.
#ifndef STATE_MAX
#define STATE_MAX 1024
#endif

using namespace std;
namespace StochHMM {
class transition;
  
class state {
  friend class model;
public:
  state();
  state(string&,stringList&,tracks&,weights*, StateFuncs*); //!Create state from string
  bool parse(string&,stringList&,tracks&,weights*,StateFuncs*);
  ~state();

  // accessors
  inline size_t getIterator(){return stateIterator;};
  inline string& getName(){return name;};
  inline string& getGFF(){return gff;};
  inline string& getLabel(){return label;};
  inline vector<transition*>* getTransitions(){return transi;};
  inline bitset<STATE_MAX>* getTo(){return &to;};
  inline bitset<STATE_MAX>* getFrom(){return &from;};
  inline transition* getTrans(size_t iter) { return (*transi)[iter]; }
  inline transition* getEnding(){return endi;};
              
  inline double emission_score(sequences& seqs, size_t iter) {
    if (seqs.n_seqs() == 2)
      return emission_.score(seqs[0], iter) + emission_.score(seqs[1], iter);
    else
      return emission_.score(seqs[0], iter);
  }
  inline double transition_score(size_t to_state) { return (*transi)[to_state]->score(); }
  double getEndTrans();
              
  void print();
  string stringify();
              
  inline void addTransition(transition* trans){transi->push_back(trans);};
  inline void setEndingTransition(transition* trans){endi=trans;};
  inline void setName(string& txt){name=txt;};
  inline void setGFF(string& txt){gff=txt;};
  inline void setLabel(string& txt){label=txt;};
  inline void addToState(state* st){to[st->getIterator()]=1;};
  inline void addFromState(state* st){from[st->getIterator()]=1;};
  inline void setIter(size_t val){stateIterator=val;};
              
  void checkLabels(set<string>& ,set<string>& ,set<string>& );
  void _finalizeTransitions(map<string,state*>& state_index);
  double mute_prob_;  // total/mean mutation probability at this base. A rough approximation -- it's main purpose is to allow multiplying by a small number when two query seqs differ in a hypothesized insert region

private:
  string name;       /* State name */
  string gff ;       /* State features description */
  string label;      /* State feature path label */
              
  vector<transition*>* transi;
  transition* endi;
  emm emission_;

  //Linking State Information (These are assigned at model finalization)
  size_t stateIterator;  //index of state in HMM
  bitset<STATE_MAX> to;
  bitset<STATE_MAX> from;
              
  bool _parseHeader(string&);
  bool _parseTransition(string&,stringList&, tracks&, weights* , StateFuncs*);
  bool _parseEmissions(string&,stringList&, tracks&, weights*, StateFuncs*);
};
}
#endif
