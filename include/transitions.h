#ifndef HAM_TRANSITIONS_H
#define HAM_TRANSITIONS_H

#include <vector> 
#include <string>
#include "text.h"
#include "track.h"
#include "mathutils.h"
#include "types.h"
#include "sequences.h"
#include <stdlib.h>

using namespace std;

namespace ham {
class State;
class transition{
  friend class model;
  friend class State;
public:
  transition(string to_state, double prob);
  
  inline void setName(string& txt) { stateName = txt; }
  inline void setState(State* st) { toState = st; }
  inline void setTransType(transType type) { transition_type = type; }
  // inline void setTransProb(double value){log_trans=value;};
  inline string& getName(){return stateName;};
  inline State* getState(){return toState;};
  inline transType getTransitionType(){return transition_type;};
  
  inline double score() { return log_trans; }
  
  void print();
private:    
  transType transition_type;  //0: standard  1:USER DISTRIBUTION  2:INTERNAL DISTRIBUTION 3:LEXICAL
  string stateName;  //What state we are transitioning to (Fill out when parsing)
  State* toState;    //pointer to the state (Filled during HMM finalization)
  double log_trans; //Log of standard transition probability
  bool _parseStandard(string&,stringList&, valueType);
};
}
#endif
