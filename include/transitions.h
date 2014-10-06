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
class Transition {
  friend class model;
  friend class State;
public:
  Transition(string to_state, double prob);
  
  inline void set_name(string& txt) { to_state_name_ = txt; }
  inline void set_state(State* st) { to_state_ = st; }
  inline void setTransType(transType type) { transition_type = type; }
  inline string& to_state_name(){return to_state_name_;};
  inline State* to_state(){return to_state_;};
  inline transType getTransitionType(){return transition_type;};
  
  inline double score() { return log_trans; }
  
  void print();
private:    
  transType transition_type;  //0: standard  1:USER DISTRIBUTION  2:INTERNAL DISTRIBUTION 3:LEXICAL
  string to_state_name_;  //What state we are transitioning to (Fill out when parsing)
  State* to_state_;    // pointer to the to-state (Filled during HMM finalization)
  double log_trans; //Log of standard transition probability
  bool _parseStandard(string&,stringList&, valueType);
};

}
#endif
