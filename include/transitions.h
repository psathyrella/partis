#ifndef STOCHHMM_TRANSITIONS_H
#define STOCHHMM_TRANSITIONS_H

#include <vector> 
#include <string>
#include "text.h"
#include "track.h"
#include "mathutils.h"
#include "stochTypes.h"
#include "sequences.h"
#include <stdlib.h>

using namespace std;

namespace stochhmm {
class state;
class transition{
  friend class model;
  friend class state;
public:
  transition();
  transition(transType);
  transition(transType,valueType,bool);
  bool parse(string&, stringList&, valueType valtyp, tracks&);
  
  inline void setName(string& txt) { stateName = txt; }
  inline void setState(state* st) { toState = st; }
  inline void setTransType(transType type) { transition_type = type; }
  inline void setTransProb(double value){log_trans=value;};
  inline string& getName(){return stateName;};
  inline state* getState(){return toState;};
  inline transType getTransitionType(){return transition_type;};
  
  inline double score() { return log_trans; }
  
  inline bool isSimple() { return true; }
  inline bool isComplex() { return false; }
  
  void print();
  string stringify();
private:    
  transType transition_type;  //0: standard  1:USER DISTRIBUTION  2:INTERNAL DISTRIBUTION 3:LEXICAL
  string stateName;  //What state we are transitioning to (Fill out when parsing)
  state* toState;    //pointer to the state (Filled during HMM finalization)
  double log_trans; //Log of standard transition probability
  bool _parseStandard(string&,stringList&, valueType);
};
}
#endif
