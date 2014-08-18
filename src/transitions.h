#ifndef TRANSITIONS_H
#define TRANSITIONS_H

#include <vector> 
#include <string>
#include "text.h"
#include "track.h"
#include "stochMath.h"
#include "stochTypes.h"
#include "weight.h"
#include "sequences.h"
#include <stdlib.h>

namespace StochHMM {
class state;
class transition{
  friend class model;
  friend class state;
public:
  transition();
  transition(transType);
  transition(transType,valueType,bool);
  bool parse(std::string&, stringList&, valueType valtyp, tracks& ,weights*, StateFuncs*);
  
  inline void setName(std::string& txt) { stateName = txt; }
  inline void setState(state* st) { toState = st; }
  inline void setTransType(transType type) { transition_type = type; }
  inline void setTransProb(double value){log_trans=value;};
  inline std::string& getName(){return stateName;};
  inline state* getState(){return toState;};
  inline transType getTransitionType(){return transition_type;};
  
  inline double score() { return log_trans; }
  
  inline bool isSimple() { return true; }
  inline bool isComplex() { return false; }
  
  void print();
  std::string stringify();
private:    
  transType transition_type;  //0: standard  1:USER DISTRIBUTION  2:INTERNAL DISTRIBUTION 3:LEXICAL
  std::string stateName;  //What state we are transitioning to (Fill out when parsing)
  state* toState;    //pointer to the state (Filled during HMM finalization)
  double log_trans; //Log of standard transition probability
  bool _parseStandard(std::string&,stringList&, valueType);
};
}
#endif
