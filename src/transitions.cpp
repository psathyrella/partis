#include "transitions.h"
namespace stochhmm {
  
// ----------------------------------------------------------------------------------------
transition::transition(transType type) : transition_type(type) {
  log_trans = -INFINITY;        
  toState = NULL;
}
  
// ----------------------------------------------------------------------------------------
transition::transition(transType type, valueType valtyp, bool survival) : transition_type(type) {
  log_trans = -INFINITY;
  toState = NULL;
}

// ----------------------------------------------------------------------------------------
bool transition::parse(string& txt,stringList& names, valueType valtyp, tracks& trks) {
  assert(_parseStandard(txt,names, valtyp));
  return true;
}
  
// ----------------------------------------------------------------------------------------
bool transition::_parseStandard(string& txt, stringList& names, valueType valtyp) {
  stringList line;
  line.splitString(txt,"\t:");
  assert(line.size() >= 2);  // cerr << "Line should contain 2 values (STATE  VALUE).  Couldn't parse:\n" << txt << endl;
  stateName = line[0];
      
  if (!names.contains(stateName) && stateName!="END"){
    cerr << "Tried to create a transition in the model to state named : " << stateName << " but there is no state with that name.  Please check the formatting. \n";
    return false;
  }
      
  if (valtyp == PROBABILITY) {
    double tempValue;
    if (!stringToDouble(line[1], tempValue)) {
      cerr << "Probability value not numeric: " << line[1] << endl;
      return false;
    }
    log_trans = log(tempValue);
  } else if (valtyp == LOG_PROB) {
    double tempValue;
    if (!stringToDouble(line[1], tempValue)) {
      cerr << "Log Probability value not numeric: " << line[1] << endl;
      return false;
    }
    log_trans = tempValue;
  }
  return true;
}
  
// ----------------------------------------------------------------------------------------
void transition::print(){
  cout << stringify() << endl;
}
  
// ----------------------------------------------------------------------------------------
string transition::stringify() {
  string transString;
  transString+="\t" + stateName + ":\t" + double_to_string(log_trans);
  return transString;
}
}
