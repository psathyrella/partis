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
  void finalize();

  inline string& getName() { return name; }
  inline size_t n_states() { return states.size(); }
  inline string& getStateName(size_t iter){ return states[iter]->getName(); }
  inline string& getStateLabel(size_t iter){ return states[iter]->getLabel(); }
  State* getState(const string&);
  inline State*  getState(size_t iter) { return states[iter]; }
  // inline State* operator[](size_t iter) { return states[iter]; }
  inline bitset<STATE_MAX>* getInitialTo() { return &(initial->to); }  //!Get vector of states that the initial state transitions to
  inline State*  getInitial() { return initial; }  //!Get pointer to the initial state
  inline State*  getEnding() { return ending; }  //!Get pointer to the ending state
  inline track* trck() { return tracks_[0]; }  // this is why we don't name classes with lowercase letters
  inline double overall_gene_prob() { return overall_gene_prob_; }

  void addState(State*);
  inline void setInit(State* st) { initial=st; }
  inline void setEnd(State* st) { ending=st; }
  bool checkTopology();
		
private:
  double overall_gene_prob_;  // prob to select this gene

  //!Flag set to tell whether the transitions bitsets have been set foreach
  //!state.  Model is also checked for correct order of states
  bool finalized;   
		
  string name;	//! Model Name
  string desc;	//! Model Description
  string date;   //! Model Creation Date
  string command; //! Model Creation Command
  string author;  //! Model Author
		
  tracks tracks_;

  vector<State*> states; //!  All the states contained in the model

  map<string,State*> stateByName; //Ptr to state stored by State name;
		
  State* initial; //!Initial state q0
  State* ending;	//!Ending state
		
  void _addStateToFromTransition(State*);
  void _checkTopology(State* st, vector<uint16_t>& visited); //!Checks to see that all states are connected and there
};

}
#endif
