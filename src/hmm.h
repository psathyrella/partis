#ifndef HMM_H
#define HMM_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
#include <list>
#include <set>
#include <stdlib.h>
#include <bitset>

#include "state.h"
#include "track.h"
#include "text.h"
#include "stochTypes.h"

using namespace std;
namespace StochHMM{
/*! \class model
  Hidden Markov Model Class
  model class combines the States, and model information together in a single unit.
  This includes the states(emissions, transitions), initial and ending states, track 
  information(alphabet and ambiguous character definitions).
	 
  Provides functions to import the model from a text file
	 
  Model is used by trellis class to evaluates sequences.  
*/
class model {
public:
  model();
  void finalize();

  inline string& getName() { return name; }
  inline size_t n_states() { return states.size(); }
  inline string& getStateName(size_t iter){ return states[iter]->getName(); }
  inline string& getStateLabel(size_t iter){ return states[iter]->getLabel(); }
  state* getState(const string&);
  inline state*  getState(size_t iter) { return states[iter]; }
  inline state* operator[](size_t iter) { return states[iter]; }
  inline bitset<STATE_MAX>* getInitialTo() { return &(initial->to); }  //!Get vector of states that the initial state transitions to
  inline state*  getInitial() { return initial; }  //!Get pointer to the initial state
  inline state*  getEnding() { return ending; }  //!Get pointer to the ending state
  inline track* trck() { return tracks_[0]; }  // this is why we don't name classes with lowercase letters
  inline double overall_gene_prob() { return overall_gene_prob_; }

  void print();
  string stringify();

  bool parse(string);
		
  void addState(state*);
  inline void setInit(state* st) { initial=st; }
  inline void setEnd(state* st) { ending=st; }
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

  vector<state*> states; //!  All the states contained in the model

  map<string,state*> stateByName; //Ptr to state stored by State name;
		
  state* initial; //!Initial state q0
  state* ending;	//!Ending state
		
  bool _parseHeader(string&);	//! Function to parse header of the model from text file
  bool _parseTracks(string&);	//! Parse Tracks definitions from text file
		
  bool _parseStates(string&); //!Parse state from text file
  bool _splitStates(string&,stringList&); //!Split the state definitions into individual states from text file
  bool _getOrderedStateNames(stringList&,stringList&); //! Gets list of states names from model
		
  string _stringifyHeader();	//!Converts Header information from model to string representation found in text file
  string _stringifyTracks(); //!Converts Tracks information from model to string representation found in text file
  string _stringifyAmbig();  //!Converts Ambiguous Character information from model to text string
  string _stringifyScaling();//!Converts Scaling definitions from model to text string
  string _stringifyStates(); //!Converts States definitions from model to text string

  void _addStateToFromTransition(state*);
  void _checkTopology(state* st, vector<uint16_t>& visited); //!Checks to see that all states are connected and there
};
}
#endif
