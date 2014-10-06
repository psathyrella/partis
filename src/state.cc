#include "state.h"
namespace ham {

// ----------------------------------------------------------------------------------------
State::State() : end_trans_(NULL), index_(SIZE_MAX) {
  transitions_ = new (nothrow) vector<Transition*>;
}
  
// ----------------------------------------------------------------------------------------
State::~State(){
  delete transitions_;
  transitions_=NULL;
}

// ----------------------------------------------------------------------------------------
void State::parse(YAML::Node node, vector<string> state_names, Tracks trks) {
  name_ = node["name"].as<string>();
  label_ = node["label"].as<string>();

  double total(0.0); // make sure things add to 1.0
  for (YAML::const_iterator it=node["transitions"].begin(); it!=node["transitions"].end(); ++it) {
    string to_state(it->first.as<string>());
    if (to_state != "end" && find(state_names.begin(), state_names.end(), to_state) == state_names.end()) {  // make sure transition is either to "end", or to a state that we know about
      cout << "ERROR attempted to add transition to unknown state \"" << to_state << "\"" << endl;
      assert(0);
    }
    double prob(it->second.as<double>());
    total += prob;
    Transition *trans = new Transition(to_state, prob);
    if (trans->to_state_name() == "end")
      end_trans_ = trans;
    else
      transitions_->push_back(trans);
  }
  // TODO use something cleverer than a random hard coded EPS
  assert(fabs(total-1.0) < EPS);  // make sure transition probs sum to 1.0

  if (name_ == "init")
    return;
      
  if (node["emissions"])
    emission_.parse(node["emissions"], "single", trks);
  if (node["pair_emissions"])
    pair_emission_.parse(node["pair_emissions"], "pair", trks);
}
  
// ----------------------------------------------------------------------------------------
void State::print() {
  cout << "state: " << name_ << " (" << label_ << ")" << endl;;

  cout << "  transitions:" << endl;;
  for(size_t i=0; i<transitions_->size(); ++i) {
    if ((*transitions_)[i]==NULL){ assert(0); continue;}  // wait wtf would this happen?
    (*transitions_)[i]->print();
  }

  if (end_trans_)
    end_trans_->print();
      
  if (name_ == "init")
    return;

  cout << "  emissions:" << endl;;
  emission_.print();  // TODO allow state to have only one or the other of these
  cout << "  pair emissions:" << endl;
  pair_emission_.print();
}
      
//! Get the log probability transitioning to end from the state.
double State::getEndTrans(){
  if (end_trans_==NULL){
    return -INFINITY;
  }
  return end_trans_->log_trans;
}
  
  
// ----------------------------------------------------------------------------------------
/* On initial import of the states they are pushed on the transi vector in
   the order written in model.   However, the analysis requires that they be
   in the particular position defined by state iterator.
   
   This function puts the transitions in the proper order for analysis
*/
void State::_finalizeTransitions(map<string,State*>& state_index){
              
  //Get size # of states, but correct by -1 because
  //initial state will be kept separate.
  size_t number_of_states = state_index.size();
  vector<Transition*>* fixed_trans = new vector<Transition*>(number_of_states-1,NULL);
      
  //Find the proper place for the transition and put it in the correct position
  for(size_t i = 0; i < transitions_->size(); i++){
    Transition* temp = (*transitions_)[i];
    string name = temp->to_state_name();
    State* st = state_index[name];
    if (st == NULL){
	cerr << "State: " << name << " was declared but not defined in the model." << endl;
	exit(2);
    }
    size_t index = st->index();
    (*fixed_trans)[index]=temp;
    (*transitions_)[i]=NULL;
  }
      
  delete transitions_;  //Don't need the old transition vector anymore
  transitions_ = fixed_trans;
  return;
}
}
