#include "state.h"
namespace stochhmm {

// ----------------------------------------------------------------------------------------
state::state():endi(NULL), stateIterator(SIZE_MAX) {
  transi = new (nothrow) vector<transition*>;
}

// ----------------------------------------------------------------------------------------
state::state(string& txt, stringList& names,tracks& trcks) : endi(NULL), stateIterator(SIZE_MAX) {
  transi = new vector<transition*>;
  parse(txt,names,trcks);
}
  
// ----------------------------------------------------------------------------------------
state::~state(){
  delete transi;
  transi=NULL;
}

// ----------------------------------------------------------------------------------------
bool state::parse(string& txt, stringList& names,tracks& trks) {
  size_t stateHeaderInfo = txt.find("STATE:");
  size_t transitionsInfo = txt.find("TRANSITION:");
  size_t emissionInfo    = txt.find("EMISSION:");  // location of start of first emission definition
  stringList lines;
      
  //Extract and Parse Header Information
  if (stateHeaderInfo==string::npos) {
    cerr << "Couldn't identify state Header Information. Please check the formatting.  The supplied text was : " << txt << endl;
    return false;
  }
  string header = txt.substr(stateHeaderInfo,transitionsInfo-stateHeaderInfo);
  assert(_parseHeader(header));
              
  //Extract and Parse Transition Information
  if (transitionsInfo == string::npos) {
    cerr << "Couldn't identify state Transitions Information.  Please check the formatting.  The supplied text was : " << txt << endl;
    return false;
  }
  string trans = (emissionInfo==string::npos) ? txt.substr(transitionsInfo) : txt.substr(transitionsInfo, emissionInfo - transitionsInfo);
  assert(_parseTransition(trans, names, trks));
      
  //Check emission's existence  (only INIT state can have no emission)
  if (emissionInfo != string::npos) {
    string emmis = txt.substr(emissionInfo);
    assert(_parseEmissions(emmis, names, trks));
  } else {
    assert(name == "INIT");
  }
      
  return true;
}
  
// ----------------------------------------------------------------------------------------
bool state::_parseHeader(string& txt){
  clear_whitespace(txt,"\t ");
  stringList lst = splitString(txt,"\n:");
      
  size_t idx;
  // parse NAME
  if (lst.contains("NAME")) {
    idx = lst.indexOf("NAME");
    if (idx+1 < lst.size()) {
	idx++;
	name = lst[idx];
    } else {
      cerr << "The states NAME couldn't be parsed from " << txt << endl;
      return false;
    }
          
    if (name == "INIT")
      return true;
  } else {
    cerr << "No NAME found in state header information. Please check formatting. Following is header information recieved: " << txt << endl;
    return false;
  }
      
  // parse PATH_LABEL
  if (lst.contains("PATH")) {
    idx = lst.indexOf("PATH");
    if (idx+1 < lst.size()) {
      idx++;
      label = lst[idx];
    } else {
      cerr <<  "The states PATH_LABEL couldn't be parsed from " << txt << endl;
      return false;
    }
  } else {
    cerr << "No PATH_LABEL found in state header information. Please check formatting. Following is header information recieved: " << txt << endl;
    return false;
  }
      
  return true;
}
  
// ----------------------------------------------------------------------------------------
bool state::_parseTransition(string& txt, stringList& names, tracks& trks) {
  stringList lst;
  lst.splitND(txt,"TRANSITION:");
      
  for (size_t iter=0; iter<lst.size(); iter++) {
    stringList line(splitString(lst[iter],"\n"));
    line.removeLWS("\t ");
    line.removeComments();
          
    clear_whitespace(line[0],"\t");
    stringList head(splitString(line[0],":"));
    assert(head.contains("STANDARD"));
    transType tp(STANDARD);
                      
    valueType valtyp(UNKNOWN);
    if (head.contains("P(X)"))
      valtyp = PROBABILITY;
    else if (head.contains("LOG"))
      valtyp = LOG_PROB;
    else if (head.contains("COUNTS"))
      valtyp = COUNTS;
    else
      assert(0);
                              
    for(size_t iter=1; iter<line.size(); iter++) {
      transition* tmp_trans = new transition(tp);
      assert(tmp_trans->parse(line[iter], names, valtyp, trks));
      if (tmp_trans->getName() == "END")
	endi = tmp_trans;
      else
	addTransition(tmp_trans);
    }
  }
  return true;
}

// ----------------------------------------------------------------------------------------
bool state::_parseEmissions(string& txt, stringList& names, tracks& trks) {
  stringList lst;
  lst.splitND(txt,"EMISSION:");
  assert(lst.size() == 1 || lst.size() == 2);
  emission_.parse(lst[0], trks);
  if (lst.size() == 2) {
    assert(!emission_.pair());  // if pair and non-pair emission present, require the latter to be first
    pair_emission_.parse(lst[1], trks);
  }
  return true;
}
  
// ----------------------------------------------------------------------------------------
void state::print(){
  cout<< stringify() << endl;
}
  
// ----------------------------------------------------------------------------------------
string state::stringify(){
  string stateString;
  stateString+="STATE:\n";
  stateString+="\tNAME:\t" + name + "\n";
      
  if (name.compare("INIT")!=0){
    stateString+="\tPATH_LABEL:\t" + label + "\n";
  }
      
  string standardString;
  string distribString;
  string lexicalString;
  string pdfString;
      
  //Get Transitions Standard, Distribution, Lexical
  for(size_t i=0;i<transi->size();i++){
    if ((*transi)[i]==NULL){ continue;}
          
    transType tp = (*transi)[i]->getTransitionType();
    assert(tp == STANDARD);
    if (standardString.empty()){
      standardString+="TRANSITION:\tSTANDARD:\tLOG\n";
    }
    standardString+=(*transi)[i]->stringify() + "\n";
  }
      
      
  //Process Transition to Ending;
  if (endi!=NULL){
    if (standardString.empty()){
	standardString+="TRANSITIONS:\tSTANDARD:\tLOG\n";
    }
    standardString+=endi->stringify();
  }
      
      
  if (name.compare("INIT")==0){
    stateString+=standardString + "\n\n";
    return stateString;
  }
  else{
    if (!standardString.empty()){
	stateString+=standardString + "\n";
    }
  }
      
  //Print Emissions
  stateString += emission_.stringify();
      
  return stateString;
}
      
// // ----------------------------------------------------------------------------------------
// //! Get the emission value from a state at a particular position in the sequence
// //! \param seqs Sequences the model is analyzing
// //! \param iter Position in the sequences
// double state::emission_score(sequences &seqs, size_t iter) {
//   // double value(-INFINITY);
//   // for (size_t iseq=0; iseq<seqs.n_seqs(); ++iseq) {
//   //   double tmp_val = emission_.score(seqs[iseq], iter);  // get_emission should really be called get_emission_log_prob
//   //   if (iseq==0)
//   // 	value = tmp_val;
//   //   else
//   // 	value += tmp_val;
//   // }

//   // // TODO make these less hacky
//   // if (seqs.n_seqs()==2 && seqs[0][iter]!=seqs[1][iter])  // for pair hmm, multiply by the total mute prob if the nukes in the two seqs are different. NOTE this is a pretty approximate way to do this
//   //   value += log(mute_prob_);  // NOTE mute prob should only be non-1.0 for insert states
//   // if (seqs.size()==2 && seqs[0][iter]==seqs[1][iter] && seqs[0].get_undigitized(iter)!=getGFF()[0] && seqs[1].get_undigitized(iter)!=getGFF()[0])  // for *non*-insert states, if *both* seqs are mutated, but they're mutated to the *same* nucleotide, only apply *one* power of the mutation prob. we approximate this number by just taking the square root of the two. TODO don't need both of the last two clauses
//   //   value /= 2;

//   // return value;
// }

// // ----------------------------------------------------------------------------------------
// double state::get_emission_prob(sequence &seq1, sequence &seq2, size_t iter) {
//   return emissions[0]->get_emission(seq1, iter) + emissions[0]->get_emission(seq2, iter);
// }

// // ----------------------------------------------------------------------------------------
// double state::get_emission_prob(sequence &seq, size_t iter) {
//   return emissions[0]->get_emission(seq, iter);
// }
  
// ----------------------------------------------------------------------------------------
//! Get the transition value from a state at a particular position in the sequence
//! \param seqs Sequences the model is analyzing
//! \param to State that transition is being calculated to
//! \param iter Position in the sequence
// double state::transition_score(sequences &seqs, size_t to, size_t iter){
//   double value;
      
//   if ((*transi)[to]==NULL){
//     return -INFINITY;
//   }
//   else if ((*transi)[to]->transition_type==STANDARD){
//     value = (*transi)[to]->log_trans;
//   }
//   else{
//     // cerr << "Need to implement this functionality" <<endl;
//     assert(0);
//     // value = (*transi)[to]->getTransition(iter,&seqs);
          
//   }
//   return value;
// }
  
//! Get the log probability transitioning to end from the state
double state::getEndTrans(){
  if (endi==NULL){
    return -INFINITY;
  }
  return endi->log_trans;
}
  
//TODO: complete the checkLabels function
//! Checks the label tags for traceback and combine identifiers
void state::checkLabels(set<string>& labels, set<string>& gff, set<string>& name){
  return;
}
  
  
// ----------------------------------------------------------------------------------------
/* On initial import of the states they are pushed on the transi vector in
   the order written in model.   However, the analysis requires that they be
   in the particular position defined by state iterator.
   
   This function puts the transitions in the proper order for analysis
*/
void state::_finalizeTransitions(map<string,state*>& state_index){
              
  //Get size # of states, but correct by -1 because
  //initial state will be kept separate.
  size_t number_of_states = state_index.size();
  vector<transition*>* fixed_trans = new vector<transition*>(number_of_states-1,NULL);
      
  //Find the proper place for the transition and put it in the correct position
  for(size_t i = 0; i < transi->size(); i++){
    transition* temp = (*transi)[i];
    string name = temp->getName();
    state* st = state_index[name];
    if (st == NULL){
	cerr << "State: " << name << " was declared but not defined in the model." << endl;
	exit(2);
    }
    size_t index = st->getIterator();
    (*fixed_trans)[index]=temp;
    (*transi)[i]=NULL;
  }
      
  delete transi;  //Don't need the old transition vector anymore
  transi = fixed_trans;
  return;
}
}
