#include "hmm.h"

namespace stochhmm {
// ----------------------------------------------------------------------------------------
model::model() : overall_gene_prob_(0),finalized(false),initial(NULL) {
  ending = new state;
}

// ----------------------------------------------------------------------------------------
//!Parses text model file
//!Splits the model into sections that are then parsed by the individiual classes
//!parse() functions.
bool model::parse(string infname) {
  YAML::Node config = YAML::LoadFile("test.yaml");
  YAML::Node tracks(config["tracks"]);
  for (size_t i=0; i<tracks.size(); ++i) {
    YAML::Node track(tracks[i]);
    for (YAML::const_iterator it=track.begin(); it!=track.end(); ++it) {
      cout << it->first << endl;
      for (size_t ic=0; ic<it->second.size(); ++ic) {
	cout << "  " << it->second[ic] << endl;
      }
    }
  }
  assert(0);
  
  string file_str = slurpFile(infname);  // string containing entire file
      
  size_t iheader = file_str.find("MODEL INFORMATION");
  size_t itrack = file_str.find("TRACK SYMBOL DEFINITIONS");
  size_t istate = file_str.find("STATE DEFINITIONS");
  size_t iblank;
  size_t nlChar;
      
  //Parse Model Informaton (Optional)
  if (iheader!=string::npos) {
    iblank = file_str.find("\n\n", iheader); //Get coordinates of splitting Model Information
          
    size_t nlCharEq = file_str.rfind("####\n",iblank);
    size_t nlCharNum= file_str.rfind("====\n",iblank);
    //Check for optional dividing line
    if (nlCharEq != string::npos) {
      nlChar = nlCharEq + 5;
    } else if (nlCharNum != string::npos) {
      nlChar = nlCharNum + 5;
    } else {  //No divider line
      nlChar = file_str.find("\n", iheader);
      nlChar++;
    }
          
    string head = file_str.substr(nlChar, iblank - nlChar);
    if (!_parseHeader(head))
      return false;
  }
      
  //Parse Tracks (Required)
  if (itrack!=string::npos) {
    iblank = file_str.find("\n\n", itrack);
          
    size_t nlCharEq = file_str.rfind("####\n", iblank);
    size_t nlCharNum = file_str.rfind("====\n", iblank);
    //Check for optional dividing line
    if (nlCharEq != string::npos) {
      nlChar = nlCharEq + 5;
    } else if (nlCharNum != string::npos) {
      nlChar = nlCharNum + 5;
    } else {  //No divider line
      nlChar = file_str.find("\n", itrack);
      nlChar++;
    }
          
    string trck(file_str.substr(nlChar, iblank-nlChar));
    if (!_parseTracks(trck))
      return false;
  } else {
    cerr << "Required section: TRACK SYMBOL DEFINITIONS missing from the model" << endl;
    return false;
  }
              
  //Parse State Definitions (Required)
  if (istate!=string::npos) {
    size_t blankNum = file_str.find("####\n", istate);
    size_t blankEq  = file_str.find("====\n", istate);
          
    iblank = file_str.find("####\n", istate);
          
    if (blankEq != string::npos) {
      iblank = blankEq + 5;
    } else if (blankNum != string::npos) {
      iblank = blankNum + 5;
    } else {
      iblank = file_str.find("\n", istate);
      iblank++;
    }
          
    nlChar = file_str.find("\n//END");
    if (nlChar == string::npos) {
      nlChar = file_str.size() - 1;
    }
          
    string stateTxt = file_str.substr(iblank, nlChar-iblank);
    if (!_parseStates(stateTxt))
      return false;
  } else {
    cerr << "Required sections <STATE DEFINITIONS> missing from the model" << endl;
    return false;
  }
      
  return true;
}
      
// ----------------------------------------------------------------------------------------
bool model::_parseHeader(string& txt) {
  stringList lst;
  size_t index;
  bool first(false);
  bool second(false);
  string headers[] = {"NAME", "DESCRIPTION", "CREATION_DATE","CREATION_COMMAND", "AUTHOR","NUM_ATTRIB","UPPER","LOWER"};
  string* head[] = {&name, &desc, &date, &command, &author};
      
  lst.fromTxt(txt);
  for(int i=0; i<5; i++) {
    if (lst.contains(headers[i])) {
      index = lst.indexOf(headers[i]);
      if (index+1 < lst.size()) {
        index++;
        (*head[i]) = lst[index];
	if (headers[i] == "DESCRIPTION") {
	  stringstream ss(lst[index]);
	  ss >> overall_gene_prob_;
	  assert(overall_gene_prob_ > 0.0 && overall_gene_prob_ < 1.0);
	}
      } else {
        cerr << "Couldn't parse " << headers[i] << " from \"MODEL INFORMATION\" section." << endl;
        return false;
      }
    }
  }
      
  //Set Numerical Attributes of Model
  if (lst.contains(headers[5])){
    index = lst.indexOf(headers[5]);
    if (index+1 < lst.size()) {
      index++;
      double tempValue;
      if (!stringToDouble(lst[index], tempValue)) {
        cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << endl;
        return false;
      }
    } else {
      cerr << "Couldn't parse " << headers[5] << " value from \"MODEL INFORMATION\" section." << endl;
      return false;
    }
  } else if (lst.contains(headers[6]) && lst.contains(headers[7])) {
    index = lst.indexOf(headers[6]);
    if (index+1<lst.size()) {
      index++;
      double tempValue;
      if (!stringToDouble(lst[index], tempValue)) {
        cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << endl;
        return false;
      }
      second=true;
    } else {
      cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << endl;
      return false;
    }
          
    index = lst.indexOf(headers[7]);
    if (index+1<lst.size()) {
      index++;
      double tempValue;
      if (!stringToDouble(lst[index], tempValue)) {
        cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << endl;
        return false;
      }
      first=true;
    } else {
      cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << endl;
      return false;
    }
          
    if (first && second) {
      assert(0);
    } else {
      cerr << "Unable to parse both UPPER and LOWER" << endl;
      return false;
    }
  } else if (lst.contains(headers[6])) {
    index = lst.indexOf(headers[6]);
    if (index+1<lst.size()) {
      index++;
      double tempValue;
      if (!stringToDouble(lst[index], tempValue)) {
        cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << endl;
        return false;
      }
    } else {
      cerr << "Couldn't parse " << headers[6] << " value from \"MODEL INFORMATION\" section." << endl;
      return false;
    }
  } else if (lst.contains(headers[7])) {
    index = lst.indexOf(headers[6]);
    if (index+1<lst.size()) {
      index++;
      double tempValue;
      if (!stringToDouble(lst[index], tempValue)) {
        cerr << "Numerical attribute couldn't be converted to numerical value: " << lst[index] << endl;
        return false;
      }
    } else {
      cerr << "Couldn't parse " << headers[7] << " value from \"MODEL INFORMATION\" section." << endl;
      return false;
    }
  }
      
  return true;
}
  
// ----------------------------------------------------------------------------------------
bool model::_parseTracks(string& txt) {
  stringList lst;
  lst.splitString(txt, "\n");
  assert(lst.size() == 1);
  track *tmp_trk = new track;
  tmp_trk->parse(lst[0]);
  tracks_.push_back(tmp_trk);
  return true;
}

// ----------------------------------------------------------------------------------------
bool model::_parseStates(string& txt) {
  //1. split sections and identify any template sections
  //2. get state names list
  //3. create and parse states
      
  stringList stats;
  _splitStates(txt,stats);
      
      
  stringList NameList;
  _getOrderedStateNames(stats,NameList);
      
  for(size_t iter=0;iter<stats.size();iter++){
    state* st = new(nothrow) state;
          
    if (st==NULL){
      cerr << "OUT OF MEMORY\nFile" << __FILE__ << "Line:\t"<< __LINE__ << endl;
      exit(1);
    }
          
    if (!st->parse(stats[iter], NameList, tracks_)){
      delete st;
      return false;
    }
          
          
    if (st->getName() == "INIT"){
      initial=st;
      stateByName[st->getName()]=st;
    }
    else{
	assert(states.size() < STATE_MAX);
      states.push_back(st);
      stateByName[st->getName()]=st;
    }
  }
      
      
  //Post process states to create final state with only transitions from filled out.
  finalize();
      
      
  return true;
}
  
  
// ----------------------------------------------------------------------------------------
void model::addState(state* st){
  assert(states.size() < STATE_MAX);
  states.push_back(st);
  stateByName[st->getName()]=st;
  return;
};
      
  
// ----------------------------------------------------------------------------------------
//Split states into individual state strings
bool model::_splitStates(string& txt ,stringList& sts){
      
  size_t start=0;
  size_t end=0;
      
  while(start!=string::npos){
    end=txt.find("STATE:",start+1);
          
    string st = txt.substr(start,end-start);
          
    clear_whitespace(st, "#");
          
    if (st.find("TEMPLATE:")!=string::npos){
              
      stringList tmpd;
              
      for(size_t i=0;i<tmpd.size();i++){
        sts.push_back(tmpd[i]);
      }
    }
    else{
      sts.push_back(st);
    }
                      
    start=txt.find("STATE:",end);
          
  }
  return true;
}

// ----------------------------------------------------------------------------------------
bool model::_getOrderedStateNames(stringList& states, stringList& names){
  for(size_t i=0;i<states.size();i++){
    size_t nameHeader=states[i].find("NAME:");
    size_t nameLineEnding=states[i].find_first_of("\n",nameHeader);
    string name = states[i].substr(nameHeader+5,nameLineEnding-(nameHeader+5));
    clear_whitespace(name, " \t\n");
    if (names.containsExact(name)){
      cerr << "State with name of: " << name << " is defined twice in the model\n";
      return false;
    }
    else{
      names.push_back(name);
    }
  }
  return true;
}
  
  
//!Get pointer to state by state name
//!\param txt String name of state
//!\return pointer to state if it exists;
//!\return NULL if state doesn't exist in model
state* model::getState(const string& txt){
  if (stateByName.count(txt)){
    return stateByName[txt];
  }
  else {
    return NULL;
  }
}
      
  
//!Finalize the model before performing decoding
//!Sets transitions, checks labels, Determines if model is basic or requires intermittent tracebacks
void model::finalize(){
  if (!finalized){
    //Add States To Transitions
    set<string> labels;
    set<string> gff;
    set<string> name;
                      
    //Create temporary hash of states for layout
    for(size_t i=0;i<states.size();i++){
      labels.insert(states[i]->getLabel());
      // gff.insert(states[i]->getGFF());
      // gff.insert(states[i]->getName());
    }
                      
    for (size_t i=0; i < states.size() ; ++i){
      states[i]->setIter(i);
    }
                      
          
          
    //Add states To and From transition
          
    for(size_t i=0;i<states.size();i++){
      states[i]->checkLabels(labels,gff,name);
      _addStateToFromTransition(states[i]);
              
    }
                      
    _addStateToFromTransition(initial);
                      
                      
    //Now that we've seen all the states in the model
    //We need to fix the States transitions vector transi, so that the state
    //iterator correlates to the position within the vector
    for(size_t i=0;i<states.size();i++){
      states[i]->_finalizeTransitions(stateByName);
    }
    initial->_finalizeTransitions(stateByName);
          
    //Check to see if model is basic model
    //Meaning that the model doesn't call outside functions or perform
    //Tracebacks for explicit duration.
    //If explicit duration exist then we'll keep track of which states
    //they are in explicit_duration_states
    //            for(size_t i=0;i<states.size();i++){
    //                vector<transition*>* transitions = states[i]->getTransitions();
    //                for(size_t trans=0;trans<transitions->size();trans++){
    //                                      if ((*transitions)[trans] == NULL){
    //                                              continue;
    //                                      }
    //                                      
    //                    if ((*transitions)[trans]->FunctionDefined()){
    //                        basicModel=false;
    //                        break;
    //                    }
    //                    else if ((*transitions)[trans]->getTransitionType()!=STANDARD || (*transitions)[trans]->getTransitionType()!=LEXICAL){
    //                                              
    //                                              if ((*transitions)[trans]->getTransitionType() == DURATION){
    //                                                      (*explicit_duration_states)[i]=true;
    //                                              }
    //                                              
    //                        basicModel=false;
    //                        break;
    //                    }
    //                }
    //            }
                      
    checkTopology();
                      
                      
    // //Assign StateInfo
    // for(size_t i=0; i < states.size();i++){
    //   string& st_name = states[i]->getName();
    //   info.stateIterByName[st_name]=i;
    //   info.stateIterByLabel[st_name].push_back(i);
    //   info.stateIterByGff[st_name].push_back(i);
    // }
    finalized = true;
  }
  return;
}
      

// ----------------------------------------------------------------------------------------
void model::_addStateToFromTransition(state* st){
  vector<transition*>* trans;
      
  //Process Initial State
  trans = st->getTransitions();
  for(size_t i=0;i<trans->size();i++){
    state* temp;
    temp=this->getState((*trans)[i]->getName());
    if (temp){
      st->addToState(temp); //Also add the ptr to state vector::to
      (*trans)[i]->setState(temp);
      if (st!=initial){
        temp->addFromState(st);
      }
    }
  }
      
  if (st->endi){
    ending->addFromState(st);
  }
}
  
  
//!Get string representation of model
//!\return string
string model::stringify(){
  string model;
  string lnSep(50,'=');
  model+="#STOCHHMM MODEL FILE\n\n";
      
  model+=_stringifyHeader();
  model+=_stringifyTracks();
  model+=_stringifyStates();
  return model;
}
  
// ----------------------------------------------------------------------------------------
string model::_stringifyHeader(){
  string headerString;
  string lnSep(50,'=');
  headerString+="MODEL INFORMATION\n" + lnSep + "\n";
  string headers[] = {"MODEL_NAME", "MODEL_DESCRIPTION", "MODEL_CREATION_DATE","MODEL_CREATION_COMMAND", "MODEL_AUTHOR","MODEL_NUM_ATTRIB","MODEL_NUM_ATTRIB_UPPER","MODEL_NUM_ATTRIB_LOWER"};
  string* head[]={&name, &desc, &date, &command,&author};
      
  for(int i=0;i<5;i++){
    if (!(*head[i]).empty()){
      headerString+=headers[i] + ":\t" + (*head[i]) + "\n";
    }
  }
      
      
  headerString+="\n";
      
  return headerString;
}
  
// ----------------------------------------------------------------------------------------
string model::_stringifyTracks(){
  return tracks_[0]->stringify() + "\n";
}
  
// ----------------------------------------------------------------------------------------
string model::_stringifyStates(){
  string stateString;
  string lnSep(50,'#');
  stateString+="STATE DEFINITIONS\n" + lnSep + "\n";
      
  if (initial){
    stateString+=initial->stringify();
  }
      
  stateString+= lnSep + "\n";
      
  for(size_t i=0;i<states.size();i++){
    if (states[i]){
      stateString+=states[i]->stringify();
      stateString+=lnSep + "\n";
    }
  }
  stateString+= "//END";
  return stateString;
}
  
  
// ----------------------------------------------------------------------------------------
//!Print the string representation to stdout
void model::print(){
  cout << stringify() << endl;
  return;
}

// ----------------------------------------------------------------------------------------
bool model::checkTopology(){
  //!Check model topology
  //!Iterates through all states to check to see if there are any:
  //! 1. Orphaned States
  //! 2. Dead end States
  //! 3. Uncompleted States
              
  vector<bool> states_visited (states.size(),false);
  vector<uint16_t> visited;
              
  bool ending_defined(false);
              
  _checkTopology(initial, visited);
              
  while (visited.size()>0){
    uint16_t st_iter = visited.back();
    visited.pop_back();
                      
    if (!states_visited[st_iter]){
      vector<uint16_t> tmp_visited;
      _checkTopology(states[st_iter],tmp_visited);
      size_t num_visited = tmp_visited.size();
                              
      //Check orphaned
      if (num_visited == 0 ){
        //No transitions
        //cerr << "Warning: State: "  << states[st_iter]->getName() << " has no transitions defined\n";
      }
      else if (num_visited == 1 && tmp_visited[0] == st_iter){
        //Orphaned
        if(states[st_iter]->getEnding() == NULL){
          cerr << "State: "  << states[st_iter]->getName() << " is an orphaned state that has only transition to itself\n";
        }
        //                                      else{
        //                                              cerr << "State: "  << states[st_iter]->getName() << " may be an orphaned state that only has transitions to itself and END state.\n";
        //                                      }
      }
                              
      for(size_t i=0; i < tmp_visited.size(); i++){
        if (!states_visited[tmp_visited[i]]){
          visited.push_back(tmp_visited[i]);
        }
      }
                              
      states_visited[st_iter] = true;
    }
  }
              
  //Check for defined ending
  for(size_t i=0; i< states.size() ; i++){
    if ( states[i]->getEnding() != NULL){
      ending_defined = true;
      break;
    }
  }
              
  if (!ending_defined){
    cerr << "No END state defined in the model\n";
  }
              
  for(size_t i=0; i< states_visited.size(); i++){
    if (!states_visited[i]){
      cerr << "State: "  << states[i]->getName() << " doesn't have valid model topology\n\
                              Please check the model transitions\n";
      return false;
    }
  }
  return true;
}
      
// ----------------------------------------------------------------------------------------
void model::_checkTopology(state* st, vector<uint16_t>& visited){
  //Follow transitions to see if every state is visited
  for(size_t i = 0 ; i < st->transi->size() ; i++){
    if (st->transi->at(i) != NULL){
      visited.push_back(st->transi->at(i)->getState()->getIterator());
    }
                      
  }
  return;
}
}
