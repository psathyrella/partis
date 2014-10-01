#include "tracebackpath.h"

namespace stochhmm {


  int WID=80;
  //!Create a TracebackPath
  //!\param modl Pointer to model file 
  TracebackPath::TracebackPath(model* modl){
    hmm=modl;
  }

  //!Pushes a state index onto the end of the path
  //!\param state Index of state to add
  void TracebackPath::push_back(int state){
    trace_path.push_back(state);
  }
    
  //!Returns the size (ie. length) of the TracebackPath
  size_t TracebackPath::size()const {
    return trace_path.size();
  }
    
  //!Clears all traceback path information
  void TracebackPath::clear(){
    trace_path.clear();
  }
    
    
  //! Get the path to vector<int>
  //! \param [out] pth vector<int> that represents path
  void TracebackPath::path(vector<int>& pth){
    pth=trace_path;
    return ;
  }

  //TODO: change assignment to lhs
  //! Get the label of the TracebackPath and assigns to vector<string> ref
  void TracebackPath::label(vector<string>& pth){
    for(size_t k=trace_path.size()-1; k!=SIZE_MAX; k--){
      state* st = hmm->getState(trace_path[k]);
      pth.push_back(st->getLabel());
    }
    return;
  }
    
  //!Get string of path label traceback
  //! \param[out] pth string
  void TracebackPath::label(string& pth){
        
    if (pth.size()>0){
      pth.clear();
    }
        
        
    for(size_t k=trace_path.size()-1; k!=SIZE_MAX; k--){
      state* st = hmm->getState(trace_path[k]);
      pth+=st->getLabel();
    }
    return;
  }
    
  //!Get names of traceback path
  //!\param [out] pth vector<string>
  void TracebackPath::name(vector<string>& pth){
                
    if ( hmm==NULL ){
      cerr << "Model is NULL.  traceback::name(...) must have valid HMM model defined.\n";
      exit(2);
    }
                
    for(size_t k = trace_path.size()-1; k != SIZE_MAX; k--){
      state* st = hmm->getState(trace_path[k]);
      pth.push_back(st->getName());
    }
    return;
  }
    
    
  //!Get GFF output of traceback path
  //!\param [out] pth
  void TracebackPath::gff(vector<gff_feature>& pth,string& sequenceName){
    string current_label="";
    long long start=0;
    size_t path_size=size();
                
    if ( hmm==NULL ){
      cerr << "Model is NULL.  traceback::gff(...) must have valid HMM model defined.\n";
      exit(2);
    }
        
    for(size_t k = path_size-1;k != SIZE_MAX; k--){
      // state* st = hmm->getState(trace_path[k]);
      string new_label="FOO";//st->getGFF();
      if (new_label.compare("")==0){
	if (start>0){
	  gff_feature ln;
	  ln.seqname=sequenceName;
	  ln.source=hmm->getName();
	  ln.feature=current_label;
	  ln.start=start;
	  ln.end=path_size-(k+1);
	  ln.score='.';
	  ln.strand='+';
	  ln.frame='.';
	  ln.attribute="";
                    
	  pth.push_back(ln);
                    
	  start=0;
	  current_label=new_label;
	}
	else {
	  continue;
	}
      }
      else {
	if(k==0){
	  gff_feature ln;
	  ln.seqname=sequenceName;
	  ln.source=hmm->getName();
	  ln.feature=current_label;
	  ln.start=start;
	  ln.end=path_size-(k+1);
	  ln.score='.';
	  ln.strand='+';
	  ln.frame='.';
	  ln.attribute="";
                    
	  pth.push_back(ln);
	}
	else if (start==0){
	  start=path_size-k;
	  current_label=new_label;
	}
	else if (new_label.compare(current_label)==0){
	  continue;
	}
	else {
	  gff_feature ln;
	  ln.seqname=sequenceName;
	  ln.source=hmm->getName();
	  ln.feature=current_label;
	  ln.start=start;
	  ln.end=path_size-(k+1);
	  ln.score='.';
	  ln.strand='+';
	  ln.frame='.';
	  ln.attribute="";
                    
	  pth.push_back(ln);
                    
	  start=path_size-k;
	  current_label=new_label;
	}
      }
            
    }
        
    return;
  }
    
    
  //!Print the path to stdout
  void TracebackPath::print_path() const{
    int line=0;
    for(size_t k = this->size()-1; k != SIZE_MAX; k--){
      // for(size_t k = 0; k < this->size(); k++) { // REVERSE
      cout << trace_path[k] << " ";
      line++;
    }
    cout << endl;
  }

  //!Print the path to file stream
  void TracebackPath::fprint_path(ofstream &file){
    int line=0;
    for(size_t k=this->size()-1;k != SIZE_MAX; k--){
      // for(size_t k = 0; k < this->size(); k++) { // REVERSE
      file << trace_path[k]<< " ";
      line++;
    }
    file << endl;
  }

  //!Check to see if paths are the same
  bool TracebackPath::operator== (const TracebackPath &rhs) const{
    if (rhs.trace_path==trace_path){
      return true;
    }
    else {
      return false;
    }
  }

  //!Comparison operators for path
  bool TracebackPath::operator<  (const TracebackPath &rhs ) const{
    if (trace_path<rhs.trace_path){
      return true;
    }
    else {
      return false;
    }
  }

  //!Comparison operators for path
  bool TracebackPath::operator>  (const TracebackPath &rhs) const{
    if (trace_path>rhs.trace_path){
      return true;
    }
    else {
      return false;
    }
  }

  //!Comparison operators for path
  bool TracebackPath::operator<=  (const TracebackPath &rhs) const{
    if (trace_path<=rhs.trace_path){
      return true;
    }
    else {
      return false;
    }
  }

  //!Comparison operators for path
  bool TracebackPath::operator>=  (const TracebackPath &rhs) const{
    if (trace_path>=rhs.trace_path){
      return true;
    }
    else {
      return false;
    }
  }


  //!Print TracebackPath labels to stdout
  void TracebackPath::print_label() const {
    assert(hmm);
    for(size_t k=trace_path.size()-1; k!=SIZE_MAX; k--) {
      state* st = hmm->getState(trace_path[k]);
      cout << st->getLabel() << " ";
      cout.flush();
    }
    cout << endl;
  }

  //!Outputs the gff formatted output for the traceback to stdout
  void TracebackPath::print_gff(string sequence_name, double score, int ranking, int times, double posterior) const {
    string current_label="";
    long long start=0;
    size_t path_size=this->size();
                
    if ( hmm==NULL ){
      cerr << "Model is NULL.  traceback::print_gff(...) must have valid HMM model defined.\n";
      exit(2);
    }
        
    for(size_t k=path_size-1;k != SIZE_MAX;k--){
      // state* st = hmm->getState(trace_path[k]);
      string new_label="FOOP";//st->getGFF();
      if (new_label.compare("")==0){
	if (start>0){
	  cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
	  start=0;
	  current_label=new_label;
	}
	else {
	  continue;
	}
      }
      else {
	if(k==0){
	  cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
	}
	else if (start==0){
	  start=path_size-k;
	  current_label=new_label;
	}
	else if (new_label.compare(current_label)==0){
	  continue;
	}
	else {
	  cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t" << score << "\t+\t.\tRank:"<<ranking<<",Counts:" << times<<",Posterior:"<<posterior<<endl;
	  start=path_size-k;
	  current_label=new_label;
	}
      }
            
    }
  }   


  //!outputs the gff formatted output for the traceback
  void TracebackPath::print_gff(string sequence_name) const {
    string current_label="";
    long long start=0;
    size_t path_size=size();
                
    if (sequence_name[0] == '>'){
      sequence_name = sequence_name.substr(1);
                        
      if ( hmm==NULL ){
	cerr << "Model is NULL.  traceback::print_gff(...) must have valid HMM model defined.\n";
	exit(2);
      }
    }
        
    for(size_t k = path_size-1; k != SIZE_MAX; k--){
      // state* st = hmm->getState(trace_path[k]);
      string new_label="FOOP";//st->getGFF();
                        
      //If no label then print 
      if (new_label.compare("")==0){
	if (start>0){
	  cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+\t."<<endl;
	  start=0;
	  current_label=new_label;
	}
	else {
	  continue;
	}
      }
      else {
	if (start==0){
	  start=path_size-k;
	  current_label=new_label;
	}
	else if (new_label.compare(current_label)==0){
	  if(k==0){
	    cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size << "\t.\t+\t."<<endl;
                                                
	  }
                                        
	  continue;
	}
	else {
	  cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size-(k+1) << "\t.\t+\t."<<endl;
                                        
	  start=path_size-k;
	  current_label=new_label;
                                        
	  if(k==0){
	    cout << sequence_name << "\tStochHMM\t" << current_label <<"\t"<< start << "\t" << path_size << "\t.\t+\t."<<endl;
	  }
                    
	}
      }
            
    }
        
    cout << endl<<endl;
  }


  /* Need to re-write to handle new HMM Types
  //fix to handle higher order model 
  double TracebackPath::path_prob (const HMM &model){
  int size= trace_path.size();
        
  int seq_size= model.seq_size();
  int alpha_size=model.alpha.size();
        
        
  if (size!=seq_size){
  cerr << "Sequence and Path different sizes\n";
  exit(1);
  }
        
  vector<vector<double> > log_emm=convert_order(model, trace_path[size-1], 0);
  double prob=model.initial.get_trans(trace_path[size-1]) + log_emm[0][seq.val(0)];
        
        
  for (unsigned int i=1;i<size;i++){
  int index=0;
            
  if (model.states[trace_path[size-i-1]].order>i){
  vector<vector<double> > emmiss=convert_order(model, trace_path[size-i-1], i);
  for(int n=i;n>=1;n--){
  index+=POWER[n-1][alpha_size-1]*seq.val(i-n);
  }
                
  //cout <<i<<"\t"<<prob <<endl;
  prob+=model.states[trace_path[size-i]].get_trans(trace_path[size-i-1])+emmiss[index][seq.val(i)];
  }
            
            
  else{
                
  for(int n=model.states[trace_path[size-i-1]].order;n>=1;n--){
  index+=POWER[n-1][alpha_size-1]*seq.val(i-n);
  //cout <<"Index\t"<< index<<"\t"<< k << endl;
  }
  //cout <<i<<"\t"<<prob <<endl;
  prob+=model.states[trace_path[size-i]].get_trans(trace_path[size-i-1])+model.states[trace_path[size-i-1]].log_emm[index][seq.val(i)];
  }
  }
        
  return prob;    
  }
  */


  //!Create multiTraceback()
  multiTraceback::multiTraceback(){
    maxSize=0;
    vectorIterator=0;
    isFinalized=false;
    table=NULL;
    isFinalized=false;
  }

  //!Destroy multiTraceback
  multiTraceback::~multiTraceback(){
    if(table!=NULL){
      delete table;
    }
  }

  //!Set position in multiTraceback to the beginning
  void multiTraceback::begin(){
    vectorIterator=0;
    return;
  }
    
  //!Set position in multiTraceback to the ending
  void multiTraceback::end(){
    vectorIterator=paths.size();
    return;
  }
    
  //!Increment the iterator to next position
  void multiTraceback::operator++(){
    if (vectorIterator<maxSize){
      vectorIterator++;
    }
    return;
  }

  //!Decrement the iterator to previous position
  void multiTraceback::operator--(){
    if (vectorIterator>0){
      vectorIterator--;
    }
    return;
  }
    
  //!Set iterator to index val
  //!\param val Index value to set
  void multiTraceback::operator=(size_t val){
    if (val<=maxSize){
      vectorIterator=val;
    }
    return;
  }
    
    
  //!Get TracebackPath at index position
  //! \param val Index position
  TracebackPath multiTraceback::operator[](size_t val){
    return (*pathAccess[val]).first;
  }
    
    
  //!Get TracebackPath at currently set index in multiTraceback
  TracebackPath multiTraceback::path(){
    return (*pathAccess[vectorIterator]).first;
  }

  //!Get the number times that TracebackPath was recorded in multiple traceback
  int multiTraceback::counts(){
    return (*pathAccess[vectorIterator]).second;
  }

    
  //!Add TracebackPath to multiTraceback
  //!\param path Traceback path to add
  void multiTraceback::assign(TracebackPath& path){
    paths[path]++;
    return;
  }


  //!Sorts the multiTraceback by the number of time a particular tracback path occurred
  void multiTraceback::finalize(){
    assert(!isFinalized);
    maxSize=paths.size();
    vectorIterator=0;

    // push a pointer to each entry in paths into pathAccess
    map<TracebackPath,int>::iterator pathsIterator;
    for(pathsIterator=paths.begin();pathsIterator!=paths.end();pathsIterator++){
      pathAccess.push_back(pathsIterator);
    }

    // then sort it
    sort(pathAccess.begin(),pathAccess.end(),sortTBVec);
    isFinalized=true;
    return;
  }
    
  //!Generate a hit table from a multiple traceback paths
  //! Hit table is 2D table describing how many times a state was called at a particular position in the sequence
  heatTable* multiTraceback::get_hit_table(){
    if (table!=NULL){
      delete table;
    }
        
    // Over the length of the sequence
    model* hmm = ((*pathAccess[0]).first).getModel();
    size_t sequenceSize=((*pathAccess[0]).first).size();
    size_t stateSize=hmm->n_states();
        
        
    vector<int> states(stateSize,0);
    table = new heatTable(sequenceSize,states);
        
    map<TracebackPath,int>::iterator it;

    for( it =paths.begin(); it!=paths.end();it++){
      int count = (*it).second;
      // for(size_t position=sequenceSize-1; position!=SIZE_MAX; position--){
      for(size_t position=0;position<sequenceSize;position++){
	int tbState=(*it).first[position];
	(*table)[sequenceSize - position - 1][tbState]+=count;
      }
    }
    return table;
  }
    
    
  void multiTraceback::print_hits(){
    if (table==NULL){
      get_hit_table();
    }
        
    string header_row = "pos";
    model* hmm = ((*pathAccess[0]).first).getModel();
    for (size_t state_iter =0; state_iter<hmm->n_states(); state_iter++){
      header_row+="\t";
      header_row+=hmm->getStateName(state_iter);
    }

    cout << header_row << endl;
        
    for(size_t position = 0; position < table->size(); position++){
      string line = join((*table)[position], '\t');
      cout << position+1 << "\t" << line << endl;
    }
        
    return;
  }
    
    
  void multiTraceback::print_path(string format, string *header){
    assert(isFinalized);
    cout << "counts" << "          path" << endl;
    for(size_t iter=0; iter<this->size(); iter++){
      cout << (*pathAccess[iter]).second << "         ";
      if (format=="path")
	(*pathAccess[iter]).first.print_path();
      else if (format=="label")
	(*pathAccess[iter]).first.print_label();
      else if (format=="gff")
	(*pathAccess[iter]).first.print_gff(*header);
      else
	assert(0);
      // cout << endl;
    }
    return;
  }
    
  void multiTraceback::print_label(){
    print_path("label");
  }
    
  void multiTraceback::print_gff(string& header){
    print_path("gff", &header);
  }
    
    
  bool sortTBVec(map<TracebackPath,int>::iterator lhs, map<TracebackPath,int>::iterator rhs)
  {
    return ((*lhs).second < (*rhs).second);
  }

}
