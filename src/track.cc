#include "track.h"
namespace stochhmm {
    
  // FIXME:  Check if ambiguous before allowing to get character.
  //!Create an ambiguous character
  //! For example in DNA N = [ACGT] = [0,1,2,3]
  //! \param tr Track to use for ambiguity code
  //! \param ambChar String representation of the character/symbol/word for ambiguous character
  //! \param defs vector of strings where each string is a non-ambiguous character defined in the track
  ambigCharacter::ambigCharacter(track* tr, string& ambChar, vector<string>& defs){
    symbol=ambChar;
    for(size_t i=0;i<defs.size();i++){
      setDefinition.push_back(tr->symbolIndex(defs[i]));
    }
    sort(setDefinition.begin(),setDefinition.end());
    return;
  }
        
    
  //!Create a track
  track::track(){
    alpha_type=UNDEFINED;
    trackIndex=numeric_limits<size_t>::max();
    ambiguous=false;
    defaultAmbiguous=-1;
    trackFunctionDefined= false;
    maxSize=0;
    max_ambiguous =0;
    max_unambiguous =0;
    complementSet=false;
    charIndices = NULL;
  }
        
  track::track(string name, size_t n_seqs, vector<string>& characters):n_seqs(n_seqs) {
    setName(name);
    assert(n_seqs == 1 || n_seqs == 2);
    trackIndex=numeric_limits<size_t>::max();
    ambiguous=false;
    defaultAmbiguous=-1;
    trackFunctionDefined= false;
    maxSize=0;
    max_ambiguous =0;
    max_unambiguous =0;
    complementSet=false;
    charIndices = NULL;
                
    addAlphabetChar(characters);
    setAlphaType(ALPHA_NUM);
  }
        
  string track::getAlpha(size_t iter) {
    assert(iter < alphabet.size());
    return alphabet[iter];
  }
    
  //FIXME: Have it check before adding value
  //! Add a letter/word symbol to the track
  //! \param character Word or symbol used in undigitized sequence
  bool track::addAlphabetChar(string character){
        
    if (alphabet.size() >= 255){
      cerr << "Alphabet limit reached.   Unable to add additional characters to the track:\t" << character << endl;
      return false;
    }
        
        
    if (character.size()>maxSize){maxSize=character.size();};
        
    alphabet.push_back(character);
        
    size_t index = alphabet.size()-1;
        
    symbolIndices[character]=index;
                
    max_unambiguous = index;
    unambiguous.push_back(index);
    setAlphaType(ALPHA_NUM);
        
    return true;
  }
    
  bool track::addAlphabetChar(const char *character){
    string string_character(character);
    setAlphaType(ALPHA_NUM);
                
    if (string_character.size()>maxSize){maxSize=string_character.size();};

    return addAlphabetChar(string_character);
  }
    
  bool track::addAlphabetChar(vector<string>& characters){
    for(size_t i=0;i<characters.size();i++){
      addAlphabetChar(characters[i]);
    }
                
    setAlphaType(ALPHA_NUM);
                
    return true;
  }
        
        
  bool track::addAlphabetChar(size_t chSize, const char * characters[]){
    for(size_t i =0; i < chSize; i++ ){
      addAlphabetChar(characters[i]);
    }
    return true;
  }
        
        
  bool track::addAlphabetChar(string& character, string& complement){
    addAlphabetChar(character);
    addComplement(character, complement);
    complementSet = true;
    setAlphaType(ALPHA_NUM);
    return true;
  }
        
  bool track::addAlphabetChar(size_t chSize, const char* characters[], const char* complements[]){
    for(size_t i=0;i<chSize;i++){
      addAlphabetChar(characters[i]);
      addComplement(characters[i], complements[i]);
    }
        
    complementSet = true;
    setAlphaType(ALPHA_NUM);
        
    return true;
  }
        
        
  //FIXME:  Need to fix the code below and test
  //Complements added by Ken
  bool track::addAlphabetChar(vector<string>& characters, vector<string>& complements){
        
    if (characters.size() != complements.size()){
      //Error Message
      cerr << "Number of Complement characters and Characters don't match.\n";
      return false;
    }
        
        
    for(size_t i=0;i<characters.size();i++){
      addAlphabetChar(characters[i]);
      addComplement(characters[i], complements[i]);
    }
        
    complementSet = true;
    setAlphaType(ALPHA_NUM);
        
    return true;
  }
    
  void track::addComplement(string& character, string& complement){
    int index = symbolIndex(character);
    int comp  = symbolIndex(complement);
    complementAlphabet[index]=comp;
    complementSet = true;
        
    return;
  }
    
  void track::addComplement(const char *character, const char *complement) {
    string string_character(character);
    string string_complement(complement);
    addComplement(string_character, string_complement);
        
    return;
  }
        
  bool track::addComplement(vector<string>& characters, vector<string>& complements){
                
    if (characters.size() != complements.size()){
      cerr << "Number of Complement characters and Characters don't match.\n";
      return false;
    }
                
    for(size_t i=0;i<characters.size();i++){
      addComplement(characters[i], complements[i]);
    }
                
    return true;
  }
    
        
  //!Add an ambiguous character/word definition to the track
  //! \param ambChar  word/symbol fore the ambiguous character
  void track::addAmbiguous(string& ambChar, vector<string>& defs){
    if (defaultAmbiguous==-1){
      defaultAmbiguous = max_unambiguous+1;
    }
    ambigCharacter amb(this,ambChar,defs);
    ambiguousSymbols.push_back(amb);
        
    int index = (int) symbolIndices.size();  //Get the index position and new digital reference value
                
    if (index >= 255){
      cerr << "Maximum number of discrete symbols reached at 255\n";
      exit(2);
    }
                
    symbolIndices[ambChar]=index;
    max_ambiguous = index;
    return;
  }
        
        
  //! Get symbol assigned integer value
  //! \param symbol word/letter/symbol that we want to get it's assigned integer value
  uint8_t track::symbolIndex(const string& symbol){
    if (symbolIndices.count(symbol)==0){  //If isn't found in the hash
      if (ambiguous){ //Return default character if ambiguous is set
	cerr << symbol << "not found in HMM definitions. Using default ambiguous character.\n"; 
	return defaultAmbiguous;
      }
      else{
	cerr << "Encountered an ambiguous character (" << symbol << ") in the sequence.  No ambiguous characters are allowed because they weren't set in the model.   To allow ambiguous characters, please add an \" Ambiguous Character Definition\" to the model" << endl;
	exit(1);
      }
    }
    else{
      return symbolIndices[symbol];
    }
  }
        
  //! Get symbol assigned integer value
  //! \param symbol word/letter/symbol that we want to get it's assigned integer value
  uint8_t track::symbolIndex(unsigned char symbol){
                
    if (maxSize != 1){
      cerr << "Track Max Symbols Size:\t" << maxSize << "\t Must use function track::symbolIndex(const string& symbol)\n";
    }
                
    if (charIndices == NULL){
      charIndices = new (nothrow) vector<uint8_t>(255,255);
      for(map<string,uint8_t>::iterator it = symbolIndices.begin(); it != symbolIndices.end(); it++){
	(*charIndices)[(it->first)[0]] = it->second;
      }
    }
                
                
    if ((*charIndices)[symbol]==255){  //If isn't found in the array
      if (ambiguous){ //Return default character if ambiguous is set
	return defaultAmbiguous;
      }
      else{
	cerr << "Encountered an ambiguous character (" << symbol << ") in the sequence.  No ambiguous characters are allowed because they weren't set in the model.   To allow ambiguous characters, please add an \" Ambiguous Character Definition\" to the model" << endl;
	exit(1);
      }
    }
    else{
      return (*charIndices)[symbol];
    }
  }
    
    
  // //FIXME: Change return value so only returns true if parse is OK
  // //! Parse a string representation of track to define a tracks parameters
  // //! \param txt Line from model that describes a track
  // //! \return true if the track was parsed properly
  // bool track::parse(string& txt) {
  //   stringList lst;
  //   stringList tag = extractTag(txt);
        
  //   lst.fromNext(txt);  // remove white space, and create the list using :, as delimiters
  //   setName(lst[0]);
  //   setDescription(lst.getComment());  // use any comments that were found as a description
  //   if (lst[1] == "pair") {  // pair hmm, i.e. we expect two seqs at a time for this track. Would be easy to make it work with more than two, but no point a.t.m.
  //     assert(0);  // er, I think this is old and not used TODO clean that shit out
  //   } else {
  //     n_seqs = 1;  // default to 1
  //   }
  //   setAlphaType(ALPHA_NUM);
    
  //   for(size_t i=1; i<lst.size(); i++)
  //     assert(addAlphabetChar(lst[i]));

  //   return true;
  // }
        
    
  //! Get the string representation of the track
  //! \return string Definition of the track as in model
  string track::stringify(){
    string output;
    output+=name + ":\t";
        
    if (alpha_type == ALPHA_NUM){
      output+=join(alphabet,',');
    }
    else{
      output+="REAL_NUMBER";
      if (trackFunctionDefined){
	output+="\t[FUNCTION:\t" + trackFunction;
	output+="\tUSE:\t" + trackToUse + "]";
      }
    }
    output+="\n";
        
    return output;
  }
    
  //!Get the string representation of the ambiguous character definitions as in model file
  //! \return string
  string track::stringifyAmbig(){
    string output;
    output+=name + ":\t";
    for (size_t i = max_unambiguous+1; i <= max_ambiguous; i++){
      if (i > (size_t) max_unambiguous+1){ output+= ",";}
            
      output+=getAmbiguousCharacter(i);
      output+="[";
            
      vector<size_t>& regChar = getAmbiguousSet(i);
      for(size_t k = 0; k<regChar.size();k++){
	if (k>0){output+=",";}
	output+=alphabet[regChar[k]];
      }
      output+="]";
    }
    return output;
  }
    
    
    
  string track::convertIndexToWord(size_t wordIndex, size_t order){
    string output="";
        
    if (order == 0){
      return "";
    }
        
    size_t currentOrder = order;
        
    while (currentOrder>1){
      double dreg=POWER[currentOrder-1][alphabet.size()-1];
      size_t temp = floor ((double) wordIndex / dreg);
      output+=alphabet[temp];
      if (maxSize!=1){
	output += ",";
      }
                        
      wordIndex-=temp*dreg;
      currentOrder--;
    }
        
    output+=alphabet[wordIndex];
        
    return output;
  }
        
        
  void track::convertIndexToDigital(size_t wordIndex, size_t order, uint8_t word[]){
    if (order == 0){
      word[0]=wordIndex;
      return;
    }
                
    cout << alphabet.size() << endl;
        
    size_t currentOrder = order;
        
    while (currentOrder>1){
      double dreg=POWER[currentOrder-1][alphabet.size()-1];
      size_t temp = floor ((double) wordIndex / dreg);
      word[currentOrder-1] = temp;
      wordIndex-=temp*dreg;
      currentOrder--;
    }
        
    word[0] = wordIndex;
    return;
  }
    

  //! \param word String representation of a subsequence
  //! \return unrolled index of this word
  size_t track::convertAlphaWordToIndex(string word) {
    size_t index(0);
    uint64_t exponent(0);
    uint64_t base(getAlphaSize());
    for (int ichar=(int)word.size()-1; ichar>=0; --ichar) {
      index += POWER[exponent][base-1] * symbolIndex(word[ichar]);  // POWER[b][a-1] = a**b
      exponent++;
    }
    return index;
  }

  //! \param word String representation of a subsequence
  //! \return unrolled index of this word
  size_t track::convertDigiWordToIndex(vector<uint8_t> word) {
    size_t index(0);
    uint64_t exponent(0);
    uint64_t base(getAlphaSize());
    for (int ichar=(int)word.size()-1; ichar>=0; --ichar) {
      index += POWER[exponent][base-1] * word[ichar];  // POWER[b][a-1] = a**b
      exponent++;
    }
    return index;
  }

  //FIXME: Change return value so only returns true if parse is OK
  //! Parse the ambiguous character definitions from model file
  //! \param txt String representation of ambiguous character definition as in model file
  //! \return true if the ambiguous characters were properly parsed
  bool track::parseAmbiguous(string& txt){
    vector<pair<string,vector<string> > > temp;
        
    _splitAmbiguousList(temp, txt);
    if (temp.size()==0){
      return false;
    }
    setAmbiguous();
    for (size_t i=0;i<temp.size();i++){
      addAmbiguous(temp[i].first, temp[i].second);
    }
        
    return true;
  }
    
  void track::_splitAmbiguousList(vector<pair<string,vector<string> > >& results, const string& text){
        
    size_t opening;
    size_t closing;
    size_t start=0;
        
    opening=text.find_first_of('[');
    while(opening!=string::npos){
      pair<string,vector<string> > amb;
      amb.first=text.substr(start,opening-start);
      clear_whitespace(amb.first, "\t ");
            
      closing=text.find_first_of(']',opening);
      if (closing!=string::npos){
	string tempString=text.substr(opening+1,closing-opening-1);
	split_line(amb.second, tempString);
      }
      start=text.find_first_not_of(',',closing+1);
      results.push_back(amb);
      opening=text.find_first_of('[',closing);
    }
        
    return;
  }
    
    
  //! Get the string representation of the ambigous character defined by integer value
  //! If ambiguous character isn't defined, return value is "*"
  //! \param val Integer value representing the ambiguous character
  string track::getAmbiguousCharacter(size_t val){
    if (getAmbiguousSize()==0){
      return "*";
    }
        
    return ambiguousSymbols[(val-max_unambiguous)-1].getSymbol();
  }
    
    
  //! Add track to tracks container
  //! \param tk Pointer to track to be added
  void tracks::push_back(track* tk){
    string& name= tk->name;
        
    if (!index.count(name)){
      index[tk->getName()]=trks.size();
      trks.push_back(tk);
    }
    else{
      cerr << "Track with name: " << name << " already exists.  Cannot add tracks with the same name\n";
      exit(1);
    }
        
  }
    
  //!Get iterator index of track by tracks name
  //! \param name Name of the track
  //! \return size_t Iterator to track within the tracks
  //! \return -1 if track doesn't exist in tracks
  size_t tracks::indexOf(const string& name){
    if (index.count(name)){
      return index[name];
    }
    else{
      return SIZE_MAX;
    }
  }
    
    
  //!Get pointer to track from the track name
  //! \param name Name of the track
  //! \return pointer to track if it exists, NULL otherwise
  track* tracks::getTrack(const string& name){
    if (index.count(name)){
      return trks[index[name]];
    }
    else{
      return NULL;
    }
  }
    
  bool tracks::isTrackDefined(const string& name){
    if (index.count(name)){
      return true;
    }
        
    return false;
  }
    
  //!Print the each track in tracks to stdout
  void tracks::print(){
    cout << stringify() << endl;
  }
    
  //! Get string representation of each track in tracks
  //! \return string String representation of tracks as in model file
  string tracks::stringify(){
    string trackString;
    string ambigString;
    string lnSep(50,'=');
    trackString+="TRACK SYMBOL DEFINITIONS\n" + lnSep + "\n";
        
    for(size_t i=0;i<trks.size();i++){
      trackString+=trks[i]->stringify();
      if (trks[i]->isAmbiguousSet()){
	ambigString+=trks[i]->stringifyAmbig();
      }
    }
    if (!ambigString.empty()){
      ambigString="AMBIGUOUS SYMBOL DEFINITIONS\n" + lnSep + "\n"+ ambigString;
      trackString+="\n" + ambigString + "\n";
    }
    else{
      trackString+="\n";
    }
        
    return trackString;
  }
    
    
  //! Get the complement alphabet character digitized value given a value
  //! \param val Value of character to get complement of
  //! \return int value of complement
  uint8_t track::getComplementIndex(uint8_t val){
    if (complementSet){
      if (complementAlphabet.count(val)){
	return complementAlphabet[val];
      }
      else{
	cerr<< "Complement of " << val << " is not set in track\n";
	return -1;
      }
    }
    else{
      cerr << "No complements are set in the track\n";
      exit(1);
    }
  }
    
    
    
  //! Get the complement alphabet character digitized value given the string
  //! \param character String of alphanumerical symbol
  //! \return int Defined complement string symbol of symbol
  uint8_t track::getComplementIndex(string& character){
    if (complementSet){
      if (symbolIndices.count(character)){
	int characterIndex = symbolIndex(character);
	return complementAlphabet[characterIndex];
      }
      else{
	cerr<< "Complement of " << character << " is not set in track\n";
	return -1;
      }
    }
    else{
      cerr << "No complements are set in the track\n";
      exit(1);
    }
  }
    
    
  //! Get the complement alphanumerical string of a given integer value;
  //! \param value Integer value of a symbol
  //! \return string Defined complement string symbol of symbol with digitized value
  string track::getComplementSymbol(uint8_t value){
    if (complementSet){
      if (complementAlphabet.count(value)){
	int complement_value = complementAlphabet[value];
	return getAlpha(complement_value);
      }
      else{
	cerr<< "Complement of " << value << " is not set in track\n";
	return " ";
      }
    }
    else{
      cerr << "No complements are set in the track\n";
      exit(1);
    }
                
  }
    
    
  //! Get the compelment alphabet character digitized value given the string
  //! \param character String of alphanumerical symbol
  //! \return string Defined complement string symbol
  string track::getComplementSymbol(string& character){
    if (complementSet){
      if (symbolIndices.count(character)){
	int characterIndex = symbolIndex(character);
	int complement_value = complementAlphabet[characterIndex];
	return getAlpha(complement_value);
      }
      else{
	cerr<< "Complement of " << character << " is not set in track\n";
	return " ";
      }
    }
    else{
      cerr << "No complements are set in the track\n";
      exit(1);
    }
  }
    
        
}
