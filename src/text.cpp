#include "text.h"
namespace stochhmm {
    
  //! Create a new stringList
  //! Remove whitespace set to true by default.
  stringList::stringList(){
    removeWS=true;
    return;
  }
    
  //! Create a new StringList
  //! \param txt String to be split
  //! \param ws Whitespace characters to remove from string before splitting
  //! \param del Delimiter characters to use to split string
  //! \param remove True will remove whitespace, False will leave whitespace
    
  stringList::stringList(string& txt, string& ws, string& del,bool remove):removeWS(remove), whitespace(ws),delimiters(del)
  {
    splitString(txt);
    return;
  }
    
  //! Searches for the string and returns bool if found or not
  //! \param txt String to search stringList
  bool stringList::contains(const string& txt){
    for (size_t i=0;i<lines.size();i++){
      if (lines[i].find(txt) != string::npos){
	return true;
      }
    }
    return false;
  }
        
  bool stringList::containsExact(const string& txt){
    for (size_t i=0;i<lines.size();i++){
      if (lines[i].compare(txt) == 0){
	return true;
      }
    }
    return false;
  }
    
  //! Searches the stringList for matching string and returns the index position of first match
  //! \param txt String to search stringList for
  size_t stringList::indexOf(const string&txt) {
    for (size_t i=0;i<lines.size();i++){
      if (lines[i].find(txt) != string::npos){
	return i;
      }
    }
    return SIZE_MAX;
  }
    
  //! Searches the stringList for matching string and returns the index position of a given string from the starting position
  //! \param txt String to search stringList for
  //! \param pos Position to start search from
  size_t stringList::indexOf(const string&txt,size_t pos){
    for (size_t i=pos;i<lines.size();i++){
      if (lines[i].find(txt) != string::npos){
	return i;
      }
    }
    return SIZE_MAX;
  }
    
  //! Removes Comments, Removes Whitespace, and splits the string
  //! Whitespace and Split delimiters are required to be previously set in stringList
  void stringList::processString(string& txt){
    lines.clear();
    comment.clear();
        
    //Extract Comments
    comment=parseComment(txt, '#');
        
    //Remove Whitespace
    if (removeWS){
      clear_whitespace(txt, whitespace);
    }
                
    splitString(txt);
        
    return;
  }
    
  //! Split string <txt> using member <delimiters> and store in member <lines>.
  //! The char delimiters are deleted from the returned value. <delimiters> contains
  //! one or more delimiters in a single string.
  //! \param txt String to be split
  void stringList::splitString(const string& txt){
        
    lines.clear();
    comment.clear();
        
    size_t found  = SIZE_MAX;
    size_t initial= SIZE_MAX;
    do {
      found++;
      initial=found;
            
      found=txt.find_first_of(delimiters.c_str(),found);
      string st=txt.substr(initial,found-initial);
            
      if (st.size()>0){
	lines.push_back(st);
      }
            
    } while (found!=string::npos);
        
    return;
  }
    
  //! Split string using delimiter
  //! Splits the string up according to character delimiters supplied in the delimiter string.
  //! The char delimiters are deleted from the returned value.  Delimiters allows you supply multiple
  //! delimiters in a single string.
  //! \param txt String to be split
  //! \param del Delimiters to use to split evaluated as single characters, not whole string
  void stringList::splitString(const string& txt, const string& del){
    delimiters=del;
    splitString(txt);
    return;
  }
    
  void stringList::splitString(const string& txt,size_t charSize){
    for(size_t i=0;i<txt.size();i+=charSize){
      lines.push_back(txt.substr(i,charSize));
    }
    return;
  }
    
  //! Non-destructive string split
  //! Splits the string up using an additional string as the delimiter. Unlike
  //! stringSplit the delimiter isn't deleted from the returned value;
  //! \param txt  String to be split
  //! \param del  Delimiter string to use to split
  void stringList::splitND(const string& txt,const string& del){
    lines.clear();
    comment.clear();
    size_t initial=txt.find(del);
        
    do {
      size_t found=txt.find(del,initial+1);
      string st=txt.substr(initial,found-initial);
      initial=found;
      if (st.size()>0){
	lines.push_back(st);
      }
            
    } while (initial!=string::npos);
        
    return;
  }
    
    
  //!Remove leading whitespace
  //!Removes the leading whitespace from the sequence
  //! \param ws String of whitespace characters to remove
  void stringList::removeLWS(const string& ws){
    for(size_t i=0;i<this->size();i++){
      removeLeadingWS((*this)[i],ws);
    }
    return;
  }
    
    
  //! Remove whitespace characters '\\t' and then splits string by delimiters ':' or '\\n'
  //! \param txt String to be split 
  bool stringList::fromTxt(string& txt){         
    setWhitespace("\t");
    setDelimiters(":\n");
    processString(txt);
        
    if (size()>0){
      return true;
    }
    else{
      return false;
    }
  }
    
  //! Remove whitespace characters '\\t' and then splits string by delimiters ':' or '\\n'
  //! Uses istream::getline(...,\\n) and splits string
  //! \param in istream reference 
  bool stringList::fromTxt(istream& in){
    string temp;
    getline(in, temp, '\n');
    return fromTxt(temp);
  }
    
    
  //! Remove whitespace characters '\t' and then splits string by delimiters ':', '\n', comma or space
  //! \param txt String to be split 
  bool stringList::fromTrack(string& txt){
    setWhitespace("\t");
    setDelimiters(":, \n");
    processString(txt);
        
    if (size()>0){
      return true;
    }
        
    return false;
  }
  //! Remove whitespace characters '\t' and then splits string by delimiters ':', '\n', comma or space
  //! \param in istream reference
  bool stringList::fromTrack(istream& in){
    string temp;
    getline(in, temp, '\n');
    return fromTrack(temp);
  }
    
  //! Remove whitespace characters '\t','\n' or space and then splits string by delimiters ':' or comma
  //! \param txt String to be split 
  bool stringList::fromNext(string& txt){
    setWhitespace("\t\n ");
    setDelimiters(":,");
    processString(txt);
        
    if (size()>0){
      return true;
    }
        
    return false;
  }
    
  //! Remove whitespace characters '\t','\n' or space and then splits string by delimiters ':' or comma
  //! \param txt istream reference
  bool stringList::fromNext(istream& in){
    string temp;
    getline(in, temp, '\n');
    return fromNext(temp);
  }
    
    
  //! Remove whitespace characters '\t' and then splits string by delimiters ':' or '\n'
  //! \param txt String to be split
  bool stringList::fromAlpha(const string& txt,size_t alpha){
    lines.clear();
    comment.clear();
                
    if (foundAlphaDelimiter(txt)){
      splitString(txt,";, \t");
    }
    else{
      splitString(txt,alpha);
    }
        
    if (this->size()>0){
      return true;
    }
        
    return false;
  }
    
  //! Returns true if text delimiters ":,[space]\\t" are found in the string
  //! \param txt String used to find delimiter
  bool stringList::foundAlphaDelimiter(const string& txt){
    for(size_t i=0;i<txt.size();i++){
      switch(txt[i]){
      case ':':
	return true;
      case ',':
	return true;
      case ' ':
	return true;
      case '\t':
	return true;
      default:
	break;
      }
            
    }
    return false;
  }
    
    

  //! Remove whitespace characters '\t' and then splits string by delimiters ':' or '\n'
  //! \param in istream stream
  //! \param alpha Size of alphabet
  bool stringList::fromAlpha(istream& in, size_t alpha){
    string temp;
    getline(in, temp, '\n');
    return fromAlpha(temp,alpha);
  }
    
    
  //! Splits a definition from model file
  //! \param txt string to split
  //! \param ws Whitespace characters to use
  //! \param del Delimiters characters to use
  bool stringList::fromDef(string& txt, string& ws, string& del){
    whitespace=ws;
    delimiters = del;
    processString(txt);
        
    if (size()>0){
      return true;
    }
    else{
      return false;
    }
  }
    
    
  //! Splits a definition from model file
  //! \param in istream to get line from
  //! \param ws Whitespace characters to use
  //! \param del Delimiters characters to use
  bool stringList::fromDef(istream& in, string& ws, string& del){
    whitespace=ws;
    delimiters = del;
    string temp;
    getline(in, temp, '\n');
    return fromDef(temp, ws, del);
  }
    
    
  string stringList::pop_ith(size_t pos){
    if (pos>=lines.size()){
      return "";
    }
        
    string temp=lines[pos];
    lines.erase(lines.begin()+pos);  // erase the <pos>th element of the (member) vector <lines>
    return temp;
  }
    
    
  //! Returns the values in the stringList as vector of doubles
  vector<double> stringList::toVecDouble(){
    vector<double> temp;
    for(size_t iter=0;iter<lines.size();iter++){
            
      double val;
            
      stringToDouble(lines[iter], val);
      temp.push_back(val);
    }
    return temp;
  }
    
  //! Returns the values in the stringList as vector of integers
  void stringList::toVecInt(vector<int>& ret_val){
        
    for(size_t iter=0;iter<lines.size();iter++){
            
      int val(0);
            
      if (!stringToInt(lines[iter], val)){
	cerr << "Couldn't convert " << lines[iter] << " to an integer\n";
      }
      else{
	ret_val.push_back(val);
      }
            
    }
        
    return;
  }

    
    
  //! Print each lines to stdout
  void stringList::print(){
    for(size_t i=0;i<lines.size();i++){
      cout << lines[i] << endl;
    }
    cout << "#" << comment << endl;
  }
    
    
  //! Joins the stringList into string using "\t" and return string
  //! \return string of stringList joined by "\\t"
  string stringList::stringify(){
    string output=join(lines,'\t');
    if (comment.size()>0){
      output+="#"+ comment;
    }
        
    return output;
  }
    
    
  //!Parses out the comments and stores the comment delimited by "#" 
  //!Comment can be accessed by command getComment();
  void stringList::removeComments(){
    for(size_t iter=0;iter<lines.size();iter++){
      comment+=parseComment(lines[iter], '#');
    }
  }
    
    
  //! Find first comment character and then return everything following the character
  string parseComment(string& txt, char commentChar){
    string comment;
    size_t commentPos=txt.find_first_of(commentChar);
    if (commentPos!=string::npos){
      comment=txt.substr(commentPos);
      txt=txt.substr(0,commentPos);
    }
    return comment;
  }

  //! Given a string, and a white space character, it will remove all the whitespace characters from the string
  //! \param input String to remove whitespace from
  //! \param white String of whitespace characters to remove
  void clear_whitespace(string &input,string white){
    size_t found;
    found=input.find_first_of(white);
    //found=input.find_first_of("\t\n ");
    while(found!=string::npos){
      input.erase(found,1);
      //found=input.find_first_of("\t\n ");
      found=input.find_first_of(white);
    }
    return;
  }
    
  //! Removes leading whitespace characters from a string
  //! \param txt String user wants to remove whitespace from
  //! \param ws  String containing whitespace characters to remove 
  void removeLeadingWS(string& txt,const string& ws){
    size_t start = txt.find_first_not_of(ws);
    if (start==string::npos){
      txt="";
    }
    else{
      txt=txt.substr(start);
    }
    return;
  }
    
    
  //! Parses key and value from a line
  //! Where key is delimited by <<KEY>> = Value
  //! \param txt String to extract key value from
  //! \param key String to assign key to
  //! \param value String to assign value to
  void getKeyValue(string& txt,string& key,string& value){
    size_t found=txt.find("<<");
        
    if (found==string::npos){
      removeLeadingWS(txt,"\t \n");
      value=txt;
      return;
    }
    else{
      size_t ending=txt.find(">>");
      if (ending!=string::npos){
	key=txt.substr(found+2,ending-(found+2));
	stringList lst;
	lst.splitString(txt,"=");
	removeLeadingWS(lst[1],"\t \n");
	value= lst[1];
	return;
      }
      else{
	cerr << "Missing closing brackets on Key\n";
	exit(1);
      }
            
    }
  }
    
    
  //! Splits string using delimiters and return stringList
  stringList& splitString(const string& txt, const string& delimiters){
    static stringList lst;
    lst.clear();
    lst.splitString(txt,delimiters);
    return lst;
  }
    
    
  //! Replace a given character with another character in a string
  //! \param txt String to use have characters replaced
  //! \param ch Character to search string for
  //! \param replaceCh Character to replace found ch with
  void replaceChar(string& txt, char ch, char replaceCh){
    size_t found = txt.find(ch);
    while(found!=string::npos){
      txt[found]=replaceCh;
      found=txt.find(ch);
    }
    return;
  }
    
  //! Converts a vector of ints into a string delimited by a character c
  //! \param input Vector of integers to be converted
  //! \param c Character to use as a delimiter
  string join(vector<int> &input, char c){
    string out;
    if (input.size()==0){
      out="";
      return out;
    }
    else if (input.size()==1){
      out=int_to_string(input[0]);
      return out;
    }
    else{
      out=int_to_string(input[0]);
      for(size_t i=1;i<input.size();i++){
	out+=c;
	out+=int_to_string(input[i]);
      }
      return out;
    }
  }
        
  //! Converts a vector of size_t into a string delimited by a character c
  //! \param input Vector of integers to be converted
  //! \param c Character to use as a delimiter
  string join(vector<size_t> &input, char c){
    string out;
    if (input.size()==0){
      out="";
      return out;
    }
    else if (input.size()==1){
      out=int_to_string(input[0]);
      return out;
    }
    else{
      out=int_to_string(input[0]);
      for(size_t i=1;i<input.size();i++){
	out+=c;
	out+=int_to_string(input[i]);
      }
      return out;
    }
  }

  //! Convert an integer to a string
  //! \param input Integer you want to convert to string;
  string int_to_string(int input){
    stringstream ss;
    ss << input;
    string s=ss.str();
    return s;
  }
    
    
  //! Convert an size_t to a string
  //! \param input Integer you want to convert to string;
  string int_to_string(size_t input){
    stringstream ss;
    ss << input;
    string s=ss.str();
    return s;
  }
    
    
        
  //! Convert a double to a string
  //! \param input Double you want to convert to a string
  string double_to_string(double input){
    stringstream ss;
    ss << input;
    string s=ss.str();
    return s;
  }
    
  //! Convert a double to a string
  //! \param input Double you want to convert to a string
  string double_to_string(float input){
    stringstream ss;
    ss << input;
    string s=ss.str();
    return s;
  }
    
    
  //!Convert string to integer
  //!\param txt Text representation of integer
  //!\param val Integer to be assigned
  //!\return true if conversion is valid
  //!\return false if conversion can't be performed
  bool stringToInt(string& txt, int& val){
    istringstream input(txt);
    if (!(input >> val)){
      return false;
    }
        
    return true;
  }
        
  //!Convert string to integer
  //!\param txt Text representation of integer
  //!\param val Integer to be assigned
  //!\return true if conversion is valid
  //!\return false if conversion can't be performed
  bool stringToInt(string& txt, size_t& val){
    istringstream input(txt);
    if (!(input >> val)){
      return false;
    }
        
    return true;
  }
    
    
  //!Convert string to double
  //!\param txt Text representation of double
  //!\param val Integer to be assigned
  //!\return true if conversion is valid
  //!\return false if conversion can't be performed
  bool stringToDouble(string& txt, double& val){
    istringstream input(txt);
    if(!(input >> val)){
      return false;
    }
    return true;
  }
    
    

        
  //! Converts a vector of shorts into a string delimited by a character c
  //! \param input Vector of shorts to be converted
  //! \param c Character to use as a delimiter
  string join(vector<short> &input, char c){
    string out;
    if (input.size()==0){
      out="";
      return out;
    }
    else if (input.size()==1){
      out=int_to_string(input[0]);
      return out;
    }
    else{
      out=int_to_string(input[0]);
      for(size_t i=1;i<input.size();i++){
	out+=c;
	out+=int_to_string(input[i]);
      }
      return out;
    }
  }

  //! Converts a vector of doubles into a string delimited by a character c
  //! \param input Vector of doubles to be converted
  //! \param c Character to use as a delimiter
  string join(vector<double> &input, char c){
    string out;
    if (input.size()==0){
      out="";
      return out;
    }
    else if (input.size()==1){
      out=double_to_string(input[0]);
      return out;
    }
    else{
      out=double_to_string(input[0]);
      for(size_t i=1;i<input.size();i++){
	out+=c;
	out+=double_to_string(input[i]);
      }
      return out;
    }
  }
    
  //! Converts a vector of strings into a string delimited by a character c
  //! \param input Vector of strings to be converted
  //! \param c Character to use as a delimiter
  string join(vector<string> &input, char c){
    string out;
    size_t sz=input.size();
    if (sz==1){
      out = input[0];
    }
    else if(sz>1){
      out+=input[0];
      for(size_t i=1;i<sz;i++){
	out+= c + input[i];
      }
    }
        
    return out;
  }
    
    
  //! Splits a line into a vector of string using delimiters ' \",[space]\\n\\t'
  //! \param line_list vector of strings to split input into
  //! \param input String to be split using delimiters ' \",[space]\\n\\t'
  void split_line(vector<string> &line_list,string &input){
    // NOTE does not push back to member variable <lines>

    //split line_list accoring to delimiters;
    size_t found=input.find_first_of("\", \n\t");
    while(found!=string::npos){
      if (found>0){
	line_list.push_back(input.substr(0,found));
	//cout << line_list.back() << endl;
	input=input.substr(found+1);
	//cout << input << endl;
      }
      else{
	input.erase(found,1);
	//cout << input <<endl;
      }
      found=input.find_first_of("\", \n\t");
      //cout <<input <<endl;
    }
    if (input.size()>0){
      line_list.push_back(input);
    }
    return;
  }
    
  //! Parse a line and extract a bracketed tag from the model file
  //! Returns a stringList which contains the tag split using ":\\t[space]"
  //! \param txt String to be have tag extracted from
  stringList extractTag(string& txt){
    stringList lst;
    pair<size_t,size_t> tagCoord = balanced_brackets(txt,"[]");
    if (tagCoord.first!=tagCoord.second){
      string tag = txt.substr(tagCoord.first+1,tagCoord.second-tagCoord.first-1);
      txt.erase(tagCoord.first,tagCoord.second-tagCoord.first+1);
      lst.splitString(tag,":\t ");
    }
    return lst;
  }

  //! Returns a pair of size_t values that describe the coordinates of between brackets
  //! \param text  String to use to search for brackets
  //! \param brackets String of length two containing opening and closing bracket
  //! \param offset Offset of where to start searching for balanced brackets
  pair<size_t,size_t> balanced_brackets(const string& text, const string& brackets, size_t offset){
    char opening = brackets[0];
    char closing = brackets[1];
        
        
    int currentTotal(0);
        
    size_t start;
    size_t found;
    found=text.find_first_of(opening,offset);
        
    if (found!=string::npos){
      start=found;
      currentTotal++;
    }
    else{
      return make_pair(0,0);
    }
        
    while(currentTotal!=0){
      found++;
      found=text.find_first_of(brackets,found);
      if (found!=string::npos){
	if (text[found]==opening){
	  currentTotal++;
	}
	else if (text[found]==closing){
	  currentTotal--;
	}
      }
      else{
	return make_pair(0,0);
      }
    }
        
    return make_pair(start,found);
  }
    
    
  //! Returns a pair of size_t values that describe the coordinates of between brackets from start of the string
  //! \param text  String to use to search for brackets
  //! \param brackets String of length two containing opening and closing bracket
  pair<size_t,size_t> balanced_brackets(const string& text, const string& brackets){
    return balanced_brackets(text,brackets,0);
  }



  //! Is the value of string numeric
  //! \param str  String to determine if it is numeric
  bool isNumeric(const string& str){
    size_t found;
    found = str.find_first_not_of("0123456789.-eE");
        
    if (found!=string::npos){
      return false;
    }
    else{
      return true;
    }
  }
    
  //!Slurp a file into a string 
  //! \param file Filename
  //! \return string that contains complete file
  string slurpFile(string& file){
    ifstream in(file.c_str(),ifstream::in);
    if (!in.good()){
      cerr << "File doesn't exist:" << file;
      exit(1);
    }
            
    stringstream sstr;
    sstr << in.rdbuf();
    return sstr.str();
  }

    
}
