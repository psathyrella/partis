#ifndef STOCHHMM_TEXT_H
#define STOCHHMM_TEXT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdlib.h>

using namespace std;
namespace stochhmm{
        
  /*! \defgroup Text  Text Handling
   *  
   *
   */
        
        
    
  /*! \class stringList
    \brief Stringlist is an list of strings that contains parsed comments and 
    can be used to split string or remove whitespace from a string
     
  */

  class stringList{
  public:
    stringList();
    stringList(string&,string&,string&,bool);
        
    //MUTATORS
    //! Set remove whitespace flag. When splitting string it will remove Whitespace first and then split the string.
    inline void setRemoveWS(){removeWS=true;}; //!<
        
        
    //! Unset remove whitespace flag. When splitting string it will remove Whitespace first and then split the string.
    inline void unsetRemoveWS(){removeWS=false;};
        
    //! Set whitespace character. When splitting string it will remove defined whitespace characters before splitting the string.
    //! \param txt String of whitespace characters to remove from sequence
    inline void setWhitespace(const string& txt){whitespace=txt;};
        
    //! Set delimiter characters.  String will be split on these characters
    //! \param txt String of text-delimiters to use in splitting string
    inline void setDelimiters(const string& txt){delimiters=txt;};
        
    void getLine(istream&);
    void processString(string&); //Remove comment, whitespace, and split string
    void splitString(const string&); //Split string with delimiter defined in stringList
    void splitString(const string&, const string&); //Split string with given delimiter
    void splitString(const string&, size_t);
    void splitND(const string&, const string&);
    void removeLWS(const string&);
    void removeComments();
        
    //! Clears all values from the stringList, including comments, whitespace, and delimiters
    inline void clear(){lines.clear();comment.clear();whitespace.clear();delimiters.clear();};
        
    bool fromTable(string&);
    bool fromTable(istream&);
    bool fromTrack(string&);
    bool fromTrack(istream&);
    bool fromAlpha(const string&,size_t);
    bool fromAlpha(istream&,size_t);
    bool fromTxt(string&);
    bool fromTxt(istream&);
    bool fromNext(string&);
    bool fromNext(istream&);
    bool fromDef(string&,string&, string&);
    bool fromDef(istream&,string&, string&);
        
    //! Add string to StringList
    //! \param txt String to add to the stringList
    inline void operator<< (string& txt){lines.push_back(txt);};
        
    //! Add string to StringList
    //! \param txt String to add to the stringList
    inline void push_back  (string& txt){lines.push_back(txt);};
        
    string pop_ith(size_t);  // erase the <pos>th 
        
    //! Copy StringList
    //! \param lst stringList to copy
    inline void operator= (stringList& lst){lines=lst.lines; comment=lst.comment; whitespace=lst.whitespace; delimiters=lst.delimiters;};
        
    //ACCESSORS
    //! Access string at iterator 
    //! \param iter Position to return;
    inline string& operator[](size_t iter){return lines[iter];};
        
    //! Return any comments parsed from the line
    inline string& getComment(){return comment;};
        
    //! Returns the list as a Vector of Doubles
    vector<double> toVecDouble();
        
    void toVecInt(vector<int>&);

        
    size_t indexOf(const string&);
    size_t indexOf(const string&,size_t);
    bool contains(const string&);
    bool containsExact(const string& txt);
    void print();
    string stringify();
        
    //! Return the amount of strings in the stringList 
    inline size_t size(){return lines.size();};
    
        
  private:
    void _splitIndividual(string&, size_t);
    bool foundAlphaDelimiter(const string&);
        
    vector<string> lines;
    string comment;
        
    bool removeWS;
    string whitespace;
    string delimiters;
  };


   
    
  stringList& splitString(const string&,const string&);
    
  void clear_whitespace(string&, string);
  void replaceChar(string&, char ,char);
  void removeLeadingWS(string&,const string&);
  string parseComment(string&,char);
  void getKeyValue(string&,string&,string&);

  string join(vector<int>&, char);
  string join(vector<size_t>&, char);
  string join(vector<short>&,char);
  string join(vector<double>&, char);
  string join(vector<string>&, char);

  string int_to_string(int);
  string int_to_string(size_t);
  string double_to_string(double);
  string double_to_string(float);
    
  bool stringToInt(string&, int&);
  bool stringToInt(string&, size_t&);
  bool stringToDouble(string&, double&);

  bool isNumeric(const string&);


  void split_line(vector<string>&, string& );

  stringList extractTag(string&);
  pair<size_t,size_t> balanced_brackets(const string&, const string& );
  pair<size_t,size_t> balanced_brackets(const string&, const string& , size_t );
    
  string slurpFile(string&);

  /*
    inline double convertToDouble (const string& s){
    istringstream i(s);
    double x;
    if (!(i>>x)){
    throw BadConversion("convertToDouble(\"" + s "\")"));
    }
    return x;
    }
  */
        
  /**@}*/ //End of Text Handling Group
    
}
#endif /*TEXT_H*/
