#ifndef HAM_TRACK_H
#define HAM_TRACK_H

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <sstream>
#include <cassert>
#include "text.h"
#include "types.h"
#include "mathutils.h"

using namespace std;

namespace ham{
  class track;
  class tracks;

  // ----------------------------------------------------------------------------------------
  // ! \class ambigCharacter
  // !  \brief Define the ambiguous characters symbol and index number for digitizing ambiguous characters in the sequence
  // !  For example in DNA N = [ACGT] = [0,1,2,3]
  class ambigCharacter {
  public:
    ambigCharacter(track*, string&, vector<string>& ); //track, ambiguous character, unambiguous characters
    inline string getSymbol(){return symbol;};
    inline vector<size_t>& getDef(){return setDefinition;};
  private:
    string symbol;  //Ambiguous character definition
    vector<size_t> setDefinition; //Set of letters by digital value that ambiguous character defines
  };
    
  // ----------------------------------------------------------------------------------------
  //! \class track
  //! Defines types of data (real-value, text-sequence) used in the model
  //! and the alphabet that a text-sequence uses.  Tracks are used to digitize
  //! the sequence before decoding in HMM
  class track {
  public:
    track();
    track(string, size_t, vector<string>&);
        
    friend class State;
    friend class model;
    friend class tracks;
        
    // bool parse(string&);
    bool parseAmbiguous(string&);
    inline void setName(string nm){name=nm;};
    inline void setDescription(string& desc){description=desc;};
    inline void setIndex(size_t indx){
      if (indx<numeric_limits<size_t>::max()){
	trackIndex=indx;
	return;
      }
      cerr << "Track index: " << indx << " is OUT_OF_RANGE\n";
      exit(1);
    }
    inline void setAlphaType(trackType typ){alpha_type=typ;};
                
    bool addAlphabetChar(string);
    bool addAlphabetChar(const char *);
    bool addAlphabetChar(vector<string>&);
    bool addAlphabetChar(size_t chSize, const char * characters[]);
    bool addAlphabetChar(string& character, string& complement);
    bool addAlphabetChar(size_t chSize, const char* characters[], const char* complements[]);
    bool addAlphabetChar(vector<string>& characters , vector<string>& complements);
                
    void addComplement(string&, string&);
    void addComplement(const char *, const char *);
    bool addComplement(vector<string>& characters, vector<string>& complements);
        
    //! Set ambiguous character flag to true
    //! This will allow ambiguous characters to be processed in sequence
    //! Without this flag, only strict track characters or values are allowed
    inline void setAmbiguous(){ambiguous=true; return;};
    void addAmbiguous(string&,vector<string>&);
        
        
    //ACCESSOR
        
    //! Get the name of the track
    //! \return string Name of the track
    inline string getName(){return name;};
        
    //! Get the description of the track
    //! \return string Description of the track
    inline string getDescription(){return description;};
        
    //! Get the index of the track
    //! \return int Index of the track
    inline size_t getIndex(){
      if (trackIndex<numeric_limits<size_t>::max()){
	return trackIndex;
      }
      cerr << "Track Index is not set in track.  Set the track index with setIndex(size_t indx) before calling.\n";
      exit(1);
    };
        
    //! Get the alphabet type of the track
    //! \return trackType Alphabet type of the track
    //! \sa enum trackType
    inline trackType getAlphaType(){return alpha_type;}
        
    //! Get the number of characters defined in the track
    //! \return size_t Number of characters/words defined in the track
    inline size_t getAlphaSize(){return alphabet.size();};
                
    //! Get alphabet size including ambiguous characters
    inline size_t getTotalAlphabetSize(){return symbolIndices.size();}
        
    //! Get the size of the largest alphabet word
    //! \return size_t
    inline size_t getAlphaMax(){return maxSize;};
        
        
    string getAlpha(size_t);
        
        
    uint8_t symbolIndex(const string&);
    uint8_t symbolIndex(unsigned char);
        
    uint8_t getComplementIndex(uint8_t val);
    uint8_t getComplementIndex(string&);
        
                
    string getComplementSymbol(string& character);
    string getComplementSymbol(uint8_t value);

    inline bool isComplementDefined(){return complementSet;}
        
        
    //! Check to see if ambiguous flag is set for the track
    //! \return true if ambiguous flag is set to handle ambig. characters
    //! \return false if not set
    inline bool isAmbiguousSet(){return ambiguous;};
        
        
    //! Check to see if the track is AlphaNum type and not a REAL Track
    //! \return true if it is AlphaNum type
    //! \return false if it is a REAL type
    inline bool isAlpha(){if (alpha_type == ALPHA_NUM){return true;} else {return false;}};
        
    //! Get the number of ambiguous characters that are defined
    //! \return size_t Number of ambiguous characters/words defined
    inline size_t getAmbiguousSize(){return ambiguousSymbols.size();};
        
    string getAmbiguousCharacter(size_t);
        
        
    //! Get the indices of characters that an ambiguous character represents
    //! \return vector<int>
    inline vector<size_t>& getAmbiguousSet(uint8_t val){return ambiguousSymbols[(val-max_unambiguous)-1].getDef();}
        
    inline vector<size_t>& getUnambiguousSet(){ return unambiguous; }
        
                
    string stringify();
    string stringifyAmbig();
                
    //! Print the string representation of the track to stdout
    inline void print(){ cout << stringify() << endl;}
        
    string convertIndexToWord(size_t,size_t);
    void convertIndexToDigital(size_t,size_t,uint8_t*);
    size_t convertAlphaWordToIndex(string);
    size_t convertDigiWordToIndex(vector<uint8_t>);
                
    size_t n_seqs;  // number of sequences for this track (eg two for a pair hmm)

    //!Check if the TrackFunction is defined for this track
    //!\return true if the track has a trackFunc defined
    inline bool isTrackFuncDefined(){return trackFunctionDefined;}
        
    //! Get name of TrackFunc defined for track
    //! \return string Name of trackFunc defined
    inline string getTrackFunction(){return trackFunction;}
        
    //! Get name of Track to use for trackFunc
    inline string getTrackToUse(){return trackToUse;}
                
    inline uint8_t getMaxUnambiguous(){return max_unambiguous;}
    inline uint8_t getMaxAmbiguous(){return max_ambiguous;}
        
  private:
    string name;       /* Track Name */
    string description;        /* Track Desc */
    size_t trackIndex;     /*track number*/
        
    trackType alpha_type;   /* Track Type 1=string, 2=real_number  0=uninitialized*/
        
    //! Track Functions for defining Real Number Tracks
    bool trackFunctionDefined;
    bool complementSet;
        
    string trackToUse;
    string trackFunction;
        
    vector<string> alphabet;  //Contains the corresponding symbol,letter, word that is referenced in the seq by index
    map<uint8_t,uint8_t> complementAlphabet;
        
    size_t maxSize;  //Maximum size of the alphabet words
                
    uint8_t max_unambiguous;
    uint8_t max_ambiguous;
    vector<size_t> unambiguous;
        
    bool ambiguous; //Are ambiguous characters set
    int defaultAmbiguous; //Default ambiguous character
        
    //Contains the ambiguous characters defined by user corresponding to position in the
    //array. Where index 0=-1, 1=-2... so on.
    vector<ambigCharacter> ambiguousSymbols;
        
    map<string,uint8_t> symbolIndices;
    //map<char,uint8_t>* charIndices;
    vector<uint8_t>* charIndices;
        
    void _splitAmbiguousList(vector<pair<string ,vector<string> > >&, const string&);
  };
        
    
  class tracks{
  public:
        
    //MUTATOR
    void push_back(track*);
        
    //ACCESSOR
    size_t indexOf(const string&);
    size_t size(){return trks.size();};
    track* getTrack(const string&);
    bool isTrackDefined(const string&);
    track* operator[](size_t i){return trks[i];};
        
    void print();
    string stringify();
        
  private:
    vector<track*> trks;
    map<string,size_t> index;
  };
        
        
        
        
}
#endif /*TRACK_H*/
