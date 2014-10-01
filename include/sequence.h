#ifndef STOCHHMM_SEQUENCE_H
#define STOCHHMM_SEQUENCE_H

#include <vector>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <sstream>
#include <cassert>
#include "track.h"

using namespace std;
namespace stochhmm {

class Sequence {
  friend class Sequences;
public:
  Sequence() { assert(0); }  // do not allow (in order to make sure things are properly initialized)
  Sequence(string& seq_str, track* trk, string name);
  Sequence(const Sequence&);
  ~Sequence();
      
  Sequence getSubSequence(size_t pos, size_t len);
      
  string* getUndigitized();
  inline string name() const { return name_; }
  inline uint8_t getValue(size_t pos) const { return (*seq_)[pos]; }  // get digitized value at <pos>
  inline string getSymbol(size_t pos) const { return track_->getAlpha((*seq_)[pos]); }  // get undigitized value at <pos>
  inline size_t size() const { return seq_->size(); }
  inline track* getTrack() const { return track_; }

  string stringify(); // Get sequence as string
  string stringifyWOHeader(); //Get sequence without [Header information
  inline void print() { cout << stringify() << endl; }
  string undigitize();
              
  bool getFasta(ifstream&, track*);
  bool getFastq(ifstream&, track*);
  bool digitize();
  inline vector<uint8_t>* getDigitalSeq(){return seq_;}
  vector<uint8_t> getDigitalSubSeq(size_t istart, size_t istop);
  
  inline uint8_t operator[](size_t index){return (*seq_)[index];}
  void clear();
private:
  bool _digitize();  //Digitize the sequence
  string name_;
  string header_;
  string undigitized_;  //Undigitized sequence
  track* track_; //Ptr to track describing alphabet and type. NOTE we don't own this pointer, i.e. we don't delete it when we die
  vector<uint8_t>* seq_; // Digitized Sequence
};
}
#endif
