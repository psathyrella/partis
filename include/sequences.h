#ifndef HAM_SEQUENCE_H
#define HAM_SEQUENCE_H

#include "track.h"

using namespace std;
namespace ham {

class Sequence {
  friend class Sequences;
public:
  Sequence(string& seq_str, Track* trk, string name);
  Sequence(const Sequence&);
  ~Sequence();
      
  Sequence GetSubSequence(size_t pos, size_t len);
  inline string name() const { return name_; }
  inline uint8_t getValue(size_t pos) const { return (*seq_)[pos]; }  // get digitized value at <pos>
  inline string getSymbol(size_t pos) const { return track_->getAlpha((*seq_)[pos]); }  // get undigitized value at <pos>
  inline size_t size() const { return seq_->size(); }
  inline Track* getTrack() const { return track_; }

  string stringify(); // Get sequence as string
  string stringifyWOHeader(); //Get sequence without [Header information
  inline void print() { cout << stringify() << endl; }
  string undigitize();
              
  bool getFasta(ifstream&, Track*);
  bool getFastq(ifstream&, Track*);
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
  Track* track_; //Ptr to track describing alphabet and type. NOTE we don't own this pointer, i.e. we don't delete it when we die
  vector<uint8_t>* seq_; // Digitized Sequence
};

// ----------------------------------------------------------------------------------------
class Sequences {
public:
  Sequences():sequence_length_(0){}
  ~Sequences();
  short getValue(int , size_t); // Get digitized value for sequence track(i) in jth position
  Sequence* getSeq(size_t); // Return ith sequence
  Sequence& operator[](size_t index) {return *seqs[index]; }
  string& getHeader(size_t iter=0);
  // string* getUndigitized(size_t iter) { return seqs[iter]->getUndigitized(); }
  size_t n_seqs() { return seqs.size(); }
  size_t GetSequenceLength() { return sequence_length_;}
  Sequences GetSubSequences(size_t pos, size_t len);
  string stringifyWOHeader(); //! Get string of sequences
  string stringify(); //! Get string of sequences
  void print() { cout << stringify() << endl; }
  string undigitize(); //! Get sequence based on alphabet
  void addSeq(Sequence* sq);
private:
  vector<Sequence*> seqs;
  size_t sequence_length_; // length of the sequences (required to be the same for all)
};

}
#endif
