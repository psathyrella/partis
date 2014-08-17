#ifndef EMM_H
#define EMM_H
    
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include "track.h"
#include "sequences.h"
#include "lexicalTable.h"

using namespace std;
namespace StochHMM {
class emm {
  friend class state;
  friend class model;
public:
  emm();
  bool parse(string&, tracks&);
  ~emm();
              
  // accessors
  inline bool isSimple() { return true; }
  inline bool isComplex() { return false; }
  bool pair() { return pair_; }
  size_t get_n_tracks() { return tracks_->size(); }
  inline double score(sequence& seq, size_t pos) { return scores.getValue(seq, pos); }
  inline double score(sequences& seqs, size_t pos) { return scores.getValue(seqs, pos); }
  
  inline lexicalTable* getTables(){ return &scores; }

  inline void print(){ cout << stringify() << endl; }
  string stringify();
private:
  bool pair_;
  lexicalTable scores;
  vector<track*>* tracks_;             //Tracks used
  vector<size_t>* track_indices; //Indices of tracks used
};
}
#endif
