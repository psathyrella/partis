#ifndef HAM_LEXICALTABLE_H
#define HAM_LEXICALTABLE_H

#include <string>
#include <vector>
#include <cassert>
#include "track.h"
#include "sequences.h"
#include "types.h"

using namespace std;
namespace ham {
      
class LexicalTable {
public:
  LexicalTable();
  void init();
  ~LexicalTable();

  void addTrack(track* trk, int order) { tracks.push_back(trk); }
  void AddColumn(vector<double> logprobs);

  inline double getValue(size_t letter) { assert(letter < (*logProb)[0].size()); return (*logProb)[0][letter]; }
  inline double getValue(size_t letter1, size_t letter2) { assert(letter1 < logProb->size()); assert(letter2 < (*logProb)[letter1].size()); return (*logProb)[letter1][letter2]; }
  inline double getValue(Sequences &seqs, size_t pos) { return (*logProb)[seqs[0][pos]][seqs[1][pos]]; }  // NOTE <seqs> *must* have length of two, and init() *must* have been called. I could check this, but I'm prematurely optimising. Good thing I'm not NASA, eh?
  inline double getValue(Sequence &seq, size_t pos) { return (*logProb)[0][seq[pos]]; }  // see NOTE above
  vector<vector<double> >* getProbabilityTable() { return prob; }
  vector<vector<double> >* getLogProbabilityTable() { return logProb; }
  vector<vector<double> >* getCountsTable() { return counts; }
  inline track* getTrack(size_t iter) { return tracks[iter]; }
  inline size_t getNTracks(){ return tracks.size(); }
  inline uint8_t getAlphaSize(size_t i){ return tracks[i]->getAlphaSize(); }
  inline size_t getNumberOfAlphabets(){ return tracks.size(); }
private:
  vector<track*> tracks;  // tracks which are used by emissions in this table
  // first index: seq 1 emission
  // second index: seq 2 emission
  // NOTE the two sequences cannot be distinguishable
  vector<vector<double> >* prob;     //p(x)
  vector<vector<double> >* counts;   //counts
  vector<vector<double> >* logProb;  //log2(P(x))
};
}
#endif
