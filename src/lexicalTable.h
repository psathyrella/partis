#ifndef StochHMM_Lexical_h
#define StochHMM_Lexical_h

#include <string>
// #include <iomanip>
#include <vector>
#include <cassert>
// #include <math.h>
// #include <ctype.h>
// #include <algorithm>
// #include <stdint.h>
// #include <stdlib.h>
// #include "Eigen/Dense"
#include "track.h"
// #include "index.h"
// #include "weight.h"
#include "sequences.h"
#include "stochTypes.h"

using namespace std;
namespace StochHMM{
      
class lexicalTable{
public:
  lexicalTable();
  ~lexicalTable();

  void addTrack(track* trk, int order) { tracks.push_back(trk); }
  void AddColumn(vector<double> logprobs);
      
  inline double getValue(sequences &seqs, size_t pos) { return (*logProb)[seqs[0][pos]][seqs[1][pos]]; }  // sequences *must* have length of two
  inline double getValue(sequence &seq, size_t pos) { return (*logProb)[0][seq[pos]]; }
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
