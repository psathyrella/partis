//
//  Lexical.h
//  StochHMM
//
//  Created by Paul Lott on 4/2/12.
//  Copyright (c) 2012 University of California, Davis. All rights reserved.
//

#ifndef StochHMM_Lexical_h
#define StochHMM_Lexical_h
#include <string>
#include <iomanip>
#include <vector>
#include <cassert>
#include <math.h>
#include <ctype.h>
#include <algorithm>
#include <stdint.h>
#include <stdlib.h>
#include "Eigen/Dense"
#include "track.h"
#include "index.h"
#include "externalFuncs.h"
#include "weight.h"
#include "sequences.h"
#include "stochTypes.h"
//#include "simpleTable.h"

namespace StochHMM{
        
  //! \class lexicalTable
  //! \brief Lexical table stores the log2 probabilities for both emissions and lexical transitions
  //!
  class lexicalTable{
  private:
                
                
  public:
    lexicalTable();
    ~lexicalTable();
    void initialize_emission_table();

    // ----------------------------------------------------------------------------------------
    // mutators
    void addTrack(track*,int);
    inline void setUnkScoreType(unknownCharScoringType type){unknownScoreType=type;};  // not used! only kept to maintain interface to emm
    inline void setUnkScore(double val){unknownDefinedScore=val;};  // not used! only kept to maintain interface to emm
        
    // ----------------------------------------------------------------------------------------
    // accessors
    inline unknownCharScoringType getAmbScoringType(){return unknownScoreType;}  // not used! only kept to maintain interface to emm
    inline double getAmbDefinedScore(){return unknownDefinedScore;}  // not used! only kept to maintain interface to emm

    double getValue(sequences&, size_t);
    double getValue(sequence& , size_t);
    std::vector<std::vector<double> >* getCountsTable();
    std::vector<std::vector<double> >* getProbabilityTable();
    std::vector<std::vector<double> >* getLogProbabilityTable();

    //!Get pointer to track at index position of emission
    //!\param iter Index iterator of position
    //!\return track* Track in emission
    inline track* getTrack(size_t iter){return tracks[iter];};
        
    //!Get the number of tracks defined in emission
    //!\return size_t
    inline size_t getNTracks(){return tracks.size();};
        
    //!Get Orders of lexical emission will use for all tracks
    //!\return std::vector<int>
    inline std::vector<uint8_t>& getOrders(){return orders;};
    inline uint8_t getOrder(size_t i){return orders[i];}
        
    //! Get Log(prob) emission table
    //! \return std::vector<std::vector<double> >
    inline std::vector<std::vector<double> >& getLogEmm(){return *logProb;}
        
    //! Get the alphabet sizes for all tracks used in emission
    //! \return std::vector<int>
    inline uint8_t getAlphaSize(size_t i){return tracks[i]->getAlphaSize();}
    inline size_t getNumberOfAlphabets(){return tracks.size();}
        
    //!Increment counts
    inline void incrementCounts(size_t word_index, size_t char_index) { if (counts != NULL) (*counts)[word_index][char_index]++; }
        
    //!Increment counts by double
    inline void incrementCountsDouble(size_t word_index, size_t char_index, double val) { if (counts != NULL) (*counts)[word_index][char_index]+= val; }
        
    std::string stringify();
                
    std::string stringifyAmbig();
        
    void print();

  private:
    unknownCharScoringType unknownScoreType;  // not used! only kept to maintain interface to emm
    double unknownDefinedScore;  // not used! only kept to maintain interface to emm

    std::vector<track*> tracks;  // tracks which are used by emissions in this table
    std::vector<uint8_t> max_unambiguous;
    std::vector<uint8_t> orders;  // orders for each emission
    uint8_t max_order;
                
    size_t y_dim;
    size_t* x_subarray;
    size_t* y_subarray;
    // NOTE these are initialized in emm.cpp
    // first index: alphabet size
    // second index: (emission order)^(alphabet size)
    std::vector<std::vector<double> >* prob;     //p(x)
    std::vector<std::vector<double> >* counts;   //counts
    std::vector<std::vector<double> >* logProb;  //log2(P(x))
                
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> log_prob_matrix;
  };
}
#endif
