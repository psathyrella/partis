#ifndef HAM_TRACEBACKPATH_H
#define HAM_TRACEBACKPATH_H

#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <algorithm>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <cassert>

#include "text.h"
#include "model.h"
#include "mathutils.h"

using namespace std;

namespace ham {

  //! \struct gff_feature
  //!Each gff_feature represents a single GFF line
  struct gff_feature{
    string seqname;    //Column 1
    string source;     //Column 2
    string feature;    //Column 3
    size_t start;              //Column 4
    size_t end;                //Column 5
    char score;             //Column 6
    char strand;            //Column 7
    char frame;             //Column 8
    string attribute;  //Column 9
  };


  //! Perform traceback of traceback table
  //! Stores one traceback path for a sequence
  class TracebackPath{
  public:
    TracebackPath(Model*);
        
    //!Add state to traceback
    void push_back(int);
                
    //!Erase the traceback
    void clear();
                
    //!Get the size of traceback
    //!\return size_t
    size_t size() const;
                
    //!Returns the state index at a given position (it) within the traceback sequence
    inline int val(size_t it){
      if (it>=trace_path.size()){
	cerr << "Out of Range index\n ";
	exit(2);
      }
      return trace_path[it];
    }
                
    //!Get the model used for the decoding
    //! \return model
    inline Model* getModel() const {return hmm;};
    inline void setModel(Model* mdl){hmm=mdl;};
        
                
    //!Print the path to file stream
    void fprint_path(ofstream&);
        
    //! Get traceback as vector of state labels
    //! \param [out] vector<string> Vector of Labels
    void label(vector<string>&);
                
    //! Get traceback as string of state labels
    void label(string&);
        
                
    //! Get traceback as vector of gff_features
    //! \param[out] pth reference to vector of gff_feature
    //! \param[in] sequenceName Name of Sequence to be used in GFF
    void gff(vector<gff_feature>&,string&);
                
    //!Get names of traceback path
    //!\param [out] pth vector<string>
    void name(vector<string>&);
                
    //! Get the path to vector<int>
    //! \param [out] pth vector<int> that represents path
    void path(vector<int>&);
        
    //! Print the traceback path as path to stdout using cout
    //! Path numbers correspond to state index in model
    void print_path() const ;
                
    //! Print the traceback path as state labels
    //! State labels
    void print_label() const ;
                
    //!Outputs the gff formatted output for the traceback to stdout
    //!Allows user to provide additional information, that may be
    //!pertinent to stochastic tracebacks
    //!\param[in] sequence_name  Name of sequence used
    //!\param[in] score score to use in the GFF output
    //!\param[in] ranking  Rank of traceback
    //!\param[in] times Number of times that traceback occurred
    //!\param[in] posterior Posterior probability score
    void print_gff(string,double,int,int,double) const ;
        
    //!Outputs the gff formatted output for the traceback
    void print_gff(string) const ;
        
    //!Get the score that is associated with the traceback
    inline double getScore(){
      return score;
    }
        
    //!Set the score for the traceback;
    inline void setScore(double scr){
      score = scr;
      return;
    };
        
    //double path_prob (const HMM&, sequence&);  //Need to rewrite function
    inline int operator[](size_t val) const {return trace_path[val];};
    bool operator== (const TracebackPath&) const;
    bool operator<  (const TracebackPath&) const;
    bool operator>  (const TracebackPath&) const;
    bool operator<= (const TracebackPath&) const;
    bool operator>= (const TracebackPath&) const;
  private:
    Model* hmm;
    vector<int> trace_path;
    double score;
  };

  // ----------------------------------------------------------------------------------------
  //Rows are the Sequence Position
  //Columns are the States
  //! \def vector<vector<int> > heatTable;
  //! Table for heat map data generated from multiple tracebacks
  typedef vector<vector<int> > heatTable;

  //! \class multiTraceback
  //! Contains multiple tracebacks.  Will store them in sorted unique list (sorted in order of number of occurances);
  class multiTraceback{
  public:
    multiTraceback();
    ~multiTraceback();
        
    //Iteratorate through map;
    void begin();
    void end();
    void operator++();
    void operator--();
    void operator=(size_t);
        
        
    void print_path(string format="path", string *header=NULL);
    void print_label();
    void print_gff(string&);
    void print_hits();
        
        
    //Access the data at a point
    TracebackPath path();
    int counts();
    TracebackPath operator[](size_t);
        
        
    //Assign
    void assign(TracebackPath&);
        
    //Finalize multiTraceback (Sort and setup Iterators);
    void finalize();
            
    inline void clear(){paths.clear();return;};
    inline size_t size(){return paths.size();};
        
    heatTable* get_hit_table();

  private:
    size_t vectorIterator;
    size_t maxSize;
    bool isFinalized;
    vector<map<TracebackPath,int>::iterator> pathAccess;
    map<TracebackPath,int> paths;
    heatTable* table;
  };
     
  bool sortTBVec(map<TracebackPath,int>::iterator, map<TracebackPath,int>::iterator);
}
#endif



