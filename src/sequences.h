//
//  sequences.h
//Copyright (c) 2007-2012 Paul C Lott 
//University of California, Davis
//Genome and Biomedical Sciences Facility
//UC Davis Genome Center
//Ian Korf Lab
//Website: www.korflab.ucdavis.edu
//Email: lottpaul@gmail.com
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of
//this software and associated documentation files (the "Software"), to deal in
//the Software without restriction, including without limitation the rights to
//use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
//the Software, and to permit persons to whom the Software is furnished to do so,
//subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
//FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
//COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
//IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef StochHMM_sequences_h
#define StochHMM_sequences_h

#include <vector>
#include <string>
#include <iostream>
#include <cassert>
#include <math.h>
#include <fstream>
#include <stdlib.h>
#include "text.h"
#include "track.h"
#include "externDefinitions.h"
#include "sequence.h"

namespace StochHMM {
class sequences{
public:
  sequences():length(0){}
  ~sequences();
      
  short seqValue(int , size_t); // Get digitized value for sequence track(i) in jth position
  sequence* getSeq(size_t); // Return ith sequence
  sequence& operator[](size_t index) {return *seqs[index]; }
  std::string& getHeader(size_t iter=0);
  std::string* getUndigitized(size_t iter) { return seqs[iter]->getUndigitized(); }
  size_t size() { return seqs.size(); }
  size_t getLength() { return length;}
  sequences getSubSequences(size_t pos, size_t len);
  bool exDefDefined() { return false; }  // removed this functionality
  bool exDefDefined(size_t) { return false; }  // removed this functionality
  
  std::string stringifyWOHeader();  //! Get string of sequences
  std::string stringify();  //! Get string of sequences
  void print() { std::cout << stringify() << std::endl; }
  std::string undigitize(); //! Get sequence based on alphabet
  
  void addSeq(sequence* sq);
  
private:
  std::vector<sequence*> seqs; 
  size_t length;  // length of the sequences (required to be the same for all)
};
}
#endif
