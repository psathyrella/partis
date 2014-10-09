#ifndef HAM_TRACEBACKPATH_H
#define HAM_TRACEBACKPATH_H

#include <iostream>
#include <fstream>

#include "text.h"
#include "model.h"
#include "mathutils.h"

using namespace std;

namespace ham {

class TracebackPath{
public:
  TracebackPath(Model* model) : hmm_(model) {}
  void push_back(int state) { path_.push_back(state); }
  void clear() { path_.clear(); }

  inline size_t size() const { return path_.size(); }
  inline Model* model() const { return hmm_; }
  inline double score() { return score_; }  // get score associated with this path
  void path_of_labels(vector<string> &str_path);  // put the path (as a vector of state labels) in <str_path>
  void path_of_labels(string &str_path); // or as a string
  void path_of_names(vector<string> &str_path);  // or use the names instead of the labels
  inline int operator[](size_t val) const {return path_[val];};
  bool operator== (const TracebackPath &rhs) const { return rhs.path_ == path_; }
  bool operator<  (const TracebackPath &rhs) const { return rhs.path_ < path_; }
  bool operator>  (const TracebackPath &rhs) const { return rhs.path_ > path_; }
  bool operator<= (const TracebackPath &rhs) const { return rhs.path_ <= path_; }
  bool operator>= (const TracebackPath &rhs) const { return rhs.path_ >= path_; }
  
  inline void set_model(Model* model) { hmm_ = model; }
  inline void set_score(double score) { score_ = score; }  // set score associated with this path
      
  void print_names() const ;
  void print_labels(string separator="") const ;  // if separator is specified, print it between each element in the sequence
private:
  Model* hmm_;
  vector<int> path_;
  double score_;
};

}
#endif
