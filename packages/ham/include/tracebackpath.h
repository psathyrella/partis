#ifndef HAM_TRACEBACKPATH_H
#define HAM_TRACEBACKPATH_H

#include <iostream>
#include <fstream>

#include "text.h"
#include "model.h"
#include "mathutils.h"

using namespace std;

namespace ham {

class TracebackPath {
public:
  TracebackPath(Model* model) : hmm_(model), score_(0), abbreviate_(false) {}
  TracebackPath() : hmm_(nullptr) {}
  void push_back(int state) { path_.push_back(state); }
  void clear() { path_.clear(); }

  inline size_t size() const { return path_.size(); }
  inline void abbreviate(bool abb = true) { abbreviate_ = abb; }
  inline Model* model() const { return hmm_; }
  inline double score() { return score_; }  // get score associated with this path
  vector<string> name_vector();
  inline int operator[](size_t val) const {return path_[val];};
  bool operator== (const TracebackPath &rhs) const { return rhs.path_ == path_; }
  bool operator< (const TracebackPath &rhs) const { return rhs.path_ < path_; }
  bool operator> (const TracebackPath &rhs) const { return rhs.path_ > path_; }
  bool operator<= (const TracebackPath &rhs) const { return rhs.path_ <= path_; }
  bool operator>= (const TracebackPath &rhs) const { return rhs.path_ >= path_; }

  inline void set_model(Model* model) { hmm_ = model; }
  inline void set_score(double score) { score_ = score; }  // set score associated with this path

  friend ostream& operator<<(ostream &os, const TracebackPath &self) {
    for(size_t k = self.path_.size() - 1; k != SIZE_MAX; --k) {
      State* st = self.hmm_->state(self.path_[k]);
      os << (self.abbreviate_ ? st->abbreviation() : st->name()) << " ";
    }
    os << endl;
    return os;
  }

private:
  Model* hmm_;
  vector<int> path_;
  double score_;
  bool abbreviate_;
};

}
#endif
