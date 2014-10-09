#include "tracebackpath.h"

namespace ham {

// ----------------------------------------------------------------------------------------
void TracebackPath::path_of_labels(vector<string> &str_path) {
  for(size_t k=path_.size()-1; k!=SIZE_MAX; --k) {
    State *st = hmm_->state(path_[k]);
    str_path.push_back(st->label());
  }
}

// ----------------------------------------------------------------------------------------
void TracebackPath::path_of_labels(string &str_path) {
  if (str_path.size() > 0)
    str_path.clear();
  for(size_t k=path_.size()-1; k!=SIZE_MAX; --k) {
    State* st = hmm_->state(path_[k]);
    str_path += st->label();
  }
}
  
// ----------------------------------------------------------------------------------------
void TracebackPath::path_of_names(vector<string> &str_path) {
  for (size_t k=path_.size()-1; k!=SIZE_MAX; --k) {
    State* st = hmm_->state(path_[k]);
    str_path.push_back(st->name());
  }
}
  
// ----------------------------------------------------------------------------------------
void TracebackPath::print_names() const {
  for(size_t k=path_.size()-1; k!=SIZE_MAX; --k) {
    State* st = hmm_->state(path_[k]);
    cout << st->name() << " ";
  }
  cout << endl;
}

// ----------------------------------------------------------------------------------------
void TracebackPath::print_labels(string separator) const {
  assert(hmm_);
  for(size_t k=path_.size()-1; k!=SIZE_MAX; --k) {
    State* st = hmm_->state(path_[k]);
    cout << st->label() << separator;
    cout.flush();
  }
  cout << endl;
}

}
