#include "tracebackpath.h"

namespace ham {

// ----------------------------------------------------------------------------------------
vector<string> TracebackPath::name_vector() {
  vector<string> str_path;
  for(size_t k = path_.size() - 1; k != SIZE_MAX; --k) {
    State *st = hmm_->state(path_[k]);
    str_path.push_back(st->name());
  }
  return str_path;
}

}
