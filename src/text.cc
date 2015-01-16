#include "text.h"
namespace ham {

// ----------------------------------------------------------------------------------------
void ClearWhitespace(string white, string *input) {
  size_t found = input->find_first_of(white);
  while(found != string::npos) {
    input->erase(found, 1);
    found = input->find_first_of(white);
  }
}

// ----------------------------------------------------------------------------------------
// split a colon-separated list in a string into a vector of strings, e.g. "a:b:c" --> {"a", "b", "c"}
vector<string> SplitString(string argstr, string separated) {
  vector<string> arglist;
  while(true) {
    size_t i_next_colon(argstr.find(":"));
    string arg = argstr.substr(0, i_next_colon); // get the next arg in the colon-separated list
    arglist.push_back(arg); // add it to arglist
    argstr = argstr.substr(i_next_colon + 1); // then excise it from argstr
    if(i_next_colon == string::npos)
      break;
  }
  return arglist;
}

}
