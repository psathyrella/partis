#include "text.h"
namespace ham {

// ----------------------------------------------------------------------------------------
// removes all instances of string <white> (nothing particular to white space)
void ClearWhitespace(string white, string *input) {
  size_t found = input->find_first_of(white);
  while(found != string::npos) {
    input->erase(found, 1);
    found = input->find_first_of(white);
  }
}

// ----------------------------------------------------------------------------------------
// split <argstr> into a vector by white space, as in python's str.split() NOTE not the same as SplitString() below
vector<string> PythonSplit(string argstr) {
  istringstream iss(argstr);
  vector<string> tokens{istream_iterator<string>{iss}, istream_iterator<string>{}};
  return tokens;
}

// ----------------------------------------------------------------------------------------
// split a <delimiter>-separated list in a string into a vector of strings, e.g. "a:b:c" --> {"a", "b", "c"}
vector<string> SplitString(string argstr, string delimiter) {
  vector<string> arglist;
  while(true) {
    size_t i_next_colon(argstr.find(delimiter));  // not necessarily colon, but it does default to that
    string arg = argstr.substr(0, i_next_colon); // get the next arg in the <delimiter>-separated list
    arglist.push_back(arg); // add it to arglist
    argstr = argstr.substr(i_next_colon + 1); // then excise it from argstr
    if(i_next_colon == string::npos)
      break;
  }
  return arglist;
}

// ----------------------------------------------------------------------------------------
// find if <query> is in the colon-separated list <liststr>
bool InString(string query, string liststr, string delimiter) {
  while(true) {
    size_t i_next_colon(liststr.find(delimiter));  // not necessarily colon, but it does default to that
    string arg = liststr.substr(0, i_next_colon); // get the next arg in the <delimiter>-separated list
    if(arg == query)
      return true;
    liststr = liststr.substr(i_next_colon + 1); // then excise it from liststr
    if(i_next_colon == string::npos)
      break;
  }
  return false;  // didn't find it
}

// ----------------------------------------------------------------------------------------
string JoinStrings(vector<string> &strlist, string delimiter) {
  string return_str;
  for(size_t is=0; is<strlist.size(); ++is) {
    if(is > 0)
      return_str += delimiter;
    return_str += strlist[is];
  }
  return return_str;
}

// ----------------------------------------------------------------------------------------
vector<int> Intify(vector<string> strlist) {
  vector<int> intlist;
  for(auto &str : strlist)
    intlist.push_back(stoi(str));
  return intlist;
}

// ----------------------------------------------------------------------------------------
vector<double> Floatify(vector<string> strlist) {
  vector<double> floatlist;
  for(auto &str : strlist)
    floatlist.push_back(stof(str));
  return floatlist;
}

}
