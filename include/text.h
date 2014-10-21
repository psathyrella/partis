#ifndef HAM_TEXT_H
#define HAM_TEXT_H

#include <sstream>
#include <string>
#include <vector>

using namespace std;
namespace ham {
  string Join(vector<int> &input, char c);
  string IntToString(int input);
  void ClearWhitespace(string white, string *input);
}
#endif
