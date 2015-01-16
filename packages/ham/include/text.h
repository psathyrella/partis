#ifndef HAM_TEXT_H
#define HAM_TEXT_H

#include <sstream>
#include <string>
#include <vector>

using namespace std;
namespace ham {
void ClearWhitespace(string white, string *input);
vector<string> SplitString(string argstr, string separator);
}
#endif
