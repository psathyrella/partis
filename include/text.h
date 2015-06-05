#ifndef HAM_TEXT_H
#define HAM_TEXT_H

#include <sstream>
#include <string>
#include <vector>

using namespace std;
namespace ham {
void ClearWhitespace(string white, string *input);
vector<string> SplitString(string argstr, string delimiter=":");
string JoinStrings(vector<string> &strlist, string delimiter=":");
vector<int> Intify(vector<string> strlist);
vector<double> Floatify(vector<string> strlist);
}
#endif
