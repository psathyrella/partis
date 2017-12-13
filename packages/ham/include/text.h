#ifndef HAM_TEXT_H
#define HAM_TEXT_H

#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <set>

using namespace std;
namespace ham {
void ClearWhitespace(string white, string *input);
vector<string> PythonSplit(string argstr);
vector<string> SplitString(string argstr, string delimiter=":");
bool InString(string query, string liststr, string delimiter=":");
string JoinStrings(vector<string> &strlist, string delimiter=":");
vector<int> Intify(vector<string> strlist);
vector<double> Floatify(vector<string> strlist);
}
#endif
