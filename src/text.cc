#include "text.h"
namespace ham {

// ----------------------------------------------------------------------------------------
void ClearWhitespace(string white, string *input) {
  size_t found = input->find_first_of(white);
  while (found != string::npos) {
    input->erase(found, 1);
    found = input->find_first_of(white);
 }
}

// ----------------------------------------------------------------------------------------
string IntToString(int input) {
  stringstream ss;
  ss << input;
  return ss.str();
}

// ----------------------------------------------------------------------------------------
string Join(vector<int> &input, char c) {
  string out;
  if (input.size()==0){
    out="";
    return out;
  }
  else if (input.size()==1){
    out=IntToString(input[0]);
    return out;
  }
  else{
    out=IntToString(input[0]);
    for(size_t i=1;i<input.size();i++){
	out+=c;
	out+=IntToString(input[i]);
    }
    return out;
  }
}

}
