#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
#include <cassert>

using namespace std;

// Normalizes counts per VDJ combination in the last column of a text
// file such that they sum to one. I.e. so you can throw probabilities from them.
// Runs on a file created with the following command:
// % bzgrep -v ',[1-9]$' 01-C-N_filtered.vdjcdr3.csv.bz2 | column -t -s, > filtered-vdj-counts.txt
// NOTE this removes all combinations with fewer than ten instances.

//----------------------------------------------------------------------------------------
void read_file(string fname, string config, double *total) {
  ifstream ifs(fname.c_str());
  assert(ifs.is_open());
  string line;
  if (config=="count") *total = 0.0;
  double check_sum(0.0);
  while (getline(ifs,line)) {
    if (line[0]=='#') continue;
    stringstream ss(line);
    string v_gene,d_gene,j_gene,cdr3_length;
    int count;
    double prob;
    if (config=="count") {
      ss >> v_gene >> d_gene >> j_gene >> cdr3_length >> count;
      *total += count;
    } else if (config=="normalize") {
      ss >> v_gene >> d_gene >> j_gene >> cdr3_length >> count;
      cout
	<< setw(22) << v_gene
	<< setw(22) << d_gene
	<< setw(22) << j_gene
	<< setw(8) << cdr3_length
	<< setw(18) << (double(count) / *total)
	<< endl;
    } else if (config=="check") {
      ss >> v_gene >> d_gene >> j_gene >> cdr3_length >> prob;
      check_sum += prob;
    } else {
      cout << "ERROR bad config" << endl;
      assert(0);
    }
  }

  if (config=="count")
    cout << "# total " << *total << endl;
  if (config=="check")
    cout << "# check_sum: " << check_sum << endl;
  ifs.close();
}
//----------------------------------------------------------------------------------------
int main() {
  double total(0.0);
  read_file("filtered-vdj-counts.txt", "count", &total);
  read_file("filtered-vdj-counts.txt", "normalize", &total);
  read_file("filtered-vdj-probs.txt", "check", &total);
}
