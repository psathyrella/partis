#include "smctc.hh"
#include "phylofunc.hh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <fstream>
#include <memory>

using namespace std;

void read_alignment( istream& in, vector< pair< string, string > >& aln ){
	string line, cur_seq, name;
	while(getline(in, line)){
		if(line[0] == '>'){
			if(cur_seq.size()>0){
				aln.push_back( make_pair( name, cur_seq ) );
			}
			name = line.substr(1);
			cur_seq = "";
		}else{
			cur_seq += line;
		}
	}
	aln.push_back( make_pair( name, cur_seq ) );
}

int main(int argc, char** argv)
{
	if(argc != 2){
		cerr << "Usage: phylo <fasta alignment>\n\n";
		return -1;
	}
	long population_size = 1000;
	long node_count = 10;

	string file_name = argv[1];
	ifstream in(file_name.c_str());
	read_alignment( in, aln );

	long lIterates = node_count;

try {
	// TODO: read alignment
	// For now just generate some simply named leaves in the
	// globally-declared leaf_nodes vector.
	leaf_nodes.resize(aln.size());
	for(int i=0; i < node_count; i++){
		leaf_nodes[i] = make_shared< phylo_node >();
		leaf_nodes[i]->id = i;
	}

	// Initialise and run the sampler
	smc::sampler<particle> Sampler(population_size, SMC_HISTORY_NONE);
	smc::moveset<particle> Moveset(fInitialise, fMove, NULL);

	Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
	Sampler.SetMoveSet(Moveset);
	Sampler.Initialise();

	for(int n=1 ; n < lIterates ; ++n) {
		Sampler.Iterate();
	}
}

catch(smc::exception  e)
{
	cerr << e;
	exit(e.lCode);
}
}


