#include "smctc.hh"
#include "phylofunc.hh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <memory>

using namespace std;

int main(int argc, char** argv)
{
	long population_size = 1000;
	long node_count = 10;
	long lIterates = node_count;

try {
	// TODO: read alignment
	// for now just generate some randomly named leaves
	leaf_nodes.resize(node_count);
	for(int i=0; i < node_count; i++){
		leaf_nodes[i] = make_shared< phylo_node >();
		stringstream ss;
		ss << "leaf_" << i;
		leaf_nodes[i]->name = ss.str();
	}

	//Initialise and run the sampler
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


