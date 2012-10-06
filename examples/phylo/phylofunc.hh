#ifndef __PHYLOFUNC_H__
#define __PHYLOFUNC_H__

#include "smctc.hh"
#include <memory>
#include <string>

// This class stores the SMC forest implicitly, by specifying the collections
// of mergers that must be made in order to get the forest from \perp (i.e. the
// completely un-merged state).

class phylo_node
{
public:
	phylo_node();
	~phylo_node();

	std::shared_ptr< phylo_node > child1;
	std::shared_ptr< phylo_node > child2;
	double dist1, dist2;
	double height;	// convenience for proposals, height must always increase
	int id;	// node id (1..n-1) for leaf nodes, corresponds to index in alignment. n..2n-1 for internal nodes.
};

class phylo_particle
{
public:
	// The merge novel to this particle. If NULL then it's \perp.
	std::shared_ptr< phylo_node > node;
	// The predecessor particles, which specify the rest of the merges for this particle.
	std::shared_ptr< phylo_particle > predecessor;
};

class particle
{
public:
	std::shared_ptr< phylo_particle > pp;
};


double logLikelihood(long lTime, const particle & X);

smc::particle<particle> fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<particle>& p, smc::rng *pRng);
void fMove(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng);

extern std::vector< std::shared_ptr< phylo_node > > leaf_nodes;
extern std::vector< std::pair< std::string, std::string > > aln;


#endif // __PHYLOFUNC_H__

