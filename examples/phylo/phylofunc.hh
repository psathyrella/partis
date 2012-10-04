#ifndef __PHYLOFUNC_H__
#define __PHYLOFUNC_H__

#include "smctc.hh"
#include <memory>
#include <string>

class phylo_node
{
public:
	std::shared_ptr< phylo_node > child1;
	std::shared_ptr< phylo_node > child2;
	double dist1, dist2;
	double height;	// convenience for proposals, height must always increase
	std::string name;
};

class phylo_particle
{
public:
	std::shared_ptr< phylo_node > node;
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


#endif // __PHYLOFUNC_H__

