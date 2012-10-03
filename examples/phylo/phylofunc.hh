#ifndef __PHYLOFUNC_H__
#define __PHYLOFUNC_H__

#include "smctc.hh"

class phylo_node
{
public:
	phylo_node* child1;
	phylo_node* child2;
	double dist1, dist2;
	double height;	// convenience for proposals, height must always increase
	char* name;
};

class phylo_particle
{
public:
	phylo_node* node;
	phylo_particle* predecessor;
};

class particle
{
public:
	phylo_particle* pp;
};


double logLikelihood(long lTime, const particle & X);

smc::particle<particle> fInitialise(smc::rng *pRng);
long fSelect(long lTime, const smc::particle<particle>& p, smc::rng *pRng);
void fMove(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng);

extern std::vector< phylo_node* > leaf_nodes;


#endif // __PHYLOFUNC_H__

