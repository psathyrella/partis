#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <vector>
#include <stack>
#include <unordered_set>
#include <algorithm>
#include "smctc.hh"
#include "phylofunc.hh"

using namespace std;

vector< phylo_node* > leaf_nodes;

///The function corresponding to the log likelihood at specified time and position (up to normalisation)

///  \param lTime The current time (i.e. the number of coalescence events so far)
///  \param X     The state to consider
double logLikelihood(long lTime, const particle& X)
{
	// TODO: implement phylo likelihood!
	return 1.0;
}

///A function to initialise particles

/// \param pRng A pointer to the random number generator which is to be used
smc::particle<particle> fInitialise(smc::rng *pRng)
{
	particle value;
	// initial particles have all sequences uncoalesced
	value.pp = new phylo_particle();
	// loglike should just be the background distribution on character state frequencies
	return smc::particle<particle>(value,logLikelihood(0,value));
}

///The proposal function.

///\param lTime The sampler iteration.
///\param pFrom The particle to move.
///\param pRng  A random number generator.
void fMove(long lTime, smc::particle<particle>& pFrom, smc::rng *pRng)
{
	particle* part = pFrom.GetValuePointer();

	phylo_particle* pp = new phylo_particle();
	pp->predecessor = part->pp;
	part->pp = pp;
	// Propose to merge two nodes.
	// Accumulate a list of all possible nodes for merging, add all
	// possible nodes to set, then remove those which are not active.
	unordered_set<phylo_node*> proposal_set;
        // Insert all of the leaf nodes into the proposal set.
	proposal_set.insert( leaf_nodes.begin(), leaf_nodes.end() );
	// Walk back to predecessor particles, adding root nodes and removing
	// coalesced nodes.
	unordered_set<phylo_node*> coalesced;
	for( phylo_particle* cur = pp->predecessor; cur != NULL; cur = cur->predecessor ){
		if(cur->node == NULL) continue;
		// Skip if we've already processed this subtree.
		if(coalesced.find(cur->node) != coalesced.end()) continue;
		// Insert this active root node to the proposal set.
		proposal_set.insert(cur->node);
		stack<phylo_node*> s;
		// Recursively add all descendants of the root nodes to the coalesced set.
		s.push(cur->node);
		while(s.size()>0){
			phylo_node* n = s.top();
			s.pop();
			if(n->child1 == NULL) continue;	// leaf node, nothing more to do.
			coalesced.insert(n->child1);
			coalesced.insert(n->child2);
			s.push(n->child1);
			s.push(n->child2);
		}
	}

	// The set difference of available (i.e. proposal_set) and coalesced
	// nodes yields the final proposal set; store it in prop_vector.
	vector< phylo_node* > prop_vector(proposal_set.size()+coalesced.size());
	vector< phylo_node* >::iterator last_ins = set_difference( proposal_set.begin(), proposal_set.end(), coalesced.begin(), coalesced.end(), prop_vector.begin() );
	prop_vector.resize( last_ins - prop_vector.begin() );

	// Pick two nodes from the prop_vector to join.
	int n1 = pRng->UniformDiscrete(0, prop_vector.size()-1);
	int n2 = pRng->UniformDiscrete(0, prop_vector.size()-1);
	pp->node = new phylo_node();
	pp->node->child1 = prop_vector[n1];
	pp->node->child2 = prop_vector[n2];
	// Propose a coalescence time.
	double h = pRng->Uniform(0, 2);
	double prev_h = pp->predecessor->node != NULL ? pp->predecessor->node->height : 0;
	pp->node->height = prev_h + h;
	pp->node->dist1 = h - pp->node->child1->height;
	pp->node->dist2 = h - pp->node->child2->height;

	pFrom.AddToLogWeight(logLikelihood(lTime, *part));
}


