#ifndef __hmsbeagle__
#define __hmsbeagle__

#include <string>
#include <cstdlib>
#include <cstdio>
#include "libhmsbeagle/beagle.h"
#include "phylofunc.hh"
#include <memory>
#include <vector>
#include <stack>

double* get_partials(const std::string& sequence);

class OnlineCalculator
{
public:
    OnlineCalculator() : stateCount(4), instance(-1), next_id(0) {};
    ~OnlineCalculator() {
        beagleFinalizeInstance(instance);
    };

    double calculate_ll( std::shared_ptr< phylo_node > node, std::vector<bool>& visited );
    void initialize(const std::vector<std::string>& seqs);
    int get_id();
    void free_id(int id);
private:
    BeagleInstanceDetails instDetails;
    int stateCount;
    int nPatterns;
    int nPartBuffs;
    int instance;
    int treeknown;

    int next_id;
    std::stack<int> free_ids;
};

int OnlineCalculator::get_id()
{
    if(free_ids.size()>0) {
        int id = free_ids.top();
        free_ids.pop();
        return id;
    }
    if(next_id == nPartBuffs) {
        // oh no! we ran out of buffer slots! need to reallocate
        throw "Ran out of beagle slots";
    }
    return next_id++;
}

void OnlineCalculator::free_id(int id)
{
    free_ids.push(id);
}

void OnlineCalculator::initialize(const std::vector<std::string>& seqs)
{
    // initialize the instance
    bool reinit = false;
    if(nPatterns != seqs[1].size()) reinit = true;
    if(seqs.size()*2 + 2 >= nPartBuffs) reinit = true;
    if(instance < 0) reinit = true;


    // create an instance of the BEAGLE library
    if( reinit ) {
        nPatterns = seqs[1].size();
        nPartBuffs = seqs.size() * 3000;
        if(instance >= 0) beagleFinalizeInstance(instance);
        instance = beagleCreateInstance(
                       0,			/**< Number of tip data elements (input) */
                       nPartBuffs,       	/**< Number of partials buffers to create (input) */
                       0,		        /**< Number of compact state representation buffers to create (input) */
                       stateCount,		/**< Number of states in the continuous-time Markov chain (input) */
                       nPatterns,		/**< Number of site patterns to be handled by the instance (input) */
                       nPartBuffs,	        /**< Number of rate matrix eigen-decomposition buffers to allocate (input) */
                       nPartBuffs,	        /**< Number of rate matrix buffers (input) */
                       1,			/**< Number of rate categories (input) */
                       nPartBuffs,		/**< Number of scaling buffers */
                       NULL,			/**< List of potential resource on which this instance is allowed (input, NULL implies no restriction */
                       0,			    /**< Length of resourceList list (input) */
                       BEAGLE_FLAG_PRECISION_DOUBLE | BEAGLE_FLAG_SCALING_AUTO,	/**< Bit-flags indicating preferred implementation charactertistics, see BeagleFlags (input) */
                       0,                /**< Bit-flags indicating required implementation characteristics, see BeagleFlags (input) */
                       &instDetails);
        if (instance < 0) {
            fprintf(stderr, "Failed to obtain BEAGLE instance\n\n");
            exit(1);
        }

        treeknown = 0;	// none of the tip sequences have been set
    }


    // add any new sequences to the beagle instance
    for(int i=treeknown; i<seqs.size(); i++) {
        double *seq_partials   = get_partials(seqs[i]);
        beagleSetPartials(instance, i, seq_partials);
//		std::cout << "Adding new sequence " << i << "\n" << seqs[i] << std::endl;
    }
    treeknown = seqs.size();
    next_id = seqs.size();

    if( reinit ) {
        // create base frequency array
        double freqs[16] = { 0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25,
                             0.25, 0.25, 0.25, 0.25
                           };

        // an eigen decomposition for the JC69 model
        double evec[4 * 4] = {
            1.0,  2.0,  0.0,  0.5,
            1.0,  -2.0,  0.5,  0.0,
            1.0,  2.0, 0.0,  -0.5,
            1.0,  -2.0,  -0.5,  0.0
        };

        double ivec[4 * 4] = {
            0.25,  0.25,  0.25,  0.25,
            0.125,  -0.125,  0.125,  -0.125,
            0.0,  1.0,  0.0,  -1.0,
            1.0,  0.0,  -1.0,  0.0
        };

        double eval[4] = { 0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333 };

        // set the Eigen decomposition
        beagleSetEigenDecomposition(instance, 0, evec, ivec, eval);

        beagleSetStateFrequencies(instance, 0, freqs);
        double rates[1] = { 1 };
        beagleSetCategoryRates(instance, &rates[0]);
        double weights[1] = { 1 };
        beagleSetCategoryWeights(instance, 0, weights);

        double* patternWeights = (double*) malloc(sizeof(double) * nPatterns);
        for (int i = 0; i < nPatterns; i++) {
            patternWeights[i] = 1.0;
        }
        beagleSetPatternWeights(instance, patternWeights);
    }
}

double OnlineCalculator::calculate_ll( std::shared_ptr< phylo_node > node, std::vector<bool>& visited )
{
    // prepare visitor vector
    if(visited.size() < nPartBuffs) {
        visited.resize(nPartBuffs);
    }

    // accumulate a vector of operations through depth first search
    std::vector< BeagleOperation > fops, ops;
    std::vector< int > nind;
    std::vector< double > lens;
    std::stack< std::shared_ptr< phylo_node > > s;
    s.push(node);
    visited[node->id]=true;
    while(s.size()>0) {
        std::shared_ptr< phylo_node > cur = s.top();
        s.pop();
        if(cur->child1 == NULL)
            continue;
        fops.push_back( {cur->id, BEAGLE_OP_NONE, BEAGLE_OP_NONE, cur->child1->id, cur->child1->id, cur->child2->id, cur->child2->id} );
        nind.push_back(cur->child1->id);
        nind.push_back(cur->child2->id);
        lens.push_back(cur->dist1);
        lens.push_back(cur->dist2);
        visited[cur->child1->id]=true;
        visited[cur->child2->id]=true;
    }
    // need to reverse the order to make post-order
    ops.insert(ops.begin(), fops.rbegin(),fops.rend());

    /*
    	int nodeIndices[tree.size()];
    	double edgeLengths[tree.size()];
    	for(int i=0; i<tree.size(); i++){
    		nodeIndices[i]=i;
    		edgeLengths[i]=tree[i].distance;
    	}
    */

    // tell BEAGLE to populate the transition matrices for the above edge lengths
    beagleUpdateTransitionMatrices(instance,     // instance
                                   0,             // eigenIndex
                                   nind.data(),   // probabilityIndices
                                   NULL,          // firstDerivativeIndices
                                   NULL,          // secondDerivativeIndices
                                   lens.data(),   // edgeLengths
                                   nind.size());            // count

    // create a list of partial likelihood update operations
    // the order is [dest, destScaling, sourceScaling, source1, matrix1, source2, matrix2]
    // TODO: make this peel only where the new node was added.
    BeagleOperation operations[ops.size()];
    int scaleIndices[ops.size()];
    for(int i=0; i<ops.size(); i++) {
        scaleIndices[i]=ops[i].destinationPartials;
        operations[i] = ops[i];
    }

    // update the partials
    beagleUpdatePartials(instance, operations, ops.size(), node->id);// cumulative scaling index
    beagleAccumulateScaleFactors(instance, scaleIndices, ops.size(), BEAGLE_OP_NONE);

    double logL = 0.0;
    int returnCode = 0;

    // calculate the site likelihoods at the root node
    int rootIndices[ 1 ] = { node->id };
    int categoryWeightsIndices[ 1 ] = { 0 };
    int stateFrequencyIndices[ 1 ] = { 0 };
    int cumulativeScalingIndices[ 1 ] = { BEAGLE_OP_NONE };
    returnCode = beagleCalculateRootLogLikelihoods(instance,               // instance
                 (const int *)rootIndices,// bufferIndices
                 (const int *)categoryWeightsIndices,                // weights
                 (const int *)stateFrequencyIndices,                  // stateFrequencies
                 cumulativeScalingIndices,// cumulative scaling index
                 1,                      // count
                 &logL);         // outLogLikelihoods


    return logL;
}




double* get_partials(const std::string& sequence)
{
    int n = sequence.size();
    double *partials = (double*)malloc(sizeof(double) * n * 4);

    int k = 0;
    for (int i = 0; i < n; i++) {
        switch (sequence[i]) {
        case 'A':
            partials[k++] = 1;
            partials[k++] = 0;
            partials[k++] = 0;
            partials[k++] = 0;
            break;
        case 'C':
            partials[k++] = 0;
            partials[k++] = 1;
            partials[k++] = 0;
            partials[k++] = 0;
            break;
        case 'G':
            partials[k++] = 0;
            partials[k++] = 0;
            partials[k++] = 1;
            partials[k++] = 0;
            break;
        case 'T':
            partials[k++] = 0;
            partials[k++] = 0;
            partials[k++] = 0;
            partials[k++] = 1;
            break;
        default:
            partials[k++] = 1;
            partials[k++] = 1;
            partials[k++] = 1;
            partials[k++] = 1;
            break;
        }
    }
    return partials;
}



#endif //  __hmsbeagle__


