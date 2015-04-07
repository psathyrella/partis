#include <map>
#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "smctc.hh"

using namespace std;

vector<string> all_uids;
map<string, float> positions;
bool debug;

class State
{
public:
    State() {}
    State(vector<string> &uids)    // initialize with the all-separate partition
    {
        for(auto &uid : uids)
            partition_.push_back(vector<string> {uid});
    }
    vector<vector<string> > partition_;
};

// return probability of <x1> and <x2> being from the same cluster
double Prob(double x1, double x2)
{
    double sigma = 1.0;  // width
    double factor = (x1 - x2) / sigma;
    return exp(-0.5 * factor * factor) / (sigma * sqrt(2 * 3.1415926)); // normal distribution
}

// return log prob of <cluster>
double ClusterLogProb(vector<string> &cluster)
{
    // first find the mean x value of the cluster
    double xmean(0.0);
    for(auto &uid : cluster) {
        if(positions.count(uid) == 0)
            throw runtime_error("ERROR " + uid + " not in positions");
        xmean += positions[uid];
    }
    xmean /= cluster.size();

    // then add up the log probs of each point with respect to <xmean>
    double logprob(0.0);
    for(auto &uid : cluster)
        logprob += log(Prob(xmean, positions[uid]));

    return logprob;
}

// return the log prob of <state>'s partition
double PartitionLogProb(vector<vector<string> > partition)
{
    double logprob(0.0);
    for(auto &cluster : partition)
        logprob += ClusterLogProb(cluster);
    return logprob;
}

// ----------------------------------------------------------------------------------------
void PrintPartition(vector<vector<string> > partition)
{
    printf("   %7.2f:", PartitionLogProb(partition));
    for(unsigned ic = 0; ic < partition.size(); ++ic) {
        if(ic > 0)
            cout << "           ";
        for(auto &uid : partition[ic])
            cout << " " << uid;
        cout << endl;
    }
}

// return a randomized particle, i.e. a <State> with associated logprob
smc::particle<State> Init(smc::rng *pRng)
{
    State init_state(all_uids);
    return smc::particle<State>(init_state, PartitionLogProb(init_state.partition_));
}

// choose a pair of clusters to merge according to the log probs
void Move(long time, smc::particle<State> &partifrom, smc::rng *rgen)
{
    State *state = partifrom.GetValuePointer();
    if(state->partition_.size() == 1) {
        if(debug)
            cout << "nothing to move (only one cluster)" << endl;
        return;
    }

    // first find the net probability of all potential merges
    map<pair<unsigned, unsigned>, double> potential_merges;
    double total(0.0);
    if(debug) {
        cout << "moving" << endl;
        cout << "  potential" << endl;
    }
    for(unsigned ic1 = 0; ic1 < state->partition_.size(); ++ic1) {
        for(unsigned ic2 = ic1 + 1; ic2 < state->partition_.size(); ++ic2) {
            vector<string> cl1(state->partition_[ic1]), cl2(state->partition_[ic2]);
            vector<string> merged_cluster(cl1);
            merged_cluster.insert(merged_cluster.begin(), cl2.begin(), cl2.end());
            // change in the log probability of this partition that would result from merging these two clusters:
            double net_logprob = ClusterLogProb(merged_cluster) - ClusterLogProb(cl1) - ClusterLogProb(cl2);
            potential_merges[pair<unsigned, unsigned>(ic1, ic2)] = exp(net_logprob);  // should probably really stay in log space here
            total += exp(net_logprob);
            if(debug)
                printf("    %3d%3d     %.9f\n", ic1, ic2, exp(net_logprob));
        }
    }

    // normalize
    double chktotal(0.0);
    for(auto &kv : potential_merges) {
        kv.second /= total;
        chktotal += kv.second;
    }
    assert(fabs(chktotal - 1.) < 1e-7);

    // then choose one at random according to the probs
    double drawpoint = rgen->Uniform(0., 1.);
    double sum(0.0);
    pair<unsigned, unsigned> icls(9999, 9999);
    double chosennetprob(0.0);
    if(debug)
        cout << "  choosing with " << drawpoint << endl;
    for(auto &kv : potential_merges) {
        double netprob(kv.second);
        sum += netprob;
        if(debug)
            cout << "    " << kv.first.first << " " << kv.first.second << "    " << sum << endl;
        if(sum > drawpoint) {
            icls = kv.first;
            chosennetprob = netprob;
            break;
        }
    }
    assert(icls.first != 9999);
    assert(chosennetprob != 0.0);
    if(debug) {
        cout << "  chose " << icls.first << " " << icls.second << "    " << chosennetprob << endl;
        cout << "  before" << endl;
        PrintPartition(state->partition_);
    }

    state->partition_[icls.first].insert(state->partition_[icls.first].begin(),
                                         state->partition_[icls.second].begin(),
                                         state->partition_[icls.second].end());
    state->partition_.erase(state->partition_.begin() + icls.second);
    partifrom.AddToLogWeight(chosennetprob);

    if(debug) {
        cout << "  after" << endl;
        PrintPartition(state->partition_);
    }
}

// ----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
    debug = false;
    string infname("data.txt");
    ifstream ifs(infname);
    if(!ifs.is_open())
        throw runtime_error("ERROR input file " + infname + " d.n.e.");
    string line;
    int iline(0);
    while(getline(ifs, line)) {
        positions[to_string(iline)] = atof(line.c_str());
        ++iline;
    }
    for(auto &kv : positions)
        all_uids.push_back(kv.first);

    long n_particles(5);
    long n_max_steps(100);

    try {
        smc::sampler<State> smp(n_particles, SMC_HISTORY_RAM);
        smc::moveset<State> mvs(Init, Move);

        smp.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
        smp.SetMoveSet(mvs);
        smp.Initialise();

        vector<vector<State> > history;
        for(int ip = 0; ip < n_particles; ++ip) {
            history.push_back(vector<State>());
        }

        for(int istep = 1 ; istep < n_max_steps ; ++istep) {
            smp.Iterate();
            for(int ip = 0; ip < n_particles; ++ip) {
                history[ip].push_back(smp.GetParticleValue(ip));
            }

            bool all_finished(true);
            for(int ip = 0; ip < n_particles; ++ip) {
                bool this_finished = smp.GetParticleValue(ip).partition_.size() == 1;  // this one is finished if we've finished merging
                all_finished &= this_finished;
            }
            if(all_finished)
                break;
        }

        for(int ip = 0; ip < n_particles; ++ip) {
            cout << " particle " << ip << endl;
            for(int is = 0; is < n_max_steps - 1; ++is) { // not sure why the *@*% it has to be minus one
                if(history[ip][is].partition_.size() < 1)
                    break;
                PrintPartition(history[ip][is].partition_);
            }
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        exit(e.lCode);
    }
}
