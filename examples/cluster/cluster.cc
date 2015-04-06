#include <map>
#include <iostream>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "smctc.hh"

using namespace std;

vector<string> all_uids;
map<string, float> positions;

class State
{
public:
  State() {}
  State(vector<string> &uids) {  // initialize with the all-separate partition
    for(auto &uid : uids)
      partition_.push_back(vector<string>{uid});
  }
  vector<vector<string> > partition_;
};

// return probability of <x1> and <x2> being from the same cluster
double Prob(double x1, double x2) {
  double sigma = 1.0;  // width
  double factor = (x1-x2)/sigma;
  return exp(-0.5*factor*factor) / (sigma * sqrt(2*3.1415926));  // normal distribution
}

// return log prob of <cluster>
double ClusterLogProb(vector<string> &cluster) {
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

  // first find the net probability of all potential merges
  map<pair<unsigned, unsigned>, double> potential_merges;
  double total(0.0);
  for(unsigned ic1=0; ic1<state->partition_.size(); ++ic1) {
    for(unsigned ic2=ic1+1; ic2<state->partition_.size(); ++ic2) {
      vector<string> cl1(state->partition_[ic1]), cl2(state->partition_[ic2]);
      vector<string> merged_cluster(cl1);
      merged_cluster.insert(merged_cluster.begin(), cl2.begin(), cl2.end());
      // change in the log probability of this partition that would result from merging these two clusters:
      double net_logprob = ClusterLogProb(merged_cluster) - ClusterLogProb(cl1) - ClusterLogProb(cl2);
      potential_merges[pair<unsigned, unsigned>(ic1, ic2)] = exp(net_logprob);  // should probably really stay in log space here
      total += exp(net_logprob);
    }
  }

  // normalize
  double chktotal(0.0);
  for(auto &kv : potential_merges) {
    kv.second /= total;
    chktotal += kv.second;
  }
  cout << "TOTAL " << chktotal << endl;

  // then choose one at random according to the probs
  double drawpoint = rgen->Uniform(0., 1.);
  double sum(0.0);
  pair<unsigned, unsigned> icls(9999, 9999);
  double chosennetprob(0.0);
  for(auto &kv : potential_merges) {
    double netprob(kv.second);
    sum += netprob;
    if(sum > drawpoint) {
      icls = kv.first;
      chosennetprob = netprob;
      break;
    }
  }
  assert(icls.first != 9999);
  assert(chosennetprob != 0.0);

  // vector<string> cl1(state->partition_[icls.first]), cl2(state->partition_[icls.second]);
  // vector<string> merged_cluster(cl1);
  // merged_cluster.insert(merged_cluster.begin(), cl2.begin(), cl2.end());
  state->partition_[icls.first].insert(state->partition_[icls.first].begin(),
				       state->partition_[icls.second].begin(),
				       state->partition_[icls.second].end());
  state->partition_.erase(state->partition_.begin() + icls.second);
  partifrom.AddToLogWeight(chosennetprob);
}

// ----------------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  long n_particles(10);
  long n_steps(10);

  try {
    smc::sampler<State> smp(n_particles, SMC_HISTORY_NONE);
    smc::moveset<State> mvs(Init, Move);

    smp.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
    smp.SetMoveSet(mvs);
    smp.Initialise();

    for(int n = 1 ; n < n_steps ; ++n) {
      smp.Iterate();
      State state(smp.GetParticleValue(0));
      cout << state.partition_.size() << endl;
    }
  }

  catch(smc::exception  e) {
    cerr << e;
    exit(e.lCode);
  }
}
