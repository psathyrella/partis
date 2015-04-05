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
  return exp(-0.5*factor*factor)) / (sigma * sqrt(2*3.1415926));  // normal distribution
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
smc::particle<State> Init()
{
  State init_state(all_uids);
  return smc::particle<State>(init_state, PartitionLogProb(init_state.partition_));
}

// choose a pair of clusters to merge according to the log probs
void Move(smc::particle<State> &partifrom, smc::rng *rgen)
{
  State *state = partifrom.GetValuePointer();

  // first find the net probability of all potential merges
  map<pair<unsigned, unsigned>, double> potential_merges;
  double total(0.0);
  for(unsigned ic1=0; ic1<state->partition_.size(), ++ic1) {
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
  for(auto &kv : potential_merges) {
    double netprob(kv.second);
    sum += netprob;
    if(sum > drawpoint) {
      icls = kv.first;
      break;
    }
  }
  assert(icls.first != 9999);

    cv_to->x_pos += cv_to->x_vel * Delta + pRng->Normal(0, sqrt(var_s));

    pFrom.AddToLogWeight(logLikelihood(lTime, *cv_to));
}

using namespace std;

///The observations
cv_obs * ydata;
long load_data(char const * szName, cv_obs** y);

double integrand_mean_x(const cv_state&, void*);
double integrand_mean_y(const cv_state&, void*);
double integrand_var_x(const cv_state&, void*);
double integrand_var_y(const cv_state&, void*);

int main(int argc, char** argv)
{
    long lNumber = 10;
    long lIterates;

    try {
        //Load observations
        lIterates = load_data("data.csv", &ydata);

        //Initialise and run the sampler
        smc::sampler<cv_state> Sampler(lNumber, SMC_HISTORY_NONE);
        smc::moveset<cv_state> Moveset(fInitialise, fMove);

        Sampler.SetResampleParams(SMC_RESAMPLE_RESIDUAL, 0.5);
        Sampler.SetMoveSet(Moveset);
        Sampler.Initialise();

        for(int n = 1 ; n < lIterates ; ++n) {
            Sampler.Iterate();

            double xm, xv, ym, yv;
            xm = Sampler.Integrate(integrand_mean_x, NULL);
            xv = Sampler.Integrate(integrand_var_x, (void*)&xm);
            ym = Sampler.Integrate(integrand_mean_y, NULL);
            yv = Sampler.Integrate(integrand_var_y, (void*)&ym);

            cout << xm << "," << ym << "," << xv << "," << yv << endl;
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        exit(e.lCode);
    }
}

long load_data(char const * szName, cv_obs** yp)
{
    FILE * fObs = fopen(szName, "rt");
    if(!fObs)
        throw SMC_EXCEPTION(SMCX_FILE_NOT_FOUND, "Error: pf assumes that the current directory contains an appropriate data file called data.csv\nThe first line should contain a constant indicating the number of data lines it contains.\nThe remaining lines should contain comma-separated pairs of x,y observations.");
    char* szBuffer = new char[1024];
    fgets(szBuffer, 1024, fObs);
    long lIterates = strtol(szBuffer, NULL, 10);

    *yp = new cv_obs[lIterates];

    for(long i = 0; i < lIterates; ++i) {
        fgets(szBuffer, 1024, fObs);
        (*yp)[i].x_pos = strtod(strtok(szBuffer, ",\r\n "), NULL);
        (*yp)[i].y_pos = strtod(strtok(NULL, ",\r\n "), NULL);
    }
    fclose(fObs);

    delete [] szBuffer;

    return lIterates;
}

double integrand_mean_x(const cv_state& s, void *)
{
    return s.x_pos;
}

double integrand_var_x(const cv_state& s, void* vmx)
{
    double* dmx = (double*)vmx;
    double d = (s.x_pos - (*dmx));
    return d * d;
}

double integrand_mean_y(const cv_state& s, void *)
{
    return s.y_pos;
}

double integrand_var_y(const cv_state& s, void* vmy)
{
    double* dmy = (double*)vmy;
    double d = (s.y_pos - (*dmy));
    return d * d;
}
