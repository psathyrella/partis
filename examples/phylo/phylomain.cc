#include "smctc.hh"
#include "phylofunc.hh"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <fstream>
#include <stack>
#include <memory>

using namespace std;

void read_alignment( istream& in, vector< pair< string, string > >& aln )
{
    string line, cur_seq, name;
    while(getline(in, line)) {
        if(line[0] == '>') {
            if(cur_seq.size()>0) {
                aln.push_back( make_pair( name, cur_seq ) );
            }
            name = line.substr(1);
            cur_seq = "";
        } else {
            cur_seq += line;
        }
    }
    aln.push_back( make_pair( name, cur_seq ) );
}

bool check_visited( vector< bool >& visited, int id )
{
    // ensure visited has enough space allocated to store the id
    // if not, resize it large enough and leave some wiggle to prevent frequent resizings
    if(id >= visited.size()){
        visited.resize(id+100);
    }
}

bool visited_id( vector< bool >& visited, int id )
{
    check_visited( visited, id );
    return visited[id];
}

bool set_visited_id( vector< bool >& visited, int id )
{
    check_visited( visited, id );
    visited[id]=true;
}

void write_tree( ostream& out, shared_ptr< phylo_node > root )
{
    vector< bool > visited;
    stack< shared_ptr< phylo_node > > s;
    s.push(root);
    while(s.size() > 0) {
        shared_ptr< phylo_node > cur = s.top();
        if(cur->child1 == NULL) {
            cout << aln[cur->id].first;
            set_visited_id( visited, cur->id );
            s.pop();
            continue;
        }
        if(!visited_id( visited, cur->child1->id )) {
            cout << "(";
            s.push(cur->child1);
            continue;
        } else if(!visited_id( visited, cur->child2->id )) {
            cout << ":" << cur->dist1 << ",";
            s.push(cur->child2);
            continue;
        }
        cout << ":" << cur->dist2 << ")";
        set_visited_id( visited, cur->id );
        s.pop();
    }
    cout << ";\n";
}

int main(int argc, char** argv)
{
    if(argc != 2) {
        cerr << "Usage: phylo <fasta alignment>\n\n";
        return -1;
    }
    long population_size = 1000;

    string file_name = argv[1];
    ifstream in(file_name.c_str());
    read_alignment( in, aln );

    long lIterates = aln.size();

    try {
        leaf_nodes.resize(aln.size());
        for(int i=0; i < aln.size(); i++) {
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

            double max_ll = -std::numeric_limits<double>::max();
            for(int i=0; i < population_size; i++) {
                particle X = Sampler.GetParticleValue(i);
                // write the log likelihood
                double ll = logLikelihood( lIterates, X );
                max_ll = max_ll > ll ? max_ll : ll;
            }
            cerr << "Iter " << n << " max ll " << max_ll << endl;
        }

        for(int i=0; i < population_size; i++) {
            particle X = Sampler.GetParticleValue(i);
            // write the log likelihood
            cout << logLikelihood( lIterates, X ) << "\t";
            // write out the tree under this particle
            write_tree( cout, X.pp->node );
        }
    }

    catch(smc::exception  e) {
        cerr << e;
        exit(e.lCode);
    }
}


