#Linking to StochHMM Static Library

##Under Construction

#Minimal Executable

```
#include <iostream>
#include <string>
#include "hmm.h"
#include "trellis.h"
#include "seqTracks.h"

int main(int argc, const char * argv[]){
    //Get Model and Sequence file arguments
    std::string model_file(argv[1]);
    std::string seq_file(argv[2]);

    //Create and Load the model
    model hmm;
    hmm.import(file,NULL);

    //Load the Sequences
    seqTracks jobs;
    jobs.loadSeqs(hmm,seq_file,FASTA);
    seqJob *job=jobs.getJob();

    //Initialize the trellis and perform viterbi algorithm
    trellis trell(hmm, job->getSeqs());    
    trell.viterbi();

    //Create traceback_path and traceback through trellis
    traceback_path pth(hmm);
    trell.traceback(pth);

    //Print the traceback in GFF format 
    pth.print_gff();

    return 0;
}
```