#### Evaluating confidence in germline inference

We unfortunately do not yet have a way to report simple confidence estimates for each gene in the inferred germline set.
There is, however, a wealth of information that can be used to get a good sense of this confidence, and here we describe how to interpret it.
The following probably will not make much sense without some familiarity with the corresponding paper:

  * [Germline inference](https://arxiv.org/abs/1711.05843) Ralph, DK, & Matsen IV, FA (in review at PLOS Computational Biology) _Per-sample immunoglobulin germline inference \[...\]_ 

First off, you should run with a few extra arguments that tell partis to output more information than by default.
Set `--plotdir <dir>` to write the mutation accumulation fit plots to disk, as well as `--debug-allele-finding` to print extra info to std out.
It's probably worth piping the resulting std out to a log file, e.g. with `<cmd> >log.txt`, and viewing afterward with `less -RS` (typically looks nicer with a dark background terminal).

The first step in germline inference is allele removal, where we group all V genes with matches into subgroups such that genes within a group are separated by a handful of point mutations (i.e. confusable by SHM), but genes in different groups are not.
This process is illustrated in the first block of output printed by `--debug-allele-finding`, a screenshot of which is included here:

![what is this, xkcd?](images/allele-removal.png)
