#### Evaluating confidence in germline inference

We unfortunately do not yet have a way to report simple confidence estimates for each gene in the inferred germline set.
There is, however, a wealth of information that can be used to get a good sense of this confidence, and here we describe how to interpret it.
The following probably will not make much sense without some familiarity with the corresponding paper:

  * [Germline inference](https://arxiv.org/abs/1711.05843) Ralph, DK, & Matsen IV, FA (in review at PLOS Computational Biology) _Per-sample immunoglobulin germline inference \[...\]_ 

First off, you should run with a few extra arguments that tell partis to output more information than by default.
Set `--plotdir <dir>` to write the mutation accumulation fit plots to disk, as well as `--debug-allele-finding` to print extra info to std out.
It's probably worth piping the resulting std out to a log file, e.g. with `<cmd> >log.txt`, and viewing afterward with `less -RS` (typically looks nicer with a dark background terminal).
The example log file that we'll describe in detail below can be found [here](images/example-log.txt).

The first step in germline inference is allele removal, where we group all V genes with matches into subgroups such that genes within a group are separated by a handful of point mutations (i.e. confusable by SHM), but genes in different groups are not.
This process is illustrated in the first block of output printed by `--debug-allele-finding`, a screenshot of which is included here:

![allele removal](images/allele-removal.png)

Each entry in the lefthand "genes to keep" column is a kept gene, i.e. the most common gene from that group, with the corresponding counts in the next column.
The alleles within that group that were removed are detailed in the last column.
This last column consists of two lists: a list of parenthesis-inclosed numbers for each removed gene (number of SNPs to the kept gene, and the counts), and the corresponding list of removed gene names.

After reducing the germline set to a minimal list of genes about whose presence we can be quite confident, partis then runs smith-waterman alignment on these genes, and performs new-allele inference.
This will re-infer any genes that were removed in the first step for which there is substantial evidence, along with any novel alleles (alleles not in the original germline set).

Because the mutation accumulation plots require independent mutation events, it then collapses clonal families into a small number of representatives.
The results are summarized in the debug output, an example of which is shown here:

![clonal collapse](images/clonal-collapse.png)

Where the first line, for instance, tells us that 2698 total sequences aligned most closely to IGHV1-46*01, that these corresponded to 1815 clonal families, and that 2034 sequences were chosen to represent these clonal families in germline inference.

Next we print a summary of the number of mutations observed in the sequences assigned to each gene in which the top line is the number of mutations:

![mutation summary](images/mutation-summary.png)

So, for instance, the top line tells us that we found 772 unmutated sequences assigned to IGHV1-46*01, 184 with one mutation, etc.

The next section documents the results of the method of excluding sequences that are missing too many bases one the 5' or 3' ends.
Details of this procedure can be found in the paper, and for most purposes this summary is probably not super interesting.

Next, it performs the actual germline inference fits for each gene.
Skipping a few sections to the one corresponding to IGHV1-69*13:

![1-69 numbers](images/1-69-numbers.png)

Within this block, there's a subsection for each candidate allele (number-of-SNPs) that was deemed interesting.
The first one (3 snps) has a row for each of the corresponding three positions (169, 219, 162; note: as everywhere, these are zero-indexed positions).
The primary measure of confidence in the inference is the ratio of the goodness of fit for the one piece vs two piece fits, which is the next column.
The first position (169) has a ratio of 26, which means the one piece fit was vastly worse than the two piece fit.
This is broken down in the next column (459. / 18.), which tells us that the one piece fit had a chi-square/dof of 459, while the two piece had 18.
Generally speaking, 18 is of course still a pretty dismal goodness of fit, even if it's vastly better than 459, so this is a good time to look at the corresponding plots.
The html summary for example is [here](http://psathyrella.github.io/partis/example-plots/germline-inference/try-0.html); each row corresponds to an inferred allele, and you can click on each plot to get the underlying svg.
Here's a screenshot of the row corresponding to the IGHV1-69*13 inference:

![1-69 fits](images/1-69-fits.png)
