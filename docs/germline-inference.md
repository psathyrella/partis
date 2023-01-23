[Up to table of contents](contents.md)

#### Evaluating confidence in germline inference

We unfortunately do not yet have a way to report simple confidence estimates for each gene in the inferred germline set.
There is, however, a wealth of information that can be used to get a good sense of this confidence, and here we describe how to interpret it.
The following probably will not make much sense without some familiarity with the corresponding paper:

  * [Germline inference](https://doi.org/10.1371/journal.pcbi.1007133) Ralph, DK, & Matsen IV, FA (2019). _Per-sample immunoglobulin germline inference \[...\]._ PLOS Computational Biology, 15(7), e1007133.

First off, you should run with a few extra arguments that tell partis to output more information than it does by default.
Set `--plotdir <dir>` to write the mutation accumulation fit plots to disk, as well as `--debug-allele-finding` to print extra info to std out.
It's probably worth piping the resulting std out to a log file, e.g. with `<cmd> >log.txt`, and viewing it afterward with `less -RS` (it typically looks nicer with a dark background terminal).
The example log file that we'll describe in detail below can be found [here](images/example-log.txt).

The first step in germline inference is allele removal, where we group all V genes with matches into subgroups such that genes within a group are separated by a handful of point mutations (i.e. confusable by SHM), but genes in different groups are not.
This process is illustrated in the first block of output printed by `--debug-allele-finding`, a screenshot of which is included here:

![allele removal](images/allele-removal.png)

Each entry in the lefthand "genes to keep" column is a kept gene, i.e. the most common gene from that group, with the corresponding counts in the next column.
The alleles within that group that were removed are detailed in the last column.
This last column consists of two lists: a list of parenthesis-enclosed numbers for each removed gene (number of SNPs from the kept gene, and the counts), and the corresponding list of removed gene names.
So for instance, the first line tells us that IGHV4-39*07 was the most common gene in this group, with 34680 sequences, and there were four less common genes in the group that were removed.
The most common of these removed genes was IGHV4-39*01, which differs from *07 by three SNPs, and was observed 174 times.

After reducing the germline set to a minimal set of genes about whose presence we can be quite confident, partis then runs smith-waterman alignment on these genes, and performs new-allele inference.
This will re-infer any genes that were removed in the first step for which there is substantial evidence, along with any novel alleles.

Because the mutation accumulation plots require independent mutation events, it first collapses clonal families into a small number of representatives.
The results are summarized in the next bit of the debug output:

![clonal collapse](images/clonal-collapse.png)

Where the first line, for instance, tells us that 2698 total sequences aligned most closely to IGHV1-46*01, that these corresponded to 1815 clonal families, and that 2034 sequences were chosen to represent these clonal families in germline inference.

Next we print a summary of the number of mutations observed in the sequences assigned to each gene in which the top line is the number of mutations:

![mutation summary](images/mutation-summary.png)

So, for instance, the first gene line tells us that we found 772 unmutated sequences assigned to IGHV1-46*01, 184 with one mutation, etc.

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
The html summary for this example is [here](http://psathyrella.github.io/partis/example-plots/germline-inference/try-0.html); each row corresponds to an inferred allele, and you can click on each plot to get the underlying svg.
Here's a screenshot of the row corresponding to the IGHV1-69*13 inference:

![1-69 fits](images/1-69-fits.png)

These are the same three positions as in the previous paragraph, and each bin in these plots corresponds to a column in the table.
For instance, the fourth column in the table for position 169 says that of the 1892 unmutated sequences, zero were mutated at position 169.
Somewhat less tautologically, the next column says that of the 1269 sequences with one V mutation, two had that mutation at position 169.
These are exactly the numbers that go into the bins in the fit plots.

The fit plots show the one-piece (null hypothesis) in green, and the two-piece (new allele hypothesis) in red.
And here we see that, indeed, the red fits much better than the green, but is still a pretty crappy fit.
In this case, the poor fit of the two-piece fit is actually expected: this sample happens to have a very low SHM rate, with a big spike near zero.
The SHM distributions are written under your plotdir at `<plotdir>/sw/mute-freqs/overall.html`, and are screenshot here:

![SHM rate](images/shm-rate.png)

This spike in overall SHM rate creates a spike in the 3-bin in these fits, and is the reason that we use the ratio of goodness of fit values.
The're many more details in the paper, but in short everything is working as it should.

Taking it all together, this is a very strong inference.
The ratio of goodness of fit values for one piece/two piece is very large, and while it's not shown in the paper, this is distributed as essentially a normal distribution.
This means that it corresponds roughly to p-values, i.e. a value between two and three corresponds to between roughly .1 and 1 percent, hence why the thresholds for calling new genes are around that range (see paper).

Continuing down, the next block in the IGHV1-69\*13 table above is for a 7-SNP candidate that is not significant.
Then, it says that the minimum ratio among positions for the 3-SNP candidate was 20, which is big, so it infers the allele.
It then checks through known genes, and discovers that this is a known gene (that, indeed, was removed in the first step): IGHV1-69\*10.
Note that removing and re-inferring it is much more accurate than try to guess a prevalence threshold in the allele removal step (many more details in the paper).

The rest of the log file has similar information for the rest of the inferences, and summary information at the end.
But for the moment, it's quicker to glance through the rest of the plots in the html file linked above.
The second row shows a slam-dunk inference of an allele with prevalence of 1 (i.e. the "template" gene isn't actually present):

![homozygous fit](images/homozygous-fit.png)

Note that since the y axis is the fraction of times that the position of interest was mutated, the height of the red line corresponds roughly to the prevalence of the inferred allele.
In cases such as the IGHV1-69 inference above where the red line is very far from horizontal, this relationship is more complicated.

Scrolling down further, the two inferences from a template of IGHV4-59*01 are interesting:

![4-59 fits](images/4-59-fits.png)

They're both pretty solid fits in that the two-piece red fit describes the data well, but they're quite different in prevalence.
The bottom one (position 260) has much higher prevalence (30-40%), and thus enough observations for very confident inference (cross-referencing to the log file, it has a one-piece/two-piece ratio of 10.
The top one, though (position 114), has a prevalence of only 1-2%, and as seen by the uncertainties in the fit plot, it is far less certain.
Looking at the log file:

![4-59 table](images/4-59-table.png)

this has a high ratio of 18, but this is entirely due to a very low goodness of fit of 0.15, which indicates an unrealistically good fit, which is almost always due to chance (well, or incorrect uncertainties. But we checked our uncertainty model pretty carefully).
In the plot, you can see this -- the points are scattered much less than you'd expect based on the uncertainties.
More importantly, the green null hypothesis doesn't actually fit that badly (goodness of fit of 2.72), indicating decent, but not overwhelming, confidence.
Another thing to notice in the table is that since it tells you the exact numbers that go into the fits, you can make judgments on the data independent of our uncertainty and fitting models.
You can for instance see that the position 260 inference is essentially based on 506 sequences being in the 1-mutation bin, rather than (mostly) in the 0-bin.
That is a lot of observations.
Whereas the position 114 inference is mostly based on 25 sequences being similarly out of place: observing something 20-odd times is enough to be pretty darn suspicious, but it's not five sigma no matter your model.
Keep in mind that these numbers are after clonal collapse, though, so a more accurate clustering algorithm (see the paper, but we don't use the most accurate one for this) might get different numbers.

Finally, the prevalence frequencies for each gene in the inferred germline set after the smith-waterman step are also worth looking at.
These are in the `--plotdir`, under `<plotdir>/sw/overall/{v,d,j}_gene.{csv,svg}`, with the svg version screenshotted here:

![v gene counts](images/v-gene-counts.png)

but there is also a csv version if you want to read for further processing.
These are also written after the hmm step (in `<plotdir>/hmm/overall/{v,d,j}_gene.{csv,svg}`), and those will generally be very similar, except that a more stringent minimum prevalence threshold is applied, so any genes with a very small number of counts will be absent.
This is probably more accurate, hence why it's done that way, but is far from assured that they're not actually present at a low level.
Note that this prevalence info is also written to the `--parameter-dir`, so is available even if `--plotdir` isn't set.
