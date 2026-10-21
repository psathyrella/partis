# SHM parameters for simulation without a sample of your own

One partis SHM parameter dir per locus (`igh/`, `igk/`, `igl/`), for simulating with per-position, per-base mutation when you don't have a parameter dir from your own data:

    partis simulate --rearrange-from-scratch --mutation-multiplier <m> ...

With `--rearrange-from-scratch` and neither `--shm-parameter-dir` nor `--mutate-from-scratch`, partis uses these by default, picking the dir for the locus. Rearrangement then comes from scratch, and mutation from these tables. The tables give:
- each position's rate and the distribution of bases it mutates to (`hmm/mute-freqs/<gene>.csv`);
- how many sequences each gene's numbers came from (`hmm/{v,d,j}_gene-probs.csv`);
- the distribution of mutation levels across families (`hmm/*-mean-{mute-freqs,n-muted}.csv`).

`--mutation-multiplier` scales the mutation level. Genes with fewer than `--min-observations-per-gene` observations here get averaged in with similar genes, as for any parameter dir.

## How they're made

`make-scratch-parameters.py` pools the per-gene, per-position mutation counts that `partis cache-parameters` already wrote for public, mutated, human samples. It runs no partis inference itself.

1. For each locus, it takes the parameter dirs listed in `datasets.csv` (public datasets only, with their sources). It skips samples whose mean mutation frequency is below `--min-mean-mute-freq` (default 0.02), which are mostly naive.
2. It goes round the datasets in turn, adding one randomly chosen sample from each, until the locus has `--max-seqs` sequences (default 1M). A sample with more than `--max-seqs-per-sample` sequences (default 50k) has its counts scaled down to that size, so many individuals contribute rather than a few very deep ones.
3. It keeps only genes that are in partis's default germline set (`data/germlines/human`) with an identical sequence. That drops each sample's inferred alleles and keeps positions aligned across samples.
4. For each kept gene and position, it sums each base's counts over samples, then recomputes the mutation frequency and its uncertainty as `mutefreqer` does. To keep the files small, it writes only the columns partis reads from them (position, mutation frequency and its uncertainty, and the count of each base), not `mutefreqer`'s per-base frequencies, which can be recomputed from the counts. It also sums gene counts, and adds the mutation-level histograms with each sample weighted by its number of sequences.
5. It writes `<locus>/hmm/` and `<locus>/samples.csv`, the list of samples used with their sizes, mean mutation frequencies and sources.

To rebuild, from a machine that has the processed samples:

    ./make-scratch-parameters.py --base-dir <dir containing the datasets in datasets.csv>

These tables are from productive repertoires, so they include the effects of selection, not just the mutation process.

They replace an earlier IGH-only dir (from 2017) that had per-gene tables for 39 IGHV genes.
