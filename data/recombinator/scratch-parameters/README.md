# SHM parameters for simulation without a sample of your own

One partis SHM parameter dir per locus (`igh/`, `igk/`, `igl/`), for simulating with per-position, per-base mutation when you don't have a parameter dir from your own data:

    partis simulate --rearrange-from-scratch --mutation-multiplier <m> ...

With `--rearrange-from-scratch` and neither `--shm-parameter-dir` nor `--mutate-from-scratch`, partis uses these by default, picking the dir for the locus. Rearrangement then comes from scratch, and mutation from these tables. The tables give:
- each position's rate and the distribution of bases it mutates to (`hmm/mute-freqs/<gene>.csv`);
- how many sequences each gene's numbers came from (`hmm/{v,d,j}_gene-probs.csv`);
- the distribution of mutation levels across families (`hmm/*-mean-{mute-freqs,n-muted}.csv`);
- family-level shm parameters (`hmm/family-shm.yaml`, see below).

`--mutation-multiplier` scales the mutation level. Genes with fewer than `--min-observations-per-gene` observations here get averaged in with similar genes, as for any parameter dir.

## How they're made

`make-scratch-parameters.py` pools the per-gene, per-position mutation counts that `partis cache-parameters` already wrote for public, mutated, human samples. It runs no partis inference itself.

1. For each locus, it takes the parameter dirs listed in a dataset list (`--dataset-fname`), a csv of public datasets with columns:
   - `locus`;
   - `dataset`: the dataset's name, which groups its samples for step 2;
   - `glob`: a pattern matching the parameter dirs of the dataset's samples, relative to `--base-dir`;
   - `source`: where the data is published, copied to the sample lists;
   - `exclude`: if not empty, parameter dirs whose path contains it are skipped.

   The list isn't in the repo, because its globs name where our processed samples are. The script skips samples whose mean mutation frequency is below `--min-mean-mute-freq` (default 0.02), which are mostly naive.
2. It goes round the datasets in turn, adding one randomly chosen sample from each, until the locus has `--max-seqs` sequences (default 1M). A sample with more than `--max-seqs-per-sample` sequences (default 10k) has its counts scaled down to that size, so no small set of very deep samples dominates.
3. It keeps only genes that are in partis's default germline set (`data/germlines/human`) with an identical sequence. That drops each sample's inferred alleles and keeps positions aligned across samples.
   Like any parameter dir inferred from real data, the tables inherit whatever each sample's germline inference missed.
4. For each kept gene and position, it sums each base's counts over samples, then recomputes the mutation frequency and its uncertainty as `mutefreqer` does. To keep the files small, it writes only the columns partis reads from them (position, mutation frequency and its uncertainty, and the count of each base), not `mutefreqer`'s per-base frequencies, which can be recomputed from the counts. It also sums gene counts, and adds the mutation-level histograms with each sample weighted by its number of sequences. It does the same for each sample's mean mutated fraction per region (`mute-freq-sums.csv`), which sets the length simulation uses to turn the number-of-mutations histograms into fractions. These samples' parameter dirs predate that file, so each sample's mean is estimated from its histograms (`utils.get_mean_mfreq()`).
5. It writes `<locus>/hmm/`, and writes the list of samples used, with their sizes, weights, mean mutation frequencies and sources, to `--sample-list-dir` (also outside the repo, since it refers to our filesystem).

To rebuild, from a machine that has the processed samples:

    ./make-scratch-parameters.py --dataset-fname <dataset list> --base-dir <dir containing the datasets in it> --sample-list-dir <dir outside the repo>

### Samples in the current tables

Number of samples, and weighted number of sequences (each sample's sequences times its weight from step 2), from each public source, from the sample lists in `--sample-list-dir` of the build that made the current tables:

| source | locus | per-position tables: samples | weighted seqs | family-level shm parameters: samples | weighted seqs |
|---|---|---|---|---|---|
| Rodriguez et al. 2023 Nat Commun (IgG 5' RACE and FW1 primer) | IGH | 98 | 475,601 | | |
| Tang et al., BioProjects PRJNA591804 and PRJNA381394 (memory, marginal zone and plasma cells) | IGH | 14 | 135,261 | | |
| GEO GSE253862 (bone marrow plasma cells) | IGH | 44 | 95,333 | 43 | 95,276 |
| | IGK | 43 | 82,243 | 42 | 82,193 |
| | IGL | 44 | 14,457 | 42 | 14,425 |
| GEO GSE211869 (germinal center, PBMC, lymph node and bone marrow) | IGH | 3 | 29,998 | 3 | 29,998 |
| | IGK | 3 | 30,000 | 3 | 30,000 |
| | IGL | 3 | 30,002 | 3 | 30,002 |
| OAS paired, E-MTAB-11174 and E-MTAB-11697 (memory and plasmablasts) | IGH | 2 | 20,001 | 1 | 9,999 |
| | IGK | 2 | 20,001 | 1 | 9,998 |
| | IGL | 2 | 20,000 | 1 | 10,001 |
| PairedAbNGS (human BioProjects) | IGH | 98 | 244,897 | 328 | 867,061 |
| | IGK | 605 | 872,936 | 501 | 876,452 |
| | IGL | 947 | 942,499 | 765 | 944,864 |

The family-level parameters use only samples with partitions, so the Rodriguez and Tang samples, which have none, aren't in them.

## Family-level shm parameters

`hmm/family-shm.yaml` (written by `partis/familyshm.py`) holds, for each family size bin (1, 2, 3-5, 6-20, 21+):
- the distribution of each family mean's percentile within its sample's per-sequence mutation distribution, which simulation uses to draw each family's depth according to its size;
- the within-family relative variance of mutation level (variance over the family's sequences divided by the squared family mean, summed over families);

and, over all families, the split-half fraction: the share of within-family variance that does not come from random mutation counts, measured by splitting each sequence into alternating blocks of 30 positions and taking the covariance between the two halves' deviations from the family mean. Simulation multiplies each tree's branches by random factors so that its leaf depths have the within-family relative variance times this fraction (bppseqgen then adds the random mutation counts).

A parameter dir only has this file if its parameters were counted on a partition (`cache-parameters --input-partition-fname`, or `partition --count-parameters`, which writes `multi-hmm/`). Simulation with a parameter dir that lacks it (or lacks families in some size bin) uses this pooled one instead. All three quantities are relative to each sample's own mutation level, so the pooled numbers apply to samples with any mutation level.

The pool uses the samples in the dataset list that have a partition file (`<sample>/partition-<locus>.yaml`), chosen and weighted as in steps 1 and 2 above. It counts each partition file's families with the same code that parameter counting on a partition uses. To remake only `<locus>/hmm/family-shm.yaml`:

    ./make-scratch-parameters.py --dataset-fname <dataset list> --base-dir <dir> --sample-list-dir <dir> --only-family-shm

These tables are from productive repertoires, so they include the effects of selection, not just the mutation process.

They replace an earlier IGH-only dir (from 2017) that had per-gene tables for 39 IGHV genes.
