  - [annotate](#annotate) find most likely annotation for each sequence
  - [partition](#partition) cluster sequences into clonally-related families, and annotate each family
  - [view-output](#view-output) print the partitions and/or annotations from an existing output file
  - [cache-parameters](#cache-parameters) write parameter values and HMM model files for a given data set (if needed, runs automatically before annotation and partitioning)
    - [germline sets](#germline-sets)
  - [simulate](#simulate) make simulated sequences
  - [miscellany](#miscellany)
    - [naive probability estimates](#naive-probability-estimates)

### Subcommands

The `partis` command has a number of actions:

```
partis annotate | partition | view-output | cache-parameters | simulate
```

For information on options for each subcommand that are not documented in this manual run `partis <subcommand> --help`.

For the sake of brevity, the commands below invoke partis with no leading path, which will work if you've either:
  - linked the binary into to a directory that's already in your path, e.g. `ln -s /path/to/<partis_dir/bin/partis ~/bin/` or
  - added the partis `bin/` to your path: `export PATH=/path/to/<partis_dir/bin/:$PATH`
The other option is to specify the full path with each call `/path/to/<partis_dir>/bin/partis`.

### annotate

In order to find the most likely annotation (VDJ assignment, deletion boundaries, etc.) for each single sequence, run

```
partis annotate --infname test/example.fa --outfname _output/example.yaml
```

For information on parsing the output file, see [here](#output-formats.md).

### partition

In order to cluster sequences into clonal families and annotate each family, run

```
partis partition --infname test/example.fa --outfname _output/example.yaml
```

For information on parsing the output file, see [here](#output-formats.md).

By default, this uses the most accurate and slowest method: hierarchical agglomeration with, first, hamming distance between naive sequences for distant clusters, and full likelihood calculation for more similar clusters. This method is typically appropriate (finishing in minutes to tens of hours) for tens of thousands of sequences on a laptop, or hundreds of thousands on a server. Because the usual optimizations that prevent clustering time from scaling quadratically with sample size do not translate well to BCR rearrangement space (at least not if you care about getting accurate clusters), run time varies significantly from sample to sample, hence the large range of quoted run times.

There are also a number of options for speeding things up either by sacrificing accuracy for speed, or by ignoring sequences you don't care about:

##### subsampling

The simplest way to go faster is to only use some of your sequences. You can choose a random subset with `--n-random-queries <n>` (random seed is set with `--seed <n>`), or take the first `n` sequences with `--n-max-queries <n>`. You can also limit it to specific sequence ids using `--queries <a:b:c>`, where `a:b:c` is a colon-separated list of sequence ids.

##### `--seed-unique-id <id>`

Ignores sequences that are not clonally related to the sequence specified by `<id>`, and is consequently vastly faster than partitioning the entire sample. You can specify the sequence rather than the id with `--seed-seq <seq>`.

##### `--naive-vsearch`

As for the default clustering method, first calculates the naive ancestor for each sequence, but then passes these to vsearch for very fast, very heuristic clustering.
Vsearch is extremely fast, because it includes a large number of optimizations to avoid all-against-all comparison, and thus scales much better than quadratically.
Because these optimizations (like any purely distance-based approach) know nothing about VDJ rearrangement or SHM, however, they significantly compromise accuracy.
This is nevertheless a thoroughly reasonably way to get a rough idea of the lineage structure of your sample.
After running vsearch clustering, you can always pass families of interest (e.g. with `--seed-unique-id`) to the more accurate clustering method.

##### `--get-tree-metrics` (and `get-tree-metrics` action)

Construct a phylogenetic tree and use it to calculate tree-based selection metrics (e.g. local branching index) for each cluster in the best partition that's larger than --min-tree-metric-cluster-size (default), and any additional clusters specificied by --write-additional-cluster-annotations, --calculate-alternative-annotations, --min-largest-cluster-size, etc.
Since we want this to be fast enough to run on all the families in a large repertoire, this uses a fairly heuristic approach to calculating trees.
If this is run in the context of partitioning (so there's a clustering path available from hierarchical agglomeration), then that clustering path is used as the starting point for the tree.
It is then refined by replacing any subtrees stemming from large multifurcations with a subtree inferred using FastTree (for instance, the first clustering step is to merge all sequences with similar inferred naive sequences, which results in such subtrees).
If no clustering path information is available, FastTree is used to infer the tree for the entire cluster.
You can also calculate just these trees, without also calculating the tree metrics, by setting `--get-trees`.
If you'd like a more accurate tree, you can infer it separately using your program of choice, and specify it as input to `get-tree-metrics` using `--treefname`.
If you'd like a modern, browser-based package for visualizing the families and their trees and annotations, have a look at our other project, [Olmsted](http://www.olmstedviz.org/).

##### ignore smaller clusters

If you're mostly interested in larger clonal families, you can tell it to cluster as normal for several partitition steps, then discard smaller families (for a description of partition steps, see the paper). Any large families will have accumulated appreciable size within the first few partition steps, and since most real repertoires are dominated by smaller clusters, this will dramatically decrease the remaining sample size. This is turned on by setting `--small-clusters-to-ignore <sizes>`, where `<sizes>` is either a colon-separated list of clusters sizes (e.g. `1:2:3`) or an inclusive range of sizes (e.g. `1-10`). The number of steps after which these small clusters are removed is set with `--n-steps-after-which-to-ignore-small-clusters <n>` (default 3).

##### limit maximum cluster size

Cases where memory is a limiting factor typically stem from a sample with several very large families. Some recent optimizations mean that this doesn't really happen any more, but limiting clonal family size with `--max-cluster-size N` nevertheless can reduce memory usage. Care must be exercised when interpreting the resulting partition, since it will simply stop clustering when any cluster reaches the specified size, rather than stopping at the most likely partition.

##### annotation uncertainties

In order to get an idea of the uncertainty on a given cluster's naive sequence and gene calls, you can specify `--calculate-alternative-annotations` during the partition step.
This will write all the annotations for intermediate sub-clusters to the output file so that it counts up how many subclusters "voted" for each alternative naive sequence or gene call.
Since most annotation uncertainty in large-ish families boils down to two sub-families disagreeing about, say, which is the correct D gene, this approach typically does a decent job of spanning the real uncertainty (despite being quite heuristic).
The resulting information can be accessed either by pulling out the resulting ['alternative-annotations' key](output-formats.md#description-of-keys) from the output file, or by running view-alternative-annotations (or running the partition step with `--debug 1`).
For instance:

```
partis partition --infname test/example.fa --outfname _output/example.yaml --calculate-alternative-annotations
partis view-alternative-annotations --outfname _output/example.yaml  # pipe this to less by adding "| less -RS"
```

If you only care about one cluster, you can print only the cluster corresponding to a given set of queries by setting `--queries <queries of interest>` in the second step.
This second command prints an ascii representation of the various naive sequences and gene calls, as well as summaries of how many unique sequences and clusters of various sizes voted for each naive sequence and gene call.

##### progress file

Since individual steps of the partition algorithm can be somewhat time consuming, during each step we periodically update a file with clustering status information so you can monitor progress.
This file is written to the working directory, which you've either set, or is under `/tmp/$USER/hmms`.
The easiest way to find the working directory is to run an appropriate `ps` command to get the complete command lines of each subprocess (for instance `ps -Ho pid,args` or `ps -eo pcpu,pmem,rss,cputime:12,stime:7,user,pid,args:100 --sort pmem,pcpu`, depending if you're running from the terminal that called partis or not).
For instance running the former in a terminal in which I've backgrounded the partis command gives (so the parent python process is currently clustering with one bcrham subprocess):

```
  PID COMMAND
 2677 -bash
 4685   python ./bin/partis partition --infname test/reference-results/test/simu.yaml
 4734     /home/dralph/work/partis/packages/ham/bcrham --algorithm forward --hmmdir /home/dralph/work/partis/_output/test_reference-results_test_simu/hmm/hmms --datadir /tmp/dralph/hmms/902992/germline-sets --infile /tmp/dralph/hmms/902992/hmm_input.csv --outfile /tmp/dralph/hmms/902992/hmm_output.csv --locus igh --random-seed 1554239923 --input-cachefname /tmp/dralph/hmms/902992/hmm_cached_info.csv --output-cachefname /tmp/dralph/hmms/902992/hmm_cached_info.csv --partition --max-logprob-drop 5.0 --hamming-fraction-bound-lo 0.015 --hamming-fraction-bound-hi 0.0711739130435 --logprob-ratio-threshold 18.0 --biggest-naive-seq-cluster-to-calculate 15 --biggest-logprob-cluster-to-calculate 15 --n-partitions-to-write 10 --write-logprob-for-each-partition --ambig-base N
 6836   ps -Ho pid,args
```

There will be a progress file for each subprocess in the clustering step, in their own subdirectories, and the name of each progress file will look like `<stuff>.progress`.
A few lines from a typical progress file might look like:

```
      Apr 02 14:19:35    1850 clusters        45160 / 32857084 kB =  0.137 %           calcd:   vtb 281   fwd 0     hfrac 154382          merged:  hfrac 1    lratio 0                  clusters:  38 30 15 7 7 7 7 7 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2  (1569 singletons)
      Apr 02 14:29:22    1841 clusters        47068 / 32857084 kB =  0.143 %           calcd:   vtb 281   fwd 0     hfrac 155609          merged:  hfrac 10   lratio 0                  clusters:  38 30 15 7 7 7 7 7 7 6 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2  (1564 singletons)
      Apr 02 14:29:53    1813 clusters        50876 / 32857084 kB =  0.155 %           calcd:   vtb 289   fwd 0     hfrac 159161          merged:  hfrac 38   lratio 0                  clusters:  38 30 15 13 11 9 9 7 7 7 7 7 6 6 6 6 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 5 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2  (1542 singletons)
```

The columns are:
  * time the line was written. There are a number of points in the clustering code that trigger a progress file write, but typically a new line is written once each time it loops through all the current clusters looking for merges.
  * total number of clusters
  * memory allocated to this clustering process (i.e. bcrham).
  * total number of viterbi, forward, and naive hamming fraction calculations that have so far been made
  * total number of cluster merges that have so far been performed, broken down into hamming fraction merges and likelihood ratio (i.e. forward calculation) merges
  * a list of all cluster sizes in the current partition

### view-output

To print the partitions and/or annotations in an existing output file, run

``` 
partis view-output --outfname _output/example.yaml
```

### cache-parameters

This is run automatically if `--parameter-dir` doesn't exist (whether this directory is specified explicitly, or is left as default).
So you do not, generally, need to run it on its own.

When presented with a new data set, the first thing we need to do is infer a set of parameters, a task for which we need a preliminary annotation.
As such, partis first runs ig-sw's Smith-Waterman algorithm on the data set.
The smith-waterman annotations are used to build and write out a parameter set, which is in turn used to make a set of HMM model files for each observed allele.
These files are then passed as input to a second, HMM-based, annotation step, which again outputs (more accurate) parameter values and HMM model files.

To explicitly run this parameter caching step by itself:

``` 
partis cache-parameters --infname test/example.fa --parameter-dir _output/example
```

The resulting parameter csvs from Smith-Waterman and the HMM are put into `/sw` and `/hmm` subdirectories of `--parameter-dir`.
Within each of these, there are a bunch of csv files with (hopefully) self-explanatory names, e.g. `j_gene-j_5p_del-probs.csv` has counts for J 5' deletions subset by J gene.
The hmm model files go in the `hmms/` subdirectory, which contains yaml HMM model files for each observed allele.

If you don't specify `--parameter-dir`, it defaults to a location in the current directory that amounts to a slight bastardization of your input file path (e.g. parameters for `path/to/seqs.fa` will go in `_output/path_to_seqs/`).
This default is designed such that with typical workflows, if your input files have different paths, their parameters will go in different places.
If you run once, and then put different sequences in the same input file, this convention obviously horribly breaks (you'll be applying the parameters from the sequences in the old file contents to the new file).

If `--parameter-dir` (whether explicitly set or left as default) doesn't exist, partis assumes that it needs to cache parameters, and does that before running the requested action.

Whether caching parameters or running on pre-existing parameters, the hmm needs smith-waterman annotations as input.
While this preliminary smith-waterman step is fairly fast, it's also easy to cache the results so you only have to do it once.
By default these smith-waterman annotations are written to a yaml file in `--parameter-dir` during parameter caching (default path is `<parameter_dir>/sw-cache-<hash>.yaml` where `<hash>` is a hash of the input sequence ids).
(Because all sequences need to be aligned and padded to the same length before partititioning, the smith-waterman annotation information for each sequence depends slightly on all the other sequences in the file, hence the hash.)
These defaults should ensure that with typical workflows, smith-waterman only runs once.
If however, you're doing less typical things (running on a subset of sequences in the file), if you want smith-waterman results to be cached you'll need to specify `--sw-cachefname` explicitly, and it'll write it if it doesn't exist, and read from it if it does.

#### germline sets

By default partis infers a germline set for each sample during parameter caching, using as a starting point the germline sets in `data/germlines/<--species>`.
The resulting per-sample germline sets are written both to the output yaml file (if you've set `--outfname`), and to `<--parameter-dir>/hmm/germline-sets` (as three fasta files and a meta-info csv).
To start from a non-default germline set, use `--initial-germline-dir <dir>.`
This dir must use the same format as `data/germlines/<--species>`, i.e. have a subdir for the `--locus` you want to run, and in that subdir have a fasta file for v, d, and j, and a csv file with the (zero-indexed!) conserved codon positions for each v/j gene.
The easiest way to construct this is to copy the files from `data/germlines/<--species>`, and then substitute the fasta files for your own.
Partis aligns any genes that are missing codon information against genes with known info, so as long as there's some overlap with the default germline, it should be able to figure out a position for each gene.

If you want to be more careful, it's a good idea to check how it's figuring out the extras.csv info.
You can do this (after modifying the new dir) by running the following lines, for example by adding them to the top of bin/example-output-processing.py (just after it imports glutils):

```
ref_glfo = glutils.read_glfo('data/germlines/human', locus='igh', debug=True)
new_glfo = glutils.read_glfo('new-germline-dir', locus='igh', template_glfo=ref_glfo, debug=2)  # adds missing codon info
print 'cleaned up glfo:'
glutils.print_glfo(new_glfo)  # print the cleaned up glfo
glutils.write_glfo('cleaned-up-new-germline-dir', new_glfo, debug=True)
sys.exit()
```

Where we pass the `ref_glfo` as a template glfo, which makes it more likely to be able to figure out the cyst/tryp/phen info. And you'll probably want to pipe to `|less -RS`, especially if on a laptop.


By default, if `--species` is set to human this only looks for alleles that are separated by point mutations from existing genes.
This is appropriate for humans, since the known germline set is fairly complete.
For species for which the known germline sets are much less complete (e.g. macaque and mice), it makes sense to also look for alleles that are separated by indels and/or many SNPs from existing genes.
This is the default for non-human species, and is equivalent to turning on `--allele-cluster`.

At the moment we only infer germline sets for V, and not for D and J (although we plan to fix this).
This isn't horribly nonsensical, since there is much less variation in D and J, but it does mean that you should treat the D germline set, in particular, with quite a bit of skepticism.

We unfortunately do not yet have a way to report simple confidence estimates for each gene in the inferred germline set.
There is, however, a wealth of information that can be used to get a good sense of this confidence [here](germline-inference.md).

### simulate

Note that in order to run simulation, you have to have R installed, along with several [extra packages](install.md#simulation).

In the simplest simulation mode, partis mimics the characteristics of a particular template data sample as closely as possible: germline set, gene usage frequencies, insertion and deletion lengths, somatic hypermutation rates and per-position dependencies, etc. (as well as the correlations between these).
This mode uses the previously-inferred parameters from that sample, located in `--parameter-dir`.
By default, for instance, if a sample at the path `/path/to/sample.fa` was previously partitioned, the parameters would have been written to `_output/_path_to_sample/`.
You could thus write the mature sequences resulting from three simulated rearrangement events to the file `simu.yaml` by running

```partis simulate --parameter-dir _output/_path_to_sample --outfname simu.yaml --n-sim-events 3 --debug 1```,

where `--debug 1` prints to stdout what the rearrangement events look like as they're being made.
The resulting output file follows regular output [format](output-formats.md), with an additional column `reco_id` to identify clonal families (it's a hash of the rearrangement parameters).
When subsequently running inference on simulation, you typically want to pass the `--is-simu` option.
During parameter caching, this will write a separate parameter directory with the true parameters (in a addition to `sw/` and `hmm/`).
During annotation and partitioning, with `--debug 1` it will print the true rearrangements and partitions along with the inferred ones.

There are a wide variety of options for manipulating how the characteristics of the simulation deviate from the template data sample (for information on defaults, run `partis simulate --help`).

**Miscellaneous:**

| option                           | description
|----------------------------------|-----------------------------------------------------------------
| `--mutation-multiplier <factor>` | multiply the observed SHM rate by `<factor>`
| `--mimic-data-read-length`       | by default the simulation creates reads that extend through all of V and J. This tells it, instead, to truncate on the 5' and 3' ends according to the lengths/frequencies seen in the template data sample.

**Tree control:**

| option                                        | description
|-----------------------------------------------|-----------------------------------------------------------------
| `--n-trees <N>`                               | Before actually generating events, we first make a set of `<N>` phylogentic trees. For each event, we then choose a tree at random from this set. Defaults to the value of --n-sim-events.
| `--n-leaf-distribution <geometric,box,zipf>`  | When generating these trees, from what distribution should the number of leaves be drawn?
| `--n-leaves <N>`                              | Parameter controlling the n-leaf distribution (e.g. for the geometric distribution, it's the mean number of leaves)
| `--constant-number-of-leaves`                 | instead of drawing the number of leaves for each tree from a distribution, force every tree to have the same number of leaves
| `--root-mrca-weibull-parameter`               | adjusts tree balance/speciation by switching to TreeSimGM (useful range: 0.3 to 1.3)
| `--n-trees <N>`                               | Before actually generating events, we first make a set of `<N>` phylogentic trees. For each event, we then choose a tree at random from this set. Defaults to the value of --n-sim-events.
| `--input-simulation-treefname <N>`            | File with list of trees as newick-formatted lines, from which to draw for each event (instead of generating with TreeSim or TreeSimGM). Note that a) the tree depths are rescaled according to the shm rates requested by other command line arguments, i.e. the depths in the tree file are ignored, and b) the resulting sequences do not use the leaf names from the trees.

**SHM indel control:**

| option                          | description
|---------------------------------|-----------------------------------------------------------------
| `--indel-frequency <f>`         | fraction of simulated sequences which will contain SHM indels (currently, insertions and deletions are generated with equal probability, although this would be easy to change)
| `--n-indels-per-indeld-seq <l>` | once we've decided a sequence will have at least one indel, we choose the actual number of indels from this colon-separated list of integers
| `--mean-indel-length <N>`       | mean length of each SHM insertion or deletion
| `--indel-location <v,cdr3,<pos>>`| if not set (default), indels are placed uniformly over the whole sequence. If set to `v` or `cdr3` indels are restricted to that portion of the sequence. If set to a `<pos>`, all indels occur at that position in the sequence.

**Scratch parameters:**

In order to deviate more profoundly from the template data sample (or to use no template sample at all), there are a number of options centered around `--rearrange-from-scratch`:

| option                            | description
|-----------------------------------|-----------------------------------------------------------------
| `--rearrange-from-scratch`        | instead of taking rearrangement-level (i.e. non-SHM) parameters from --parameter-dir, make up some reasonable values from scratch (e.g. geometric insertion lengths, etc.)
| `--scratch-mute-freq-dir <path>`  | parameter directory with only SHM-level information, which allows --rearrange-from-scratch to create realistic mutation distributions for any specified germline set (the default, in data/recombinator/scratch-parameters/ shouldn't need to be changed)
| `--mutate-from-scratch`           | instead of taking realistic SHM-level (i.e. non-rearrangement level) parameters from --parameter-dir or --scratch-mute-freq-dir, use a simple flat rate across positions and genes (value is specified by --flat-mute-freq)
| `--flat-mute-freq`                | see --mutate-from-scratch

**Germline set control:**

By default, the germline set (set of germline V, D, and J genes) used for simulation, and their prevalence frequencies, are taken from the template data sample (i.e. from `<--parameter-dir>/hmm/germline-sets/<locus>`.
However, if you have a particular germline set that you want to use, that can be specified with `--initial-germline-dir` (the format, of three fastas and a csv, should mimic that in data/germlines/human/igh).
You can restrict the genes that are then actually used for simulation with `--only-genes`:

| option                        | description
|-------------------------------|-----------------------------------------------------------------
| `--initial-germline-dir <dir>`| simulate with the germline set in `<dir>`, instead of the one from --parameter-dir
| `--only-genes <gene-list>`    | restrict the germline set to `<gene-list>`, specified as a colon-separated list of genes, for instance `IGHV3-53*03:IGHJ3*02` (any regions that have no genes in the list, like the D region in this example, will be unrestricted).

Instead of modifying an existing per-sample germline set with the options above, you can also direct partis to generate the germline set by choosing genes from an existing species-wide set with `--generate-germline-set` (use --help to see default values).
The species-wide germline set defaults to the imgt set, but can be set with `--initial-germline-dir`.

| option                                 | description
|----------------------------------------|-----------------------------------------------------------------
| `--generate-germline-set`              | generate a realistic germline set from scratch, rather than mimicking an existing germline set (`--rearrange-from-scratch` must also be set)
| `--n-genes-per-region <m:n:q>`         | number of genes to choose for each of the V, D, and J regions (colon-separated list ordered like v:d:j)
| `--n-sim-alleles-per-gene <stuff>`     | mean number of alleles to choose for each gene, for each region (colon-separated list, e.g. '1.3:1.2:1.1' will choose either 1 or 2 alleles for each gene with the proper probabilities such that the mean alleles per gene is 1.3 for V, 1.2 for D, and 1.1 for J)
| `--min-sim-allele-prevalence-freq <f>` | minimum prevalence ratio between any two alleles in the germline set. I.e., the prevalence for each allele is chosen such that the ratio of any two is between `<f>` and 1
| `--allele-prevalence-fname`            | not really designed to be modified or used by hand, but used by `--allele-prevalence-freqs` option to `bin/test-germline-inference.py` (see below)


Details of the generated germline set will be printed to stdout, and after simulation the prevalence frequencies are also checked and printed.

**Generating novel alleles:**

There are also several ways to generate novel alleles, i.e. alleles that are not in an existing germline set.
Because this is somewhat complicated, we describe these options using the helper script `bin/test-germline-inference.py`, which automates a simulation run and subsequent partis germline inference on that simulation (run `./bin/test-germline-inference.py --help` for examples).

You first need to either give it an explicit list of genes to use, or tell it to generate a germline set from scratch:

| option                   | description
|--------------------------|-----------------------------------------------------------------
| `--sim-v-genes <v-list>` | colon-separated list of V genes to use for simulation
| `--inf-v-genes <v-list>` | start from this list of V genes, and try to infer the genes from --sim-v-genes
| `--dj-genes <dj-list>`   | D and J genes used for both simulation and inference
| `--gls-gen`              | instead of using explicit gene lists, generate a full germline set from scratch (see --generate-germline-set in `partis simulate --help` for details)

You can then add novel alleles to the germline set by telling it how many novel alleles, with how many SNPs and/or indels, and where to introduce the SNPs/indels:

| option                                  | description
|-----------------------------------------|-----------------------------------------------------------------
| `--nsnp-list <m:n:...>`                 | list of the number of SNPs to generate for each novel allele (each at a random position in the sequence). If --gls-gen is not set, length must equal length of <--sim-v-genes>.  E.g. '0:1:3' will generate two novel alleles, separated by 1 and 3 SNPs from the second and third genes in --sim-v-genes
| `--nindel-list <m:n:...>`               | same as --nsnp-list, but for indels
| `--snp-positions <stuff>`               | colon-separated list of comma-separated SNP positions for each gene, e.g. '3,71:45' will generate two novel alleles, separated by two SNPs (at zero-indexed sequence positions 3 and 71) and one SNP (at 45) from the two genes in --sim-v-genes.
| `--indel-positions <m:n:...>`           | same as --snp-positions, but for indels
| `--allele-prevalence-freqs <f1:f2:...>` | colon-separated list of allele prevalence frequencies, including newly-generated snpd genes (ordered alphabetically)


### Miscellany

##### naive probability estimates

The script `bin/get-naive-probabilities.py` is designed to answer the question "What is the probability of a particular rearrangement in this sample?".
It uses the `all-probs.csv` file that is written to each parameter directory, which contains a line counting how many times we saw each unique rearrangement in that sample.
Using a yaml config file, you can have the script use a variety of criteria, for instance, find the probability of rearrangements that use any IGHV1 family gene, and that have CDR3 length of 60, 63, or 66.
For more details run `./bin/get-naive-probabilities.py --help`.

