[Up to table of contents](contents.md)

### output file formatting

  * [output file overview](#output-file-overview)
  * [paired loci](paired-loci.md#output-directory)
  * [extracting simplified fasta or csv/tsv files](#extracting-simplified-files)
  * [N-padding, read length, etc](#N-padding)
  * [description of keys](#description-of-keys)
  * [output file example](#output-file-example)

All output is written to a unified yaml file (for documentation on the old csv formats, see [here](https://github.com/psathyrella/partis/blob/9d78600ba6b7d511825f27724de1f8e937a5ada3/docs/output-formats.md)).
The `partition` action writes a list of the most likely partitions and their relative likelihoods, as well as annotations for each cluster in the most likely partition.
You can write additional less-likely partitions with `--n-partitions-to-write`, as well as annotations for clusters in less-likely partitions with `--write-additional-cluster-annotations`.
Note that you should always access clusters using first the partition list and then looking for the cluster's annotation in the annotation list, and not by first looking in the annotation list.
In some cases the annotation list may correspond to the most likely partition, but there are many cases where it does not (e.g. if `--calculate-alternative-annotations` or `--write-additional-cluster-annotations` are set).
The `annotate` action, on the other hand, only writes single-sequence annotations for each sequence in the input.

If you want to print the results of existing output files to the terminal, use the partis [`view-output`](subcommands.md#view-output) action.
While by default a fairly minimal set of annotation information is written to file, many more keys are present in the dictionary in memory (see below).
Any of these keys, together with several additional ones, can be added to the output file by setting `--extra-annotation-columns key_a:key_b` (for all the choices, see below, or run `partis annotate --help|grep -C5 extra-annotation`).

An example parsing script can be found [here](../bin/parse-output.py).

To have partis write to an AIRR-format tsv file, set the partis option `--airr-output`.
To convert existing partis output to AIRR tsv, pass the same option to `bin/parse-output.py`.
Both `bin/partis` and `bin/parse-output.py` also take AIRR files as input with the flag `--airr-input`.

For more information on all options, run `partis <action> --help`.

#### output file overview

The yaml output file contains four top-level headers: 

|   name         |  description
|----------------|----------------------------------------------------------------------------
|  version-info  |  output file format version
|  germline-info |  germline sequence, names, and conserved codon positions
|  events        |  list of annotations for each rearrangement event (i.e. group of clonally-related sequences). Can have multiple annotations for an event, and can have annotations for events not in the best partition, i.e. do not "guess" the partition based on the events here.
|  partitions    |  list of partitions and their associated log probabilities. The most likely partition has the highest logprob.

#### extracting simplified files

In order to quickly extract sequences (plus limited other info) from partis output files to fasta or csv/tsv, you can use `bin/parse-output.py`.
For example
```
parse-output.py test/reference-results/partition-new-simu.yaml <tmp.csv|tmp.fa> --extra-columns cdr3_length:naive_seq
```
will write input sequences, together with inferred naive sequences and cdr3 lengths, to `tmp.csv` or `tmp.fa`.
See `parse-output.py --help` for details.

The ClusterPath class, which represents a series of partitions, is useful for handling partitions.
This snippet reads a cluster path from a file, prints an ascii summary, and gets the best partition:

```
from clusterpath import ClusterPath
_, _, cpath = utils.read_output('test/ref-results/partition-new-simu.yaml')
cpath.print_partitions()
best_ptn = cpath.best()  # return best partition
```

#### N-padding, read length, and "non-biological" insertions and deletions

As detailed elsewhere, partis first annotates all sequences with a Smith-Waterman (SW) algorithm, then uses these annotations to build and run an HMM to get final HMM annotations.
During SW annotation, anything in the observed sequences that is 5' of V or 3' of J is trimmed off, with the removed bits saved in `leader_seqs` and `c_gene_seqs` (see table below).
Observed sequences that, on the other hand, do not extend to the 5' end of V or 3' end of J are recorded as having "non-biological" deletions (`v_5p_del` and `j_3p_del`).
Before beginning the HMM, the SW annotations are N-padded such that all sequences with the same CDR3 length have the same number of bases both to 5' and 3' of the CDR3.
This is necessary since the HMM can only operate on same-length sequences (and only considers sequences as potentially clonal that have the same CDR3 length).
These N pads are recorded as "non-biological" insertions (`fv_insertion` and `jf_insertion`), and they mean that the final HMM annotations won't have V 5' or J 3' deletions (they're instead filled with Ns).

#### description of keys

Keys in the annotation dictionary are either per-family keys (that have one value for the entire rearrangement event) or per-sequence keys (that consist of a list of values, one for each sequence).
The latter are marked with `[per-seq]` below.

The following keys are written to output by default:

|   name         |  description
|----------------|----------------------------------------------------------------------------
| unique_ids     |  list of sequence identification strings `[per-seq]`
| reco_id        |  simulation only: hash of rearrangement parameters that is the same for all clonally-related sequences
| v_gene         |  V gene in most likely annotation
| d_gene         |  see v_gene
| j_gene         |  see v_gene
| cdr3_length    |  nucleotide CDR3 length of most likely annotation, but note that this __includes__ both conserved codons in their entirety, i.e. is what IMGT calls the ["junction length"](https://docs.airr-community.org/en/stable/datarep/rearrangements.html#junction-versus-cdr3)
| mut_freqs      |  list of sequence mutation frequencies `[per-seq]`
| input_seqs     |  list of input sequences (with constant regions (fv/jf insertions) removed, unless `--dont-remove-framework-insertions` was set) `[per-seq]`
| naive_seq      |  naive (unmutated ancestor) sequence corresponding to most likely annotation
| v_3p_del       |  length of V 3' deletion
| d_5p_del       |  length of D 5' deletion
| d_3p_del       |  length of D 3' deletion
| j_5p_del       |  length of J 5' deletion
| v_5p_del       |  length of a non-biological "effective" V 5' deletion, corresponding to a read that doesn't extend to the 5' end of V (only present in smith-waterman annotations; it's replaced with Ns in the hmm annotations)
| j_3p_del       |  non-biological "effective" deletion on 3' end of J (see v_5p_del)
| vd_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the V and D segments
| dj_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the D and J segments
| leader_seqs        |  sequence to 5' of V `[per-seq]` (this is trimmed off of `seqs` during smith-waterman alignment)
| c_gene_seqs        |  sequence to 3' of J `[per-seq]` (this is trimmed off of `seqs` during smith-waterman alignment)
| leaders            |  leader gene that was the best alignment to each seq in leader_seqs (only filled if --align-constant-regions is set) `[per-seq]`
| c_genes            |  constant gene that was the best alignment to each seq in c_gene_seqs (only filled if --align-constant-regions is set) `[per-seq]`
| fv_insertion       |  N-padded sequence added to 5' of V (if necessary, see above)
| jf_insertion       |  N-padded sequence added to 3' of J (if necessary, see above)
| codon_positions	 |  zero-indexed indel-reversed-sequence positions of the conserved cyst and tryp/phen codons that define the start/end of the CDR3 region, e.g. `{'v': 285, 'j': 336}`
| mutated_invariants |  true if either of the conserved codons corresponding to the start and end of the CDR3 code for a different amino acid than their original germline (cyst and tryp/phen, in IMGT numbering) `[per-seq]`
| in_frames          |  true if the net effect of VDJ rearrangement and SHM indels leaves both the start and end of the CDR3 (IMGT cyst and tryp/phen) in frame with respect to the start of the germline V sequence `[per-seq]`
| stops              |  true if there's a stop codon in frame with respect to the start of the germline V sequence `[per-seq]`
| v_per_gene_support |  approximate probability supporting the top V gene matches, as a list of lists (or ordered dict) of gene:probability pairs. Only includes the V genes that the Smith-Waterman step decided to pass to the hmm (which is quite heuristic/non-probabilistic) so should only be used to compare genes that appear in it (i.e. the absence of a gene means it wasn't passed to the hmm, not necessarily that its probability was zero). Entirely separate from 'alternative-annotations' below, which is probably more accurate.
| d_per_gene_support |  see v_per_gene_support
| j_per_gene_support |  see v_per_gene_support
| indel_reversed_seqs  |  list of input sequences with indels reversed/undone, and with constant regions (fv/jf insertions) removed. Empty string if there are no indels, i.e. if it's the same as 'input_seqs' `[per-seq]`
| gl_gap_seqs        |  list of germline sequences with gaps at shm indel positions (alignment matches qr_gap_seqs) `[per-seq]`
| qr_gap_seqs        |  list of query sequences with gaps at shm indel positions (alignment matches gl_gap_seqs) `[per-seq]`
| duplicates     |  list of "duplicate" sequences for each sequence. If --collapse-duplicate-sequences is set, then after trimming fv/jf insertions, any identical sequences are collapsed during the smith-waterman step (see also the input meta info [multiplicity key](subcommands.md#input-meta-info), as well as --also-remove-duplicate-sequences-with-different-lengths and --dont-remove-framework-insertions). `[per-seq]`
| tree           |  simulation only: newick-formatted string of the true phylogenetic tree (inferred trees are included in `tree-info`, in order to accomodate multiple trees inferred by different methods)
| tree-info      |  inferred tree-related information from various methods (e.g. local branching index/ratio, cons-dist-aa), including associated inferred trees. Written for 'get-selection-metrics' action or when --get-selection-metrics or --get-trees are set. Also can include distance to consensus sequence, since this is used similarly to the actual tree metrics.
| alternative-annotations | summary of alternative annotation information (a.t.m. naive sequences and gene calls), where counts of unique sequences have been normalized to give a number between 0 and 1 that can be interpreted as a (heuristically-derived!) probability of each potential naive sequence and gene call. See --calculate-alternative-annotations and the [view-alternative-annotations](subcommands.md#annotation-uncertainties-alternative-annotations) action for details.

The following keys are available in the dictionary in memory, but not written to disk by default (can be written by setting `--extra-annotation-columns key_a:key_b`):

|   key                   |  value
|-------------------------|----------------------------------------------------------------------------
| v_gl_seq				  |  portion of v germline gene aligned to the indel-reversed sequence (i.e. with 5p and 3p deletions removed)
| d_gl_seq				  |  see v_gl_seq
| j_gl_seq				  |  see v_gl_seq
| v_qr_seqs               |  portion of indel-reversed sequence aligned to the v region `[per-seq]`
| d_qr_seqs				  |  see v_qr_seqs `[per-seq]`
| j_qr_seqs				  |  see v_qr_seqs `[per-seq]`
| lengths				  |  lengths aligned to each of the v, d, and j regions, e.g. `{'j': 48, 'd': 26, 'v': 296}` (equal to lengths of `[vdj]_qr_seqs` and `[vdj]_gl_seq`)
| regional_bounds		  |  indices in the observed sequences corresponding to the boundaries of the v, d, and j regions (python slice conventions), e.g. `{'j': (322, 370), 'd': (296, 322), 'v': (0, 296)}`
| aligned_v_seqs		  |  list of indel-reversed sequences aligned to germline sequences given by `--aligned-germline-fname`. Only used for presto output `[per-seq]`
| aligned_d_seqs          |  see aligned_v_seqs `[per-seq]`
| aligned_j_seqs		  |  see aligned_v_seqs `[per-seq]`
| invalid				  |  indicates an invalid rearrangement event

The following keys can also be added to the output file using `--extra-annotation-columns key_a:key_b`:

|   name                  |  description
|-------------------------|----------------------------------------------------------------------------
| cdr3_seqs				  |  nucleotide CDR3 sequence, including bounding conserved codons `[per-seq]`
| full_coding_naive_seq	  |  in cases where the input reads do not extend through the entire V and J regions, the input_seqs and naive_seq keys will also not cover the whole coding regions. In such cases full_coding_naive_seq and full_coding_input_seqs can be used to tack on the missing bits.
| full_coding_input_seqs  |  see full_coding_naive_seq `[per-seq]`
| cons_dists_nuc          |  nucleotide distance to clonal family consensus sequence
| cons_dists_aa           |  amino acid distance to clonal family consensus sequence
| consensus_seq           | nucleotide consensus sequence for the family (calculated from "indel_reversed_seqs", *not* from "input_seqs", i.e. any shm indels are reversed; i.e. assumes that indels are in a minority of the family, which could be incorrect)
| consensus_seq_aa        | amino acid consensus sequence for the family (calculated from "indel_reversed_seqs", *not* from "input_seqs", i.e. any shm indels are reversed; i.e. assumes that indels are in a minority of the family, which could be incorrect) 
| seqs_aa                 |  amino acid translations of the nucleotide sequence under the 'indel_reversed_seqs' key
| naive_seq_aa            |  amino acid translation of 'naive_seq'

Partitioning results in a list of partitions, with one line for the most likely partition (the one with the highest logprob), as well as a number of lines for the surrounding less-likely partitions.
The number of partitions surrounding the best partition that are written can be configured with `--n-partitions-to-write N` (many other aspects of partitioning can also be configured, see `partis partition --help`).
It also writes the annotation for each cluster in the most likely partition (you can tell it to also write the annotations for clusters in other, less-likely, partitions by setting `--write-additional-cluster-annotations m:n`, where m (n) are integers specifying the number of partitions before (after) the best partition.
Reading these partitions is best accomplished using the ClusterPath class, as in the example parsing script [here](../bin/example-parse-output.py).
The following keys describe the partitions:

|   column header  |  description
|------------------|----------------------------------------------------------------------------
| logprob          |  Total log probability of this partition
| n_clusters       |  Number of clusters (clonal families in this partition)
| partition        |  String representing the clusters, where clusters are separated by ; and sequences within clusters by :, e.g. 'a:b;c:d:e'
| n_procs          |  Number of processes which were simultaneously running for this clusterpath. In practice, final output is usually only written for n_procs = 1

#### output file example

The following file contains the partitions and annotations for three sequences with ids 'a', 'b', and 'c'.
There are two partitions, one with 'a' by itself and 'b' and 'c' together; then the most likely partition where all three are together.
The annotations are for the most likely partition, and thus describe a single rearrangement event with three sequences.
Note that while this examples is in full yaml (since it's more human readble), by default we read and write output files using the json subset of yaml because it's much faster.
To instead write full yaml output files, set `--write-full-yaml-output`.

```
version-info: {partis-yaml: 0.1}
germline-info:
  cyst-positions: {IGHV3-48*04: 285, IGHV3-74*01: 285, IGHV4-31*10: 288}
  functionalities: {}
  locus: igh
  seqs:
    d: !!python/object/apply:collections.OrderedDict
    - - [IGHD1-20*01, GGTATAACTGGAACGAC]
      - [IGHD2-2*01, AGGATATTGTAGTAGTACCAGCTGCTATGCC]
      - [IGHD5-18*01, GTGGATACAGCTATGGTTAC]
    j: !!python/object/apply:collections.OrderedDict
    - - [IGHJ3*02, TGATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG]
      - [IGHJ4*01, ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG]
      - [IGHJ6*01, ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG]
    v: !!python/object/apply:collections.OrderedDict
    - - [IGHV3-48*04, GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTACAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTTCATACATTAGTAGTAGTAGTAGTACCATATACTACGCAGACTCTGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGA]
      - [IGHV3-74*01, GAGGTGCAGCTGGTGGAGTCCGGGGGAGGCTTAGTTCAGCCTGGGGGGTCCCTGAGACTCTCCTGTGCAGCCTCTGGATTCACCTTCAGTAGCTACTGGATGCACTGGGTCCGCCAAGCTCCAGGGAAGGGGCTGGTGTGGGTCTCACGTATTAATAGTGATGGGAGTAGCACAAGCTACGCGGACTCCGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACACGCTGTATCTGCAAATGAACAGTCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCAAGAGA]
      - [IGHV4-31*10, CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTACCATATCAGTAGACCCGTCCAAGAACCAGTTCTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGAGA]
  tryp-positions: {IGHJ3*02: 16, IGHJ4*01: 14, IGHJ6*01: 29}
partitions:
- logprob: -174.72939717720863
  n_clusters: 2
  n_procs: 1
  partition:
  - [a]
  - [c, b]
- logprob: -159.79270720880572
  n_clusters: 1
  n_procs: 1
  partition:
  - [a, c, b]
events:
- unique_ids: [a, c, b]
  cdr3_length: 45
  codon_positions: {j: 330, v: 288}
  d_3p_del: 1
  d_5p_del: 0
  d_gene: IGHD5-18*01
  d_per_gene_support: !!python/object/apply:collections.OrderedDict
  - - [IGHD5-18*01, 1.0]
  dj_insertion: A
  duplicates:
  - []
  - []
  - []
  fv_insertion: ''
  gl_gap_seqs: ['', '', CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTACCATATCAGTAGACCCGTCCAAGAACCAGTTCTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG]
  has_shm_indels: [false, false, true]
  in_frames: [true, true, false]
  indel_reversed_seqs: ['', '', CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCGTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGGGGTGGTTACTACTGGAGCTGGATCCGCCAGTACCCAGCGAAGTGCCTGGAGTGGGTTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTTCCATATCTGTAGACCCGTCCAAGAACCAGTTTTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG]
  input_seqs: [!!python/unicode CAGGTGCAGATGCAGGAGTCGGGCCCAGGACTATTGAAGCCTACACAGACCCTGTCCCTCACCTGCACTGTCTTTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGTACCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTGCATCTATTACAGTGGGAGCACGTACTACAACCCGTCCCTCAAGAGTCTAGTTACCATACCAGTAGACCCGTCCAAGAACCAGTTCTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGACGTGGATACAACTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGATCACCGTCTATTCAG, !!python/unicode CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCGTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGGAGTGGTTACTACTGGAACTGGATCCGCCAGTACCCAGCGAAGTGCCTGGAGTGGATTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAATTACCATATCAGTAGACTCGTCCAAGAACCATTTTTCCCTGAAGCCGAGCTCTGTAACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGGCTCTTCAG,
    !!python/unicode CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCGTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGGGGTGGTTACTACTGGAGCTGGATCCGCCAGTACCCAGCGAAGTGCCTGGAGTGGGTTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTTCCATATCTGTAGACCCGTCCAAGAACAGTTTTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG]
  invalid: false
  j_3p_del: 0
  j_5p_del: 2
  j_gene: IGHJ3*02
  j_per_gene_support: !!python/object/apply:collections.OrderedDict
  - - [IGHJ3*02, 1.0]
  jf_insertion: ''
  mut_freqs: [0.03571428571428571, 0.03571428571428571, 0.024725274725274724]
  mutated_invariants: [false, false, false]
  n_mutations: [13, 13, 9]
  naive_seq: CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCCTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGTGGTGGTTACTACTGGAGCTGGATCCGCCAGCACCCAGGGAAGGGCCTGGAGTGGATTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTACCATATCAGTAGACCCGTCCAAGAACCAGTTCTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG
  qr_gap_seqs: ['', '', !!python/unicode CAGGTGCAGCTGCAGGAGTCGGGCCCAGGACTGTTGAAGCCTTCACAGACCGTGTCCCTCACCTGCACTGTCTCTGGTGGCTCCATCAGCAGGGGTGGTTACTACTGGAGCTGGATCCGCCAGTACCCAGCGAAGTGCCTGGAGTGGGTTGGGTGCATCTATTACAGTGGGAGCACCTACTACAACCCGTCCCTCAAGAGTCGAGTTTCCATATCTGTAGACCCGTCCAAGAA.CAGTTTTCCCTGAAGCCGAGCTCTGTGACTGCCGCGGACACGGCCGTGGATTACTGTGCGAGGTGGATACAGCTATGGTTAAATGCTTTTGATATCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG]
  stops: [false, false, true]
  v_3p_del: 3
  v_5p_del: 0
  v_gene: IGHV4-31*10
  v_per_gene_support: !!python/object/apply:collections.OrderedDict
  - - [IGHV4-31*10, 1.0]
  vd_insertion: ''
```
