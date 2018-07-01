### output format

  * [output file overview](#output-file-overview)
  * [description of keys](#description-of-keys)
  * [output file example](#output-file-example)


All output is written to a unified yaml file (for documentation on the old csv formats, see [here](https://github.com/psathyrella/partis/blob/9d78600ba6b7d511825f27724de1f8e937a5ada3/docs/output-formats.md)).
The `annotate` action writes annotations for each single sequence in the input.
The `partition` action writes a list of the most likely partitions and their relative likelihoods, as well as annotations for each cluster in the most likely partition.

If you want to print the results of existing output files to the terminal, use the partis [`view-output`](subcommands.md#view-output) action.
While by default a fairly minimal set of annotation information is written to file, many more keys are present in the dictionary in memory (see below).
Any of these keys, together with several additional ones, can be added to the output file by setting `--extra-annotation-columns key_a:key_b` (for all the choices, see below, or run `partis annotate --help|grep -C5 extra-annotation`).

An example parsing script can be found [here](../bin/example-parse-output.py).

For more information on all options, run `partis <action> --help`.

#### output file overview

The yaml output file contains four top-level headers: 

|   name         |  description
|----------------|----------------------------------------------------------------------------
|  version-info  |  output file format version
|  germline-info |  germline sequence, names, and conserved codon positions
|  events        |  list of annotations for each rearrangement event (i.e. group of clonally-related sequences)
|  partitions    |  list of partitions, including the most likely partition (only set if running the partition action)

#### description of keys

The following keys are written to output by default:

|   name         |  description
|----------------|----------------------------------------------------------------------------
| unique_ids     |  colon-separated list of sequence identification strings
| reco_id        |  simulation only: hash of rearrangement parameters that is the same for all clonally-related sequences
| v_gene         |  V gene in most likely annotation
| d_gene         |  see v_gene
| j_gene         |  see v_gene
| cdr3_length    |  CDR3 length of most likely annotation (IMGT scheme, i.e. including both codons in their entirety)
| mut_freqs      |  colon-separated list of sequence mutation frequencies
| input_seqs     |  colon-separated list of input sequences, with constant regions (fv/jf insertions) removed
| naive_seq      |  naive (unmutated ancestor) sequence corresponding to most likely annotation
| v_3p_del       |  length of V 3' deletion in most likely annotation
| d_5p_del       |  see v_3p_del
| d_3p_del       |  see v_3p_del
| j_5p_del       |  see v_3p_del
| v_5p_del       |  length of an "effective" V 5' deletion in the most likely annotation, corresponding to a read which does not extend through the entire V segment
| j_3p_del       |  see v_5p_del
| vd_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the V and D segments
| dj_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the D and J segments
| fv_insertion       |  constant region on the 5' side of the V (accounts for reads which extend beyond the 5' end of V)
| jf_insertion       |  constant region on the 3' side of the J (accounts for reads which extend beyond the 3' end of J)
| codon_positions	 |  zero-indexed indel-reversed-sequence positions of the conserved cyst and tryp/phen codons, e.g. `{'v': 285, 'j': 336}`
| mutated_invariants |  true if the conserved codons corresponding to the start and end of the CDR3 code for the same amino acid as in their original germline (cyst and tryp/phen, in IMGT numbering)
| in_frames          |  true if the net effect of VDJ rearrangement and SHM indels leaves both the start and end of the CDR3 (IMGT cyst and tryp/phen) in frame with respect to the start of the germline V sequence
| stops              |  true if there's a stop codon in frame with respect to the start of the germline V sequence
| v_per_gene_support |  approximate probability supporting the top V gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| d_per_gene_support |  see v_per_gene_support
| j_per_gene_support |  see v_per_gene_support
| indel_reversed_seqs  |  colon-separated list of input sequences with indels "reversed" (i.e. undone), and with constant regions (fv/jf insertions) removed. Empty string if there are no indels, i.e. if it's the same as 'input_seqs'
| gl_gap_seqs        |  colon-separated list of germline sequences with gaps at shm indel positions (alignment matches qr_gap_seqs)
| qr_gap_seqs        |  colon-separated list of query sequences with gaps at shm indel positions (alignment matches gl_gap_seqs)
| duplicates     |  colon-separated list of "duplicate" sequences for each sequence, i.e. sequences which, after trimming fv/jf insertions, were identical and were thus collapsed.

While the following keys are available in the dictionary in memory, but not written to disk by default (can be written by setting `--extra-annotation-columns key_a:key_b`):

|   key                   |  value
|-------------------------|----------------------------------------------------------------------------
| v_gl_seq				  |  portion of v germline gene aligned to the indel-reversed sequence (i.e. with 5p and 3p deletions removed). Colon-separated list.
| d_gl_seq				  |  see v_gl_seq
| j_gl_seq				  |  see v_gl_seq
| v_qr_seqs               |  portion of indel-reversed sequence aligned to the v region. Colon-separated list.
| d_qr_seqs				  |  see v_qr_seqs
| j_qr_seqs				  |  see v_qr_seqs
| lengths				  |  lengths aligned to each of the v, d, and j regions, e.g. `{'j': 48, 'd': 26, 'v': 296}`
| regional_bounds		  |  indices corresponding to the boundaries of the v, d, and j regions (python slice conventions), e.g. `{'j': (322, 370), 'd': (296, 322), 'v': (0, 296)}`
| aligned_v_seqs		  |  colon-separated list of indel-reversed sequences aligned to germline sequences given by `--aligned-germline-fname`. Only used for presto output
| aligned_d_seqs          |  see aligned_v_seqs
| aligned_j_seqs		  |  see aligned_v_seqs
| invalid				  |  indicates an invalid rearrangement event

The following keys can also be added to the output file using `--extra-annotation-columns key_a:key_b`:

|   name                  |  description
|-------------------------|----------------------------------------------------------------------------
| cdr3_seqs				  |  nucleotide CDR3 sequence, including bounding conserved codons
| full_coding_naive_seq	  |  in cases where the input reads do not extend through the entire V and J regions, the input_seqs and naive_seq keys will also not cover the whole coding regions. In such cases full_coding_naive_seq and full_coding_input_seqs can be used to tack on the missing bits.
| full_coding_input_seqs  |  see full_coding_naive_seq


Partitioning results in a list of partitions, with one line for the most likely partition (the one with the lowest logprob), as well as a number of lines for the surrounding less-likely partitions.
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
Note that while theis examples is in full yaml (since it's more human readble), in practice we read and write output files using the json subset of yaml because it's much faster.

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
