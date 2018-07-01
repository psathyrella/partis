### output formats

  * [output file overview](#output-file-overview)
  * [default annotation file headers](#default-annotation-file-headers)


All output is written to a unified yaml file (for documentation on the old csv formats, see [here](https://github.com/psathyrella/partis/blob/0f1d8969af888a343d04524c3b8f21075896d8e4/manual.md)).
The `annotate` action writes annotations for each single sequence in the input.
The `partition` action, on the other hand, writes a list of the most likely partitions and their relative likelihoods, as well as annotations for each cluster in the most likely partition.

If you just want to print the results to the terminal, use the partis [`view-output`](subcommands.md#view-output) action.
While by default a fairly minimal set of annotation information is written to file, many more keys are present in the dictionary in memory (see below).
Any of these keys, together with several additional ones, can be added to the output file by setting `--extra-annotation-columns` (for all the choices, see below, or run `partis annotate --help|grep -C5 extra-annotation`).

An example parsing script can be found [here](../bin/example-parse-output.py).

#### output file overview

The yaml output file contains four top-level headers.
To keep track of the output file version, we have `version-info`:

`version-info: {partis-yaml: 0.1}`

The set of germline genes corresponding to the annotations is in `germline-info`:

```
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
```

The annotations are stored under `events`, for instance here is an annotation for three clonally-related sequences:

```
events:
- cdr3_length: 45
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
  unique_ids: [a, c, b]
  v_3p_del: 3
  v_5p_del: 0
  v_gene: IGHV4-31*10
  v_per_gene_support: !!python/object/apply:collections.OrderedDict
  - - [IGHV4-31*10, 1.0]
  vd_insertion: ''
```

and finally, the list of partitions is stored under `partitions`:

```
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
```

#### default annotation file headers

|   name         |  description
|----------------|----------------------------------------------------------------------------
| unique_ids     |  colon-separated list of sequence identification strings
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
| codon_positions    |  zero-indexed V and J conserved codon positions in the indel-reversed sequence[s]
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

while the following (together with any in-memory key) can be added using `--extra-annotation-columns`:

|   name                  |  description
|-------------------------|----------------------------------------------------------------------------
| cdr3_seqs				  |  nucleotide CDR3 sequence, including bounding conserved codons
| full_coding_naive_seq	  |  in cases where the input reads do not extend through the entire V and J regions, the input_seqs and naive_seq keys will also not cover the whole coding regions. In such cases full_coding_naive_seq and full_coding_input_seqs can be used to tack on the missing bits.
| full_coding_input_seqs  |  see full_coding_naive_seq

All columns listed as "colon-separated lists" are trivial/length one for single sequence annotation, i.e. are only length greater than one (contain actual colons) when the multi-hmm has been used for simultaneous annotation on several clonally-related sequences (typically in the cluster annotation file from partitioning, but can also be set using `--n-simultaneous-seqs` and `--simultaneous-true-clonal-seqs`).

deprecated keys (only present in old files):

|   column header        |  description
|------------------------|----------------------------------------------------------------------------
| indelfos       |  colon-separated list of information on any SHM indels that were inferred in the Smith-Waterman step

#### in-memory annotation dictionary keys

Extra information available in the dictionary in memory, but not written to disc by default (can be written by setting `--extra-annotation-columns`):

|   key                   |  value
|-------------------------|----------------------------------------------------------------------------
| codon_positions		  |  zero-indexed indel-reversed-sequence positions of the conserved cyst and tryp/phen codons, e.g. `{'v': 285, 'j': 336}`
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

#### partition file headers

The partition action writes a list of partitions, with one line for the most likely partition (with the lowest logprob), as well as a number of lines for the surrounding less-likely partitions.
It also writes the annotation for each cluster in the most likely partition to a separate file (by default `<--outfname>.replace('.csv', '-cluster-annotations.csv')`, you can change this file name with `--cluster-annotation-fname`).

|   column header  |  description
|------------------|----------------------------------------------------------------------------
| logprob          |  Total log probability of this partition
| n_clusters       |  Number of clusters (clonal families in this partition)
| partition        |  String representing the clusters, where clusters are separated by ; and sequences within clusters by :, e.g. 'a:b;c:d:e'
| n_procs          |  Number of processes which were simultaneously running for this clusterpath. In practice, final output is usually only written for n_procs = 1

