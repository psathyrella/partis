### output formats

  * [annotation file headers](#annotation-file-headers)
  * [in-memory annotation dictionary keys](#in-memory-annotation-dictionary-keys)
  * [partition file headers](#partition-file-headers)

The `annotate` action writes a single csv file with annotations for each sequence in the input.
The `partition` action, on the other hand, writes two csv files: one with a list of the most likely partitions and their relative likelihoods, and another with annotations for each cluster in the most likely partition.
This cluster annotation file is by default written to the same path as the partition file, but with `-cluster-annotations` inserted before the suffix (you can change this with `--cluster-annotation-fname`).
These two files will eventually be combined into a single yaml file.

If you just want to view an ascii representation of the results, use the partis [`view-annotations`](subcommands.md#view-annotations) and [`view-partitions`](subcommands.md#view-partitions) actions.
If you want to directly access the csv columns, the [partition](#partition-file-headers) and [annotation](#annotation-file-headers) headers are listed below.
While by default a fairly minimal set of annotation information is written to file, many more keys are added to the in-memory dictionary when the file is read back in.
Any of these keys, together with several additional ones, can be added to the output file by setting `--extra-annotation-columns` (for all the choices, see `partis annotate --help|grep -C5 extra-annotation`).

If, on the other hand, you need to do some additional calculations, there is a lot of existing code to facilitate this.
To get you started, there are simple example parsing scripts for [annotation](../bin/example-parse-annotations.py) and [partition](../bin/example-parse-partitions.py) output files.
If more detailed calculations are necessary, you can of course just add the necessary code to the example scripts.

#### annotation file headers

The annotation csv contains the following columns by default:

|   column header        |  description
|------------------------|----------------------------------------------------------------------------
| unique_ids             |  colon-separated list of sequence identification strings
| v_gene         |  V gene in most likely annotation
| d_gene         |  D gene in most likely annotation
| j_gene         |  J gene in most likely annotation
| cdr3_length    |  CDR3 length of most likely annotation (IMGT scheme, i.e. including both codons in their entirety)
| mut_freqs      |  colon-separated list of sequence mutation frequencies
| input_seqs     |  colon-separated list of input sequences, with constant regions (fv/jf insertions) removed
| naive_seq      |  naive (unmutated ancestor) sequence corresponding to most likely annotation
| v_3p_del       |  length of V 3' deletion in most likely annotation
| d_5p_del       |  length of D 5' deletion in most likely annotation
| d_3p_del       |  length of D 3' deletion in most likely annotation
| j_5p_del       |  length of J 5' deletion in most likely annotation
| v_5p_del       |  length of an "effective" V 5' deletion in the most likely annotation, corresponding to a read which does not extend through the entire V segment
| j_3p_del       |  length of an "effective" J 3' deletion in the most likely annotation, corresponding to a read which does not extend through the entire J segment
| vd_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the V and D segments
| dj_insertion       |  sequence of nucleotides corresponding to the non-templated insertion between the D and J segments
| fv_insertion       |  constant region on the 5' side of the V (accounts for reads which extend beyond the 5' end of V)
| jf_insertion       |  constant region on the 3' side of the J (accounts for reads which extend beyond the 3' end of J)
| mutated_invariants |  true if the conserved codons corresponding to the start and end of the CDR3 code for the same amino acid as in their original germline (cyst and tryp/phen, in IMGT numbering)
| in_frames          |  true if the net effect of VDJ rearrangement and SHM indels leaves both the start and end of the CDR3 (IMGT cyst and tryp/phen) in frame with respect to the start of the germline V sequence
| stops              |  true if there's a stop codon in frame with respect to the start of the germline V sequence
| v_per_gene_support |  approximate probability supporting the top V gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| d_per_gene_support |  approximate probability supporting the top D gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| j_per_gene_support |  approximate probability supporting the top J gene matches, as a semicolon-separated list of colon-separated gene:probability pairs (approximate: monotonically related to the actual probability, but not exactly one-to-one)
| indelfos       |  colon-separated list of information on any SHM indels that were inferred in the Smith-Waterman step. Written as a literal python dict; can be read in python with `ast.literal_eval(line['indelfo'])`
| indel_reversed_seqs  |  colon-separated list of input sequences with indels "reversed" (i.e. undone), and with constant regions (fv/jf insertions) removed. Empty string if there are no indels, i.e. if it's the same as 'input_seqs'
| duplicates     |  colon-separated list of "duplicate" sequences for each sequence, i.e. sequences which, after trimming fv/jf insertions, were identical and were thus collapsed.


All columns listed as "colon-separated lists" are trivial/length one for single sequence annotation, i.e. are only length greater than one (contain actual colons) when the multi-hmm has been used for simultaneous annotation on several clonally-related sequences (typically through the cluster annotation output but partition action).
You can view a colored ascii representation of the rearrangement events with the `view-annotations` action (see below).
Additional columns (for instance, cdr3_seqs) can be specified with the `--extra-annotation-columns` option (run `partis annotate --help` to see the choices).

#### in-memory annotation dictionary keys

aligned_d_seqs
aligned_j_seqs
aligned_v_seqs
cdr3_seqs
codon_positions
d_gl_seq
d_qr_seqs
full_coding_input_seqs
full_coding_naive_seq
invalid
j_gl_seq
j_qr_seqs
lengths
regional_bounds
v_gl_seq
v_qr_seqs


#### partition file headers

The partition action writes a list of partitions, with one line for the most likely partition (with the lowest logprob), as well as a number of lines for the surrounding less-likely partitions.
It also writes the annotation for each cluster in the most likely partition to a separate file (by default `<--outfname>.replace('.csv', '-cluster-annotations.csv')`, you can change this file name with `--cluster-annotation-fname`).

|   column header  |  description
|------------------|----------------------------------------------------------------------------
| logprob          |  Total log probability of this partition
| n_clusters       |  Number of clusters (clonal families in this partition)
| partition        |  String representing the clusters, where clusters are separated by ; and sequences within clusters by :, e.g. 'a:b;c:d:e'
| n_procs          |  Number of processes which were simultaneously running for this clusterpath. In practice, final output is usually only written for n_procs = 1

