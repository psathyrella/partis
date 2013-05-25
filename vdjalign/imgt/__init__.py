from Bio import SeqIO
from pkg_resources import resource_stream

_names = ('acc', 'allele', 'species', 'function', 'region', 'start_end',
    'num_nuc', 'codon_start', 'add_5prime', 'delta_3prime', 'delta_seqerr',
    'num_aa', 'num_char', 'is_partial', 'is_revcomp')

def _populate_imgt_sequence(s):
    #The FASTA header contains 15 fields separated by '|':
    # 1. IMGT/LIGM-DB accession number(s)
    # 2. gene and allele name
    # 3. species
    # 4. functionality
    # 5. exon(s), region name(s), or extracted label(s)
    # 6. start and end positions in the IMGT/LIGM-DB accession number(s)
    # 7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
    # 8. codon start, or 'NR' (not relevant) for non coding labels and out-of-frame pseudogenes
    # 9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
    # 10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
    # 11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
    # 12. number of amino acids (AA): this field indicates that the sequence is in amino acids
    # 13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
    # 14. partial (if it is)
    # 15. reverse complementary (if it is)
    pieces = [i.strip() for i in s.description.split('|')]
    assert len(pieces) == 16, 'Unexpected split length. {0}'.format(len(pieces))
    s.annotations.update(dict(zip(_names, pieces)))
    s.id = s.annotations['allele']
    return s

def _imgt_fasta_of_data(file_name):
    with resource_stream('.'.join((__name__, 'data')), file_name) as fp:
        for s in SeqIO.parse(fp, 'fasta'):
            yield _populate_imgt_sequence(s)
