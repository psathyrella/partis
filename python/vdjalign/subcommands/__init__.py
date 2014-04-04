COMMANDS = ('adaptive_tsv',
            'count_mutations',
            'collapse_identical_reads',
            'deduplicate_fastq',
            'export_targets',
            'extract_best_hits',
            'align_fastq',
            'adaptive_vj_uncertainty')


def itermodules(root=__name__):
    for command in sorted(COMMANDS):
        name = command.replace('_', '-')
        mod = __import__('.'.join((root, command)), fromlist=[command])
        yield name, mod
