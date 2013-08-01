COMMANDS = ('adaptive_tsv', 'count_mutations', 'translate_bam', 'collapse_identical_reads')


def itermodules(root=__name__):
    for command in sorted(COMMANDS):
        name = command.replace('_', '-')
        mod = __import__('.'.join((root, command)), fromlist=[command])
        yield name, mod
