COMMANDS = ('adaptive_csv', 'count_mutations', 'translate_bam')


def itermodules(root=__name__):
    for command in sorted(COMMANDS):
        name = command.replace('_', '-')
        mod = __import__('.'.join((root, command)), fromlist=[command])
        yield name, mod
