COMMANDS = ('adaptive_csv', 'count_mutations')

def itermodules(root=__name__):
    for command in sorted(COMMANDS):
        yield command.replace('_', '-'), \
                __import__('.'.join((root, command)), fromlist=[command])
