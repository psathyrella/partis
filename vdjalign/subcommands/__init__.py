COMMANDS = ('adaptive_csv',)

def itermodules(root=__name__):
    for command in sorted(COMMANDS):
        yield command.replace('_', '-'), \
                __import__('.'.join((root, command)), fromlist=[command])
