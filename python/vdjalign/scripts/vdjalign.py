"""
Tools for analyzing rearranged B cell receptors
"""

import argparse
import sys
import logging

from .. import subcommands, __version__ as version

DESCRIPTION = __doc__.strip()


def main(argv=sys.argv[1:]):
    arguments = parse_arguments(argv)

    loglevel = {
        0: logging.ERROR,
        1: logging.WARNING,
        2: logging.INFO,
        3: logging.DEBUG,
    }.get(arguments.verbosity, logging.DEBUG)

    # set up logging
    logging.basicConfig(file=sys.stdout, level=loglevel)

    return arguments.func(arguments)


def parse_arguments(argv):
    """
    """
    # Create the argument parser
    parser = argparse.ArgumentParser(description=DESCRIPTION)

    parser.add_argument('-V', '--version', action='version',
        version='taxit v' + version,
        help='Print the version number and exit')

    base_parser = argparse.ArgumentParser(add_help=False)
    base_parser.add_argument('-v', '--verbose',
        action='count', dest='verbosity', default=2,
        help='Increase verbosity of screen output (eg, -v is verbose, '
             '-vv more so)')
    base_parser.add_argument('-q', '--quiet',
        action='store_const', dest='verbosity', const=0,
        help='Suppress output')

    ##########################
    # Setup all sub-commands #
    ##########################

    subparsers = parser.add_subparsers(dest='subparser_name')

    # Help sub-command
    parser_help = subparsers.add_parser(
        'help', help='Detailed help for actions using `help <action>`')
    parser_help.add_argument('action', nargs=1)

    for name, mod in subcommands.itermodules():
        # set up subcommand help text. The first line of the dosctring
        # in the module is displayed as the help text in the
        # script-level help message (`script -h`). The entire
        # docstring is displayed in the help message for the
        # individual subcommand ((`script action -h`)).
        subparser = subparsers.add_parser(
            name,
            help=mod.__doc__.lstrip().split('\n', 1)[0],
            description=mod.__doc__,
            parents=[base_parser])

        mod.build_parser(subparser)

    # Determine we have called ourself (e.g. "help <action>")
    # Set arguments to display help if parameter is set
    #           *or*
    # Set arguments to perform an action with any specified options.
    arguments = parser.parse_args(argv)
    # Determine which action is in play.

    # Support help <action> by simply having this function call itself and
    # translate the arguments into something that argparse can work with.
    if arguments.subparser_name == 'help':
        return parse_arguments([str(arguments.action[0]), '-h'])

    return arguments


if __name__ == '__main__':
    main()
