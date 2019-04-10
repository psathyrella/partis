#!/usr/bin/env python
import argparse
import os
import sys

# ----------------------------------------------------------------------------------------
def plot_bounds(args):
    for ntl in [5]:  #range(15):
        # lbvals = treeutils.get_min_lbi(args.seq_len, args.lb_tau)
        lbvals = treeutils.get_max_lbi(args.seq_len, args.lb_tau, n_tau_lengths=ntl)

        assert args.plotdir is not None
        import plotting
        cmdfos = [plotting.get_lb_tree_cmd(lbvals['tree'], '%s/foo.svg' % args.plotdir, 'lbi', 'affinities', args.ete_path, args.workdir, metafo=lbvals)]
        utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--plotdir')
parser.add_argument('--workdir')  # default set below
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
args = parser.parse_args()

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work'))

# import plotting
# xvals, yvals = zip(*((1, 0.0034),
#                 (1, 0.0034),
#                 (2, 0.0065),
#                 (3, 0.0100),
#                 (4, 0.0130),
#                 (6, 0.0180),
#                 (7, 0.0202),
#                 (8, 0.0224),
#                 (9, 0.0242),
#                 (10, 0.0258),
#                 (12, 0.0284),
#                 (13, 0.0296),
#                 (14, 0.0307),
#                 (15, 0.0316),
#                 (16, 0.0324),
# ))
# fig, ax = plotting.mpl_init()
# ax.plot(xvals, yvals)
# plotting.mpl_finish(ax, XXX plotdir, 'foo')
# sys.exit()
plot_bounds(args)
