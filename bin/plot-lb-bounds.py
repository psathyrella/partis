#!/usr/bin/env python
import argparse
import os
import sys
import yaml

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

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--lb-tau-list', required=True)
parser.add_argument('--n-tau-lengths-list', help='set either this or --n-generations-list')
parser.add_argument('--n-generations-list', help='set either this or --n-tau-lengths-list')
parser.add_argument('--seq-len', default=400, type=int)
parser.add_argument('--outdir')
parser.add_argument('--workdir')  # default set below
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
args = parser.parse_args()

sys.path.insert(1, args.partis_dir + '/python')
try:
    import utils
    import treeutils
    import plotting
except ImportError as e:
    print e
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.lb_tau_list = utils.get_arg_list(args.lb_tau_list, floatify=True)
args.n_tau_lengths_list = utils.get_arg_list(args.n_tau_lengths_list, floatify=True)
args.n_generations_list = utils.get_arg_list(args.n_generations_list, intify=True)
if [args.n_tau_lengths_list, args.n_generations_list].count(None) != 1:
    raise Exception('have to set exactly one of --n-tau-lengths, --n-generations')
if args.n_tau_lengths_list is None:
    args.n_tau_lengths_list = [None for _ in args.n_generations_list]
if args.n_generations_list is None:
    args.n_generations_list = [None for _ in args.n_tau_lengths_list]

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

# ----------------------------------------------------------------------------------------
for lbt in args.lb_tau_list:
    print '   lb tau %.4f' % lbt
    if args.n_tau_lengths_list is not None:
        print '         tau lengths'
    elif args.n_generations_list is not None:
        print '         generations'
    for ntl, n_gen in zip(args.n_tau_lengths_list, args.n_generations_list):
        print '           %s' % (('%.2f' % ntl) if ntl is not None else ('%d' % n_gen))
        # lbvals = treeutils.get_min_lbi(args.seq_len, args. XXX lb_tau)
        # max_name, max_lbi, lbvals = treeutils.get_max_lbi(args.seq_len, lbt, n_tau_lengths=ntl, n_generations=n_gen)

# outdir=$fs/partis/lb-tau-optimization/test-2d/n_gen-$n_gen-lbt-$lbt
	# limit_procs python
    # if args.plotdir is not None:
    #     utils.prep_dir(args.plotdir, wildlings='*.svg', allow_other_files=True)
    #     cmdfos = [plotting.get_lb_tree_cmd(lbvals['tree'], '%s/tree.svg' % args.plotdir, 'lbi', 'affinities', args.ete_path, args.workdir, metafo=lbvals)]
    #     utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

    # if args.outfname is not None:
    #     if not os.path.exists(os.path.dirname(args.outfname)):
    #         os.makedirs(os.path.dirname(args.outfname))
    #     with open(args.outfname, 'w') as outfile:
    #         yaml.dump({'max' : {'name' : max_name, 'lbi' : max_lbi}})
