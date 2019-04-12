#!/usr/bin/env python
import argparse
import operator
import os
import sys
import yaml
import colored_traceback.always

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
parser.add_argument('--base-outdir', required=True)
parser.add_argument('--label', default='test')
parser.add_argument('--parse-output', action='store_true')
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

if args.workdir is None:
    args.workdir = utils.choose_random_subdir('/tmp/%s/hmms' % (os.getenv('USER', default='partis-work')))

# ----------------------------------------------------------------------------------------
def get_n_generations(ntl, tau):  # NOTE duplicates code in treeutils.get_max_lbi()
    return max(1, int(args.seq_len * tau * ntl))

# ----------------------------------------------------------------------------------------
def get_outfname(outdir):
    return '%s/vals.yaml' % outdir

# ----------------------------------------------------------------------------------------
parsed_info = {}
for lbt in args.lb_tau_list:

    gen_list = args.n_generations_list
    if gen_list is None:
        gen_list = [get_n_generations(ntl, lbt) for ntl in args.n_tau_lengths_list]
    print '   lb tau %.4f' % lbt
    print '      n gen: %s' % ' '.join(str(n) for n in gen_list)
    for n_gen in gen_list:

        # if ntl is not None:
        #     this_outdir = '%s/%s/n-tau-%.2f-lbt-%.4f' % (args.base_outdir, args.label, ntl, lbt)
        # elif n_gen is not None:
        this_outdir = '%s/%s/n_gen-%d-lbt-%.4f' % (args.base_outdir, args.label, n_gen, lbt)

        if os.path.exists(get_outfname(this_outdir)):
            if args.parse_output:
                with open(get_outfname(this_outdir)) as outfile:
                    info = yaml.load(outfile)
                if lbt not in parsed_info:
                    parsed_info[lbt] = {}
                parsed_info[lbt][n_gen] = info['max']['lbi']
            else:
                print '         output exists, skipping: %s' % get_outfname(this_outdir)
            continue

        if not os.path.exists(this_outdir):
            os.makedirs(this_outdir)

        # lbvals = treeutils.get_min_lbi(args.seq_len, args. XXX lb_tau)
        max_name, max_lbi, lbvals = treeutils.get_max_lbi(args.seq_len, lbt, n_generations=n_gen)

        with open(get_outfname(this_outdir), 'w') as outfile:
            yaml.dump({'max' : {'name' : max_name, 'lbi' : max_lbi}}, outfile)

        plotdir = this_outdir + '/plots'
        utils.prep_dir(plotdir, wildlings='*.svg')
        cmdfos = [plotting.get_lb_tree_cmd(lbvals['tree'], '%s/tree.svg' % plotdir, 'lbi', 'affinities', args.ete_path, args.workdir, metafo=lbvals)]
        utils.run_cmds(cmdfos, clean_on_success=True, shell=True, debug='print')

if args.parse_output:
    fig, ax = plotting.mpl_init()
    for lbt in sorted(parsed_info, reverse=True):
        n_gen_list, max_lbi_list = zip(*sorted(parsed_info[lbt].items(), key=operator.itemgetter(0)))
        ax.plot(n_gen_list, max_lbi_list, label='%.4f' % lbt, alpha=0.7, linewidth=4)
    plotting.mpl_finish(ax, args.base_outdir, 'tau-vs-n-gen-vs-max-lbi', log='y', xlabel='N generations', ylabel='Max LBI', leg_title='tau', leg_prop={'size' : 12}, leg_loc=(0.04, 0.67))
    
