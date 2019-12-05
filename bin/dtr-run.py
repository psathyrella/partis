#!/usr/bin/env python
import sys
import colored_traceback.always
import os
import yaml
import argparse
import numpy

sys.path.insert(1, './python')
import utils
import treeutils

# NOTE not just to run dtr, also to run other non-lb metrics

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--infname', required=True)
parser.add_argument('--base-plotdir', required=True)
parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--action', choices=['train', 'test'])
parser.add_argument('--dtr-path')
parser.add_argument('--metric-method', default='dtr')
parser.add_argument('--dtr-cfg')
parser.add_argument('--only-csv-plots', action='store_true')
parser.add_argument('--n-max-queries', type=int)
parser.add_argument('--max-family-size', type=int, help='subset each family down to this size before passing to treeutils')
# tree plots turned off in the treeutils fcn atm
# parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
# parser.add_argument('--workdir')  # only required to make ete trees
args = parser.parse_args()

if args.n_max_queries is not None:
    print '    --n-max-queries set to %d' % args.n_max_queries
glfo, true_lines, _ = utils.read_output(args.infname, n_max_queries=args.n_max_queries)

# numpy.random.seed(1)
if args.max_family_size is not None:
    for line in [l for l in true_lines if len(l['unique_ids']) > args.max_family_size]:
        iseqs_to_keep = numpy.random.choice(range(len(line['unique_ids'])), args.max_family_size)
        utils.restrict_to_iseqs(line, iseqs_to_keep, glfo)

train_dtr = args.metric_method == 'dtr' and args.action == 'train'
treeutils.calculate_non_lb_tree_metrics(args.metric_method, true_lines, base_plotdir=args.base_plotdir, train_dtr=train_dtr, dtr_path=args.dtr_path, dtr_cfg=args.dtr_cfg,
                                        only_csv=args.only_csv_plots, lb_tau=args.lb_tau, min_cluster_size=args.max_family_size)  # ete_path=args.ete_path, workdir=args.workdir,
