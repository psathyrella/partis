#!/usr/bin/env python
import sys
import colored_traceback.always
import os
import yaml
import argparse

sys.path.insert(1, './python')
import utils
import treeutils

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('action', choices=['train', 'test'])
parser.add_argument('--infname', required=True)
parser.add_argument('--base-plotdir', required=True)
parser.add_argument('--lb-tau', required=True, type=float)
parser.add_argument('--dtr-path', required=True)
parser.add_argument('--dtr-cfg')
parser.add_argument('--only-csv-plots', action='store_true')
# tree plots turned off in the treeutils fcn atm
# parser.add_argument('--ete-path', default=('/home/%s/anaconda_ete/bin' % os.getenv('USER')) if os.getenv('USER') is not None else None)
# parser.add_argument('--workdir')  # only required to make ete trees
args = parser.parse_args()

_, true_lines, _ = utils.read_output(args.infname) #, n_max_queries=10000)
treeutils.calculate_non_lb_tree_metrics('dtr', true_lines, base_plotdir=args.base_plotdir, train_dtr=args.action=='train', dtr_path=args.dtr_path, dtr_cfg=args.dtr_cfg,
                                        only_csv=args.only_csv_plots, lb_tau=args.lb_tau)  # ete_path=args.ete_path, workdir=args.workdir, 
