#!/usr/bin/env python
import os
import argparse
import sys
import random
import re
from collections import OrderedDict
import time
import csv
from subprocess import check_call, Popen
sys.path.insert(1, './python')

from humans import humans
import plotting
import seqfileopener
import utils
from clusterplot import ClusterPlot
from glomerator import Glomerator

parser = argparse.ArgumentParser()
parser.add_argument('-b', action='store_true')
parser.add_argument('--dataset', choices=['stanford', 'adaptive'], default='adaptive')
parser.add_argument('--only-run')  # colon-separated list of human,subset pairs to run, e.g. A,3:C,8
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'write-plots']
parser.add_argument('--actions', required=True)  #default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)

legends = {'vollmers-0.9' : 'VJ CDR3 0.9',
           'partition partis' : 'full partis',
           'naive-hamming-partition partis' : 'naive partis',
           'vsearch-partition partis' : 'vsearch partis',
           'changeo' : 'Change-O'
           }


if args.dataset == 'stanford':
    datadir = '/shared/silo_researcher/Matsen_F/MatsenGrp/data/stanford-lineage/2014-11-17-vollmers'
    files = [ datadir + '/' + f for f in os.listdir(datadir)]
elif args.dataset == 'adaptive':
    datadirs = [ '/shared/silo_researcher/Matsen_F/MatsenGrp/data/bcr/output_sw/' + h for h in humans['adaptive'] ]
    files = []
    for datadir in datadirs:
        files += [ datadir + '/' + fname for fname in os.listdir(datadir) if '-M_merged.tsv.bz2' in fname ]  # if you switch to naive (N), be careful 'cause A is split in pieces

# ----------------------------------------------------------------------------------------
def leafmutstr(n_leaves, mut_mult):
    return 'simu-' + str(n_leaves) + '-leaves-' + str(mut_mult) + '-mutate'

# ----------------------------------------------------------------------------------------
def write_latex_table(adj_mis):
    for mut_mult in mutation_multipliers:
  # \textbf{multiplier} & \textbf{program}  &  5     &  10     &  25     &  50   \\
  # \multirow{3}{*}{$\times 1$} & partis   &     0.94   &     0.93   &     0.89   &     0.81 \\
        for n_leaves in n_leaf_list:
            if n_leaves == n_leaf_list[0]:
                print '\\textbf{multiplier} & \\textbf{program}  ',
            print ' & %s  ' % n_leaves,
        print '\\\\'
        print '\\hline'
        iname = 0
        if 5 not in adj_mis:  # TODO remove
            adj_mis[5] = {}
        if mut_mult not in adj_mis[5]:  # TODO remove
            adj_mis[5][mut_mult] = {}
            for name in adj_mis[10]['1']:
                adj_mis[5][mut_mult][name] = -1
        for name in adj_mis[n_leaf_list[0]][mutation_multipliers[0]]:
            if name == 'vollmers-0.5':
                continue
            if iname == 0:
                print '\\multirow{3}{*}{$\\times %s$} & ' % mut_mult,
            else:
                print '& ',
            print '%25s' % legends.get(name, name),
            iname += 1
            for n_leaves in n_leaf_list:
                if n_leaves == 5 and mut_mult == '1' or n_leaves == 25 and mut_mult == '4':
                    print '  &    %5.2f' % -1,
                else:
                    print '  &    %5.2f' % adj_mis[n_leaves][mut_mult][name],
            print '\\\\'

# ----------------------------------------------------------------------------------------
def parse_partis(action, these_hists, these_adj_mis, simfname, histfname):
    args.infnames = [simfname.replace('.csv', '-' + action + '.csv'), ]  # NOTE make sure not to add any args here that conflict with the real command line args
    args.is_data = False
    args.use_all_steps = False
    args.normalize_axes = []
    args.xbounds, args.adjmi_bounds, args.logprob_bounds = None, None, None
    cplot = ClusterPlot(args)
    cplot.tmp_cluster_size_hist.write(histfname)
    these_hists[action + ' partis'] = cplot.tmp_cluster_size_hist
    these_adj_mis[action + ' partis'] = cplot.adj_mi_at_max_logprob

# ----------------------------------------------------------------------------------------
def parse_vollmers(these_hists, these_adj_mis, simfname, histfbase, true_histfname):
    vollmers_fname = simfname.replace('.csv', '-run-viterbi.csv')
    with open(vollmers_fname) as vfile:
        vreader = csv.DictReader(vfile)
        for line in vreader:
            vhist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['clusters']))
            vhist.write(histfbase + line['threshold'] + '-hist.csv')
            these_hists['vollmers-' + line['threshold']] = vhist
            these_adj_mis['vollmers-' + line['threshold']] = float(line['adj_mi'])

            truehist = plotting.get_cluster_size_hist(utils.get_partition_from_str(line['true_clusters']))  # true partition is also written here, for lack of a better place (note that it's of course the same for all thresholds)
            truehist.write(true_histfname)  # will overwite itself a few times
            these_hists['true'] = truehist

# ----------------------------------------------------------------------------------------
def parse_changeo(these_hists, these_adj_mis, simfname, simfbase):
    input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
    
    indir = fsdir.replace('/partis-dev/_output', '/changeo')
    infname = indir + '/' + simfbase.replace('-', '_') + '_db-pass_parse-select_clone-pass.tab'
    
    id_clusters = {}  # map from cluster id to list of seq ids
    with open(infname) as chfile:
        reader = csv.DictReader(chfile, delimiter='\t')
        for line in reader:
            clid = line['CLONE']
            uid = line['SEQUENCE_ID']
            if clid not in id_clusters:
                id_clusters[clid] = []
            id_clusters[clid].append(uid)
    
    partition = [ids for ids in id_clusters.values()]
    these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
    adj_mi_fname = indir + '/' + simfbase.replace('-', '_') + '-adj_mi.csv'
    with open(adj_mi_fname) as adj_mi_file:
        reader = csv.DictReader(adj_mi_file, fieldnames=['adj_mi'])
        for line in reader:
            these_adj_mis['changeo'] = float(line['adj_mi'])
            break

# ----------------------------------------------------------------------------------------
def write_all_plot_csvs(label):
    hists, adj_mis = {}, {}
    for n_leaves in n_leaf_list:
        for mut_mult in mutation_multipliers:
            if n_leaves == 5 and mut_mult == '1':
                continue
            if n_leaves == 25 and mut_mult == '4':
                continue
            write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis)

    write_latex_table(adj_mis)

# ----------------------------------------------------------------------------------------
def write_each_plot_csvs(label, n_leaves, mut_mult, hists, adj_mis):
    if n_leaves not in hists:
        hists[n_leaves] = {}
        adj_mis[n_leaves] = {}
    if mut_mult not in hists[n_leaves]:
        hists[n_leaves][mut_mult] = OrderedDict()
        adj_mis[n_leaves][mut_mult] = OrderedDict()
    these_hists = hists[n_leaves][mut_mult]
    these_adj_mis = adj_mis[n_leaves][mut_mult]

    simfname = fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code
    simfbase = leafmutstr(n_leaves, mut_mult)
    csvdir = os.path.dirname(simfname) + '/plots'

    # then vollmers annotation (and true hists)
    parse_vollmers(these_hists, these_adj_mis, simfname, csvdir + '/' + simfbase + '-vollmers-', csvdir + '/' + simfbase + '-true-hist.csv')

    # then changeo
    parse_changeo(these_hists, these_adj_mis, simfname, simfbase)

    # first do partis stuff
    for ptype in ['vsearch-', 'naive-hamming-', '']:
        parse_partis(ptype + 'partition', these_hists, these_adj_mis, simfname, csvdir + '/' + simfbase + '-' + ptype + 'partis-hist.csv')

    plotdir = os.getenv('www') + '/partis/clustering/' + label
    plotting.plot_cluster_size_hists(plotdir + '/plots/' + simfbase + '.svg', these_hists, title='%s leaves, %sx mutation' % (n_leaves, mut_mult), legends=legends, xmax=n_leaves*3.01)
    check_call(['./bin/makeHtml', plotdir, '3', 'null', 'svg'])
    check_call(['./bin/permissify-www', plotdir])

# ----------------------------------------------------------------------------------------
def execute(action, label, datafname, n_leaves=None, mut_mult=None):
    real_action = action
    if 'partition' in action:
        real_action = 'partition'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + real_action
    if n_leaves is not None:
        simfname = fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code
    extras = []

    def output_exists(outfname):
        if os.path.exists(outfname):
            print '                      partition output exists, skipping (%s)' % outfname
            return True
        else:
            return False

    if action == 'cache-data-parameters':
        if os.path.exists(fsdir + '/' + label + '/data'):
            print '                      data parametes exist in %s, skipping ' % os.path.dirname(simfname)
            return
        cmd += ' --datafname ' + datafname
        extras += ['--n-max-queries', + n_data_to_cache]
        n_procs = n_data_to_cache / 500
    elif action == 'simulate':
        if os.path.exists(simfname):
            print '                      simfile exists, skipping (%s)' % simfname
            return
        cmd += ' --simfname ' + simfname
        extras += ['--n-sim-events', int(float(n_sim_seqs) / n_leaves)]
        extras += ['--n-leaves', n_leaves, '--mutation-multiplier', mut_mult]
        n_procs = 10
    elif action == 'cache-simu-parameters':
        if os.path.exists(simfname.replace('.csv', '')):
            print '                      simu parameters exist in %s, skipping ' % simfname.replace('.csv', '')
            return
        cmd += ' --simfname ' + simfname
        n_procs = 20
    elif action == 'partition':
        if os.path.exists(simfname.replace('.csv', '-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + simfname.replace('.csv', '-' + action + '.csv')
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 15  # something like 15 seqs/process to start with
    elif action == 'naive-hamming-partition':
        if os.path.exists(simfname.replace('.csv', '-naive-hamming-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-naive-hamming-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + simfname.replace('.csv', '-' + action + '.csv')
        extras += ['--n-max-queries', n_to_partition, '--auto-hamming-fraction-bounds']
        n_procs = n_to_partition / 30
    elif action == 'vsearch-partition':
        outfname = simfname.replace('.csv', '-' + action + '.csv')
        if output_exists(outfname):
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', n_to_partition, '--naive-vsearch']
        n_procs = n_to_partition / 100  # only used for ighutil step
    elif action == 'run-viterbi':
        if os.path.exists(simfname.replace('.csv', '-run-viterbi.csv')):
            print '                      vollmers output exists, skipping (%s)' % simfname.replace('.csv', '-run-viterbi.csv')
            return
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = n_to_partition / 50
    elif action == 'run-changeo':
        changeodir = '/home/dralph/work/changeo/changeo'
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo'
        _simfbase = leafmutstr(n_leaves, mut_mult).replace('-', '_')
        if os.path.isdir(changeo_fsdir + '/' + _simfbase):
            print '                      changeo txz dir exists %s' % changeo_fsdir + '/' + _simfbase
            return
        tar_cmd = 'mkdir ' + changeo_fsdir + '/' + _simfbase + ';'
        tar_cmd += ' tar Jxvf ' + changeo_fsdir + '/' + _simfbase + '.txz --exclude=\'IMGT_HighV-QUEST_individual_files_folder/*\' -C ' + changeo_fsdir + '/' + _simfbase
        # 'tar Jxvf /fh/fast/matsen_e/dralph/work/changeo/simu_10_leaves_1_mutate.txz --exclude='IMGT_HighV-QUEST_individual_files_folder/*' -C foop/'
        # print tar_cmd
        # sys.exit()
        check_call(tar_cmd, shell=True)

        def run(cmdstr):
            print 'RUN %s' % cmdstr
            check_call(cmdstr.split())

        
        cmd = changeodir + '/MakeDb.py imgt -i ' + changeo_fsdir + '/' + _simfbase + ' -s ' + changeo_fsdir + '/' + _simfbase.replace('_', '-') + '.fasta'
        run(cmd)
        cmd = changeodir + '/ParseDb.py select -d ' + changeo_fsdir + '/' + _simfbase + '_db-pass.tab -f FUNCTIONAL -u T'
        run(cmd)
        cmd = changeodir + '/DefineClones.py bygroup -d ' + changeo_fsdir + '/' + _simfbase + '_db-pass_parse-select.tab --act first --model m1n --dist 7'
        run(cmd)

        # read changeo's output and toss it into a csv
        input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
        infname = changeo_fsdir + '/' + _simfbase.replace('-', '_') + '_db-pass_parse-select_clone-pass.tab'
        id_clusters = {}  # map from cluster id to list of seq ids
        with open(infname) as chfile:
            reader = csv.DictReader(chfile, delimiter='\t')
            for line in reader:
                clid = line['CLONE']
                uid = line['SEQUENCE_ID']
                if clid not in id_clusters:
                    id_clusters[clid] = []
                id_clusters[clid].append(uid)
        
        partition = [ids for ids in id_clusters.values()]
        # these_hists['changeo'] = plotting.get_cluster_size_hist(partition)
        adj_mi_fname = changeo_fsdir + '/' + _simfbase.replace('-', '_') + '-adj_mi.csv'
        with open(adj_mi_fname, 'w') as adj_mi_file:
            writer = csv.DictWriter(adj_mi_file, ['adj_mi'])
            print 'calcing...'
            writer.writerow({'adj_mi' : utils.mutual_information(partition, reco_info, debug=True)})
            print '  done'
        return
    else:
        raise Exception('bad action %s' % action)

    cmd +=  ' --plotdir ' + os.getenv('www') + '/partis'
    n_proc_str = str(n_procs)
    if n_procs > 10:
        n_proc_str += ':10'
        extras += ['--slurm', '--workdir', fsdir + '/_tmp/' + str(random.randint(0,99999))]
        
    cmd += ' --n-procs ' + n_proc_str
    # extras += ['--print-partitions', ]

    cmd += utils.get_extra_str(extras)
    print '   ' + cmd
    # return
    # check_call(cmd.split())
    # return
    logbase = os.path.dirname(simfname) + '/_logs/' + os.path.basename(simfname).replace('.csv', '') + '-' + action
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    time.sleep(5)

# ----------------------------------------------------------------------------------------
n_to_partition = 5000
n_data_to_cache = 50000
mutation_multipliers = ['4']  #['1', '4']
n_leaf_list = [5, 10, 25, 50]
n_sim_seqs = 10000
fsdir = '/fh/fast/matsen_e/' + os.getenv('USER') + '/work/partis-dev/_output'
procs = []

# ----------------------------------------------------------------------------------------
for datafname in files:
    if args.dataset == 'stanford':
	    human = os.path.basename(datafname).replace('_Lineages.fasta', '')
    elif args.dataset == 'adaptive':
	    human = re.findall('[ABC]', datafname)[0]
    print 'run', human
    label = human

    if args.only_run is not None and human not in args.only_run:
        continue

    performance_fname = fsdir + '/' + label + '/partition-performance.csv'
    performance_file_header = ['method', 'mut_mult', 'n_leaves', 'adj_mi', 'n_clusters', 'n_true_clusters']

    for action in args.actions:
        if action == 'write-plots':
            continue
        print '      ----> ', action
        if action == 'cache-data-parameters':
            execute(action, label, datafname)
            continue
        for n_leaves in n_leaf_list:
            print '  ----> ', n_leaves, ' leaves'
            for mut_mult in mutation_multipliers:
                if n_leaves == 5 and mut_mult == '1':
                    continue
                if n_leaves == 25 and mut_mult == '4':
                    continue
                print '         ----> mutate', mut_mult
                execute(action, label, datafname, n_leaves, mut_mult)
                # sys.exit()

    if 'write-plots' in args.actions:
        write_all_plot_csvs(label)

    break
