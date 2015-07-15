#!/usr/bin/env python
import os
import argparse
import glob
import sys
import random
import re
from collections import OrderedDict
import time
import csv
from subprocess import check_call, Popen, check_output
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
parser.add_argument('--subset', type=int)
parser.add_argument('--n-subsets', type=int)
parser.add_argument('--istartstop')  # NOTE usual zero indexing
all_actions = ['cache-data-parameters', 'simulate', 'cache-simu-parameters', 'partition', 'naive-hamming-partition', 'vsearch-partition', 'run-viterbi', 'run-changeo', 'write-plots', 'compare-sample-sizes']
parser.add_argument('--actions', required=True)  #default=':'.join(all_actions))
args = parser.parse_args()
args.only_run = utils.get_arg_list(args.only_run)
args.actions = utils.get_arg_list(args.actions)
args.istartstop = utils.get_arg_list(args.istartstop, intify=True)

assert args.subset is None or args.istartstop is None  # dosn't make sense to set both of them

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
def get_simfname(label, n_leaves, mut_mult):
    return fsdir + '/' + label + '/' + leafmutstr(n_leaves, mut_mult) + '.csv'  # NOTE duplicate code

# ----------------------------------------------------------------------------------------
def generate_incorrect_partition(true_partition, n_misassigned, error_type, debug=False):
    new_partition = list(true_partition)
    if debug:
        print '  before', new_partition
    for _ in range(n_misassigned):
        iclust = random.randint(0, len(new_partition) - 1)  # choose a cluster from which to remove a sequence
        iseq = random.randint(0, len(new_partition[iclust]) - 1)  # and choose the sequence to remove from this cluster
        uid = new_partition[iclust][iseq]
        new_partition[iclust].remove(uid)  # remove it
        if [] in new_partition:
            new_partition.remove([])
        if error_type == 'singletons':  # put the sequence in a cluster by itself
            new_partition.append([uid, ])
            if debug:
                print '    %s: %d --> singleton' % (uid, iclust)
        elif error_type == 'reassign':  # choose a different cluster to add it to
            inewclust = iclust
            while inewclust == iclust:
                inewclust = random.randint(0, len(new_partition) - 1)
            new_partition[inewclust].append(uid)
            if debug:
                print '    %s: %d --> %d' % (uid, iclust, inewclust)
        else:
            assert False
    if debug:
        print '  after', new_partition
    return new_partition

# ----------------------------------------------------------------------------------------
def write_latex_table(adj_mis):
    for mut_mult in mutation_multipliers:
        for n_leaves in n_leaf_list:
            if n_leaves == n_leaf_list[0]:
                print '\\textbf{multiplier} & \\textbf{program}  ',
            print ' & %s  ' % n_leaves,
        print '\\\\'
        print '\\hline'
        iname = 0
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
    # indir = indir.replace('changeo', 'changeo.bak')
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

    simfname = get_simfname(label, n_leaves, mut_mult)
    if args.subset is not None:  # TODO not yet functional
        assert False
        subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
    if args.istartstop is not None:  # TODO not yet functional
        assert False
        # subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
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

    # for name, adj_mi in these_adj_mis.items():
    #     with open(os.path.dirname(simfname) + '/adj_mis/' + name.replace(' ', '_') + '.csv', 'w') as adjmifile:
    #         adjmifile.write(adj_mi + '\n')

# ----------------------------------------------------------------------------------------
def compare_subsets(label, n_leaves, mut_mult):
    adj_mis = {}
    # for fname in glob.glob(os.path.dirname(get_simfname(label, n_leaves, mut_mult)) + '/adj_mis/*.csv'):
    # for fname in glob.glob('_tmp/adj_mis/*.csv'):
    #     with open(fname

# ----------------------------------------------------------------------------------------
def compare_sample_sizes(label, n_leaves, mut_mult):
    n_reps = 3
    misassign_fractions = [0.1, ]
    simfname = get_simfname(label, n_leaves, mut_mult)
    input_info, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)
    uid_list = input_info.keys()
    for nseqs in [100, ]:
        for irep in range(n_reps):  # repeat <nreps> times
            istart = irep * nseqs
            istop = istart + nseqs
            uids = uid_list[istart : istop]
            true_partition = utils.get_true_clusters(uids, reco_info).values()
            for mfrac in misassign_fractions:
                n_misassigned = int(mfrac * nseqs)
                new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type='reassign')
                # new_partition = generate_incorrect_partition(true_partition, n_misassigned, error_type='singletons')
                print utils.mutual_information(new_partition, reco_info, debug=True)

# ----------------------------------------------------------------------------------------
def execute(action, label, datafname, n_leaves=None, mut_mult=None):
    real_action = action
    if 'partition' in action:
        real_action = 'partition'
    cmd = './bin/run-driver.py --label ' + label + ' --action ' + real_action
    if n_leaves is not None:
        simfname = get_simfname(label, n_leaves, mut_mult)
        if args.subset is not None:
            ntot = int(check_output(['wc', '-l', simfname]).split()[0]) - 1
            subsimfname = simfname.replace(label + '/', label + '/subset-' + str(args.subset) + '/')
            if os.path.exists(subsimfname):
                print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d / %d' % (args.subset, args.n_subsets)
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                n_per_subset = int(float(ntot) / args.n_subsets)  # NOTE ignores remainders, i.e. last few sequences
                istart = args.subset * n_per_subset + 2  # NOTE sed indexing (one-indexed with inclusive bounds). Also note extra +1 to avoid header
                istop = istart + n_per_subset - 1
                check_call('sed -n \'' + str(istart) + ',' + str(istop) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)
            simfname = subsimfname
        if args.istartstop is not None:
            subsimfname = simfname.replace(label + '/', label + '/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
            if os.path.exists(subsimfname):
                print '      subset file exists %s' % subsimfname
            else:
                print '      subsetting %d seqs with indices %d --> %d' % (args.istartstop[1] - args.istartstop[0], args.istartstop[0], args.istartstop[1])
                if not os.path.exists(os.path.dirname(subsimfname)):
                    os.makedirs(os.path.dirname(subsimfname))
                check_call('head -n1 ' + simfname + ' >' + subsimfname, shell=True)
                check_call('sed -n \'' + str(args.istartstop[0] + 2) + ',' + str(args.istartstop[1] + 1) + ' p\' ' + simfname + '>>' + subsimfname, shell=True)  # NOTE conversion from standdard zero indexing to sed inclusive one-indexing (and +1 for header line)
            simfname = subsimfname
    extras = []

    def output_exists(outfname):
        if os.path.exists(outfname):
            print '                      partition output exists, skipping (%s)' % outfname
            return True
        else:
            return False

    if action == 'cache-data-parameters':
        if os.path.exists(fsdir + '/' + label + '/data'):
            print '                      data parametes exist in %s, skipping ' % (fsdir + '/' + label + '/data')
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
        n_procs = max(1, n_to_partition / 15)  # something like 15 seqs/process to start with
    elif action == 'naive-hamming-partition':
        if os.path.exists(simfname.replace('.csv', '-naive-hamming-partition.csv')):
            print '                      partition output exists, skipping (%s)' % simfname.replace('.csv', '-naive-hamming-partition.csv')
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + simfname.replace('.csv', '-' + action + '.csv')
        extras += ['--n-max-queries', n_to_partition, '--auto-hamming-fraction-bounds']
        n_procs = max(1, n_to_partition / 30)
    elif action == 'vsearch-partition':
        outfname = simfname.replace('.csv', '-' + action + '.csv')
        if output_exists(outfname):
            return
        cmd += ' --simfname ' + simfname
        cmd += ' --outfname ' + outfname
        extras += ['--n-max-queries', n_to_partition, '--naive-vsearch']
        n_procs = max(1, n_to_partition / 100)  # only used for ighutil step
    elif action == 'run-viterbi':
        if os.path.exists(simfname.replace('.csv', '-run-viterbi.csv')):
            print '                      vollmers output exists, skipping (%s)' % simfname.replace('.csv', '-run-viterbi.csv')
            return
        extras += ['--annotation-clustering', 'vollmers', '--annotation-clustering-thresholds', '0.5:0.9']
        cmd += ' --simfname ' + simfname
        extras += ['--n-max-queries', n_to_partition]
        n_procs = max(1, n_to_partition / 50)
    elif action == 'run-changeo':
        changeodir = '/home/dralph/work/changeo/changeo'
        changeo_fsdir = '/fh/fast/matsen_e/dralph/work/changeo'  #.bak'
        _simfbase = leafmutstr(n_leaves, mut_mult).replace('-', '_')
        if os.path.isdir(changeo_fsdir + '/' + _simfbase):
            print '                      already untar\'d into %s' % changeo_fsdir + '/' + _simfbase
        else:
            tar_cmd = 'mkdir ' + changeo_fsdir + '/' + _simfbase + ';'
            tar_cmd += ' tar Jxvf ' + changeo_fsdir + '/' + _simfbase + '.txz --exclude=\'IMGT_HighV-QUEST_individual_files_folder/*\' -C ' + changeo_fsdir + '/' + _simfbase
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
    elif action == 'compare-sample-sizes':
        compare_sample_sizes(label, n_leaves, mut_mult)
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
    logbase = fsdir + '/' + label + '/_logs/' + leafmutstr(n_leaves, mut_mult) + '-' + action
    if args.subset is not None:
        logbase = logbase.replace('_logs/', '_logs/subset-' + str(args.subset) + '/')
    if args.istartstop is not None:
        logbase = logbase.replace('_logs/', '_logs/istartstop-' + '-'.join([str(i) for i in args.istartstop]) + '/')
    if not os.path.exists(os.path.dirname(logbase)):
        os.makedirs(os.path.dirname(logbase))
    proc = Popen(cmd.split(), stdout=open(logbase + '.out', 'w'), stderr=open(logbase + '.err', 'w'))
    procs.append(proc)
    time.sleep(5)

# ----------------------------------------------------------------------------------------
n_to_partition = 1300
if args.istartstop is not None:
    n_to_partition = args.istartstop[1] - args.istartstop[0]
n_data_to_cache = 50000
mutation_multipliers = ['1']  #['1', '4']
n_leaf_list = [10]  #[5, 10, 25, 50]
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
    # label += '.bak'

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
                print '         ----> mutate', mut_mult
                execute(action, label, datafname, n_leaves, mut_mult)
                # sys.exit()

    if 'write-plots' in args.actions:
        write_all_plot_csvs(label)

    break
