#!/usr/bin/env python
import csv
import time
import sys
from subprocess import check_output, check_call, Popen
sys.path.insert(1, './python')
import random
from collections import OrderedDict
import math
import os

import utils
import seqfileopener

n_set_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 20, 25, 50]
n_set_max = n_set_list[-1]
baseplotdir = os.getenv('www') + '/partis/n-sets'
parameter_dir = '/fh/fast/matsen_e/dralph/work/partis-dev/_output/021-018/simu-2.3-leaves-1.0-mutate-zipf'
outputdir = parameter_dir + '/n-sets'
simfname = outputdir + '/' + str(n_set_max) + '-leaves.csv'

# ----------------------------------------------------------------------------------------
def get_subfname(n_set):
    return outputdir + '/' + str(n_set) + '-sub-leaves.csv'

# ----------------------------------------------------------------------------------------
def divide_simulation():
    allinfo = {}
    headers = None
    with open(simfname) as csvfile:
        reader = csv.DictReader(csvfile)
        for line in reader:
            if headers is None:
                headers = line.keys()
            if line['reco_id'] not in allinfo:
                allinfo[line['reco_id']] = []
            allinfo[line['reco_id']].append(line)
    
    for n_to_keep in n_set_list:
        outfname = get_subfname(n_to_keep)
        outfile = open(outfname, 'w')
        writer = csv.DictWriter(outfile, headers)
        writer.writeheader()
        for recoid in allinfo:
            already_used = set()
            for _ in range(n_to_keep):
                iseq = None
                while iseq is None or iseq in already_used:
                    iseq = random.randint(0, len(allinfo[recoid]) - 1)
                already_used.add(iseq)
                writer.writerow(allinfo[recoid][iseq])

# ----------------------------------------------------------------------------------------
def run_inference(algorithm):
    base_cmd = './bin/partis run-' + algorithm + ' --parameter-dir ' + parameter_dir + '/hmm --plot-performance --is-simu --slurm'
    for n_set in n_set_list:
        # if n_set == 5:
        #     continue
        n_procs = max(5, int(20. * (n_set / 25.)))
        cmd = base_cmd + ' --n-sets ' + str(n_set) + ' --seqfile ' + get_subfname(n_set)
        cmd += ' --plotdir ' + baseplotdir + '/' + str(n_set) + ' --n-procs ' + str(n_procs)
        cmd += ' --workdir /fh/fast/matsen_e/dralph/work/partis-dev/_tmp/' + str(random.randint(0, 99999))
        cmd += ' --outfname ' + outputdir + '/' + str(n_set) + '-' + algorithm + '.csv'
        print cmd
        Popen(cmd.split())
        time.sleep(0.5)

# ----------------------------------------------------------------------------------------
def peruse_naive_seqs():
    from hist import Hist
    # hall = Hist(n_set_list[-1], n_set_list[0] - 0.5, n_set_list[-1] + 0.5)
    means = []
    for n_set in n_set_list:
        plotdir = baseplotdir + '/' + str(n_set)
        hist = Hist(fname=plotdir + '/hmm/hamming_to_true_naive.csv')
        print '%2d   %.2f' % (n_set, hist.get_mean()),
        # hall.set_ibin(hall.find_bin(n_set), hist.get_mean())
        means.append(hist.get_mean())
    
    import plotting
    fig, ax = plotting.mpl_init()
    # hall.mpl_plot(ax)
    ax.plot(n_set_list, means, marker='.')
    plotting.mpl_finish(ax, baseplotdir, 'means', xlabel='N simultaneous seqs', ylabel='mean hamming to true naive', ybounds=(0, None))

# ----------------------------------------------------------------------------------------
def get_deviations(lps, i_baseline, signed=False):
    deviations = []
    for n_set in n_set_list:
        vals = []
        # print n_set
        for reco_id in lps[n_set]:
            val = lps[n_set][reco_id]
            ref_val = lps[n_set_list[i_baseline]][reco_id]
            deviation = 0.
            if ref_val != 0.:
                deviation = (val - ref_val) / ref_val
            if not signed:
                deviation = abs(deviation)
            # print '  %10.4f  %10.4f  %8.4f' % (val, ref_val, deviation)
            vals.append(deviation)
        deviations.append(sum(vals) / float(len(vals)))
    return deviations

# ----------------------------------------------------------------------------------------
def func(n, a, b):
    # return n / (1. - a / pow(n, b))
    return a * n  + b

# ----------------------------------------------------------------------------------------
def peruse_forward_scores():
    _, reco_info = seqfileopener.get_seqfile_info(simfname, is_data=False)  #, n_max_queries=10000)
    logprobs, partialcorr_logprobs, corr_logprobs = OrderedDict(), OrderedDict(), OrderedDict()
    for n_set in n_set_list:
        print n_set
        # if n_set != 5:
        #     continue
        logprobs[n_set], partialcorr_logprobs[n_set], corr_logprobs[n_set] = OrderedDict(), OrderedDict(), OrderedDict()
        with open(outputdir + '/' + str(n_set) + '-forward.csv') as csvfile:
            reader = csv.DictReader(csvfile)
            for line in reader:
                uidlist = line['unique_ids'].split(':')
                assert utils.from_same_event(reco_info, uidlist)
                reco_id = reco_info[uidlist[0]]['reco_id']
                if reco_id in logprobs[n_set]:
                    raise Exception('already had %s' % reco_id)

                logprobs[n_set][reco_id] = float(line['logprob'])

                factor = 1. / n_set
                partialcorr_logprobs[n_set][reco_id] = factor * float(line['logprob'])

                factor = (1. - 0.24 / pow(float(n_set), 0.9)) / n_set
                # factor = 1. / (0.77547824*n_set + 0.20327936)
                corr_logprobs[n_set][reco_id] = factor * float(line['logprob'])


    i_baseline = -1
    deviations = get_deviations(logprobs, i_baseline)
    # fit_stuff(n_set_list, deviations)
    partialcorr_deviations = get_deviations(partialcorr_logprobs, i_baseline)
    signed_partialcorr_deviations = get_deviations(partialcorr_logprobs, i_baseline, signed=True)
    corr_deviations = get_deviations(corr_logprobs, i_baseline)
    signed_corr_deviations = get_deviations(corr_logprobs, i_baseline, signed=True)

    import plotting
    fig, ax = plotting.mpl_init()
    ax.plot(n_set_list, deviations, marker='.')
    plotting.mpl_finish(ax, baseplotdir, 'forwards', xlabel='N simultaneous seqs', ylabel='log prob deviation to ' + str(n_set_list[i_baseline]))  #, ybounds=(-0.02, 0.02))

    # fig, ax = plotting.mpl_init()
    # ax.plot(n_set_list, partialcorr_deviations, marker='.')
    # ax.plot([n_set_list[0], n_set_list[-1]], [0, 0])
    # plotting.mpl_finish(ax, baseplotdir, 'partially-corrected-forwards', xlabel='N simultaneous seqs', ylabel='log prob deviation to ' + str(n_set_list[i_baseline])) #, ybounds=(-0.02, 0.02))

    fig, ax = plotting.mpl_init()
    ax.plot(n_set_list, partialcorr_deviations, marker='.', label='1/n (abs)')
    ax.plot(n_set_list, signed_partialcorr_deviations, marker='.', label='1/n')
    ax.plot(n_set_list, corr_deviations, marker='.', label='1/crap (abs)')
    ax.plot(n_set_list, signed_corr_deviations, marker='.', label='1/crap')
    ax.plot([n_set_list[0], n_set_list[-1]], [0, 0])
    plotting.mpl_finish(ax, baseplotdir, 'corrected-forwards', xlabel='N simultaneous seqs', ylabel='log prob deviation to ' + str(n_set_list[i_baseline])) #, ybounds=(-0.02, 0.02))

    fig, ax = plotting.mpl_init()
    ax.plot(n_set_list, signed_corr_deviations, marker='.')
    ax.plot([n_set_list[0], n_set_list[-1]], [0, 0])
    plotting.mpl_finish(ax, baseplotdir, 'signed-corrected-forwards', xlabel='N simultaneous seqs', ylabel='log prob deviation to ' + str(n_set_list[i_baseline])) #, ybounds=(-0.02, 0.02))

# ----------------------------------------------------------------------------------------
def fit_stuff(xvals, yvals):
    print xvals
    print yvals
    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit
    import numpy
    
    # coeffs, matcov = curve_fit(func, xvals, yvals, p0=[1., -3.])
    # print coeffs
    # print matcov

    # yaj = func(xvals, coeffs[0], coeffs[1])
    # plt.plot(xvals, yvals, 'x', xvals, yaj, 'r-')
    # plt.savefig(os.getenv('www') + '/partis/tmp/foo.png')

    # calculate polynomial
    z = numpy.polyfit(xvals, yvals, 1)
    print z
    f = numpy.poly1d(z)
    
    # calculate new x's and y's
    x_new = numpy.linspace(xvals[0], xvals[-1], 50)
    y_new = f(x_new)
    
    # plt.plot(xvals, yvals, 'o', x_new, y_new)
    # plt.xlim([xvals[0]-1, xvals[-1] + 1 ])

    residuals = [(f(xvals[i]) - yvals[i]) / yvals[i] for i in range(len(xvals))]
    # for i in range(len(xvals)):
    #     print xvals[i], yvals[i], f(xvals[i]), residuals[i]
    plt.plot(xvals, residuals, 'o')

    plt.savefig(os.getenv('www') + '/partis/tmp/foo.png')
    
# ----------------------------------------------------------------------------------------
def simulate():
    cmd = './bin/partis simulate --constant-number-of-leaves --n-sim-events 2000 --n-leaves ' + str(n_set_max) + ' --outfname ' + simfname + ' --parameter-dir ' + parameter_dir + '/hmm --n-procs 10'
    Popen(cmd.split())

# simulate()
# divide_simulation()
# run_inference('viterbi')
peruse_naive_seqs()
# run_inference('forward')
# peruse_forward_scores()
# xvals = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 20, 25, 50]
# yvals = [0.0, 0.7653082061900935, 1.5329842153963764, 2.311362454137212, 3.121540811338117, 3.90620744902059, 4.6167805465227065, 5.365305056511266, 6.232965310537493, 6.915898187309928, 7.781688535833085, 10.749505939689282, 14.569377524709266, 18.64614056231728, 38.01829568891039]
# fit_stuff(xvals, yvals)
