#!/usr/bin/env python
# Builds the per-locus SHM parameter dirs in this directory (<locus>/hmm/) by pooling the per-gene, per-position
# mutation counts that partis cache-parameters already wrote for a set of public, mutated human samples. See README.md.
#   ./make-scratch-parameters.py --base-dir <dir containing the datasets in datasets.csv>
import argparse, csv, glob, os, random, sys, collections
partis_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(1, partis_dir)
from partis import utils, glutils, fraction_uncertainty
from partis.hist import Hist

this_dir = os.path.dirname(os.path.abspath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('--base-dir', required=True, help='dir under which the globs in --dataset-fname are found')
parser.add_argument('--dataset-fname', default=this_dir + '/datasets.csv')
parser.add_argument('--outdir', default=this_dir)
parser.add_argument('--loci', default='igh:igk:igl')
parser.add_argument('--max-seqs', type=int, default=1000000, help='stop adding samples once a locus has this many (weighted) sequences')
parser.add_argument('--max-seqs-per-sample', type=int, default=10000, help='scale down the counts of any sample larger than this, so it contributes as if it had this many sequences (so many individuals contribute, rather than a few very deep ones)')
parser.add_argument('--min-mean-mute-freq', type=float, default=0.02, help='skip samples whose mean mutation frequency is below this (mostly naive)')
parser.add_argument('--random-seed', type=int, default=1)
parser.add_argument('--germline-dir', default=partis_dir + '/data/germlines/human')
args = parser.parse_args()
args.loci = utils.get_arg_list(args.loci)
hist_fnames = ['%s-mean-%s.csv' % (r, t) for r in ['all'] + utils.regions for t in ['mute-freqs', 'n-muted']]

# ----------------------------------------------------------------------------------------
def choose_samples(locus):  # go round robin over datasets, adding one random sample from each in turn, until we have --max-seqs sequences
    dsets = collections.OrderedDict()
    for line in csv.DictReader(open(args.dataset_fname)):
        if line['locus'] != locus:
            continue
        for hdir in sorted(glob.glob('%s/%s' % (args.base_dir, line['glob']))):
            if line['exclude'] and line['exclude'] in hdir:
                continue
            if not all(os.path.exists('%s/%s' % (hdir, f)) for f in ['v_gene-probs.csv', 'all-mean-mute-freqs.csv', 'mute-freqs']):
                continue
            n_seqs = sum(int(l['count']) for l in csv.DictReader(open(hdir + '/v_gene-probs.csv')))
            mean_mfreq = Hist(fname=hdir + '/all-mean-mute-freqs.csv').get_mean()
            if mean_mfreq < args.min_mean_mute_freq:
                continue
            dsets.setdefault(line['dataset'], []).append({'hdir' : hdir, 'n_seqs' : n_seqs, 'mean_mfreq' : mean_mfreq, 'source' : line['source']})
    rng = random.Random(args.random_seed)
    for dname in dsets:
        rng.shuffle(dsets[dname])
    chosen, n_tot = [], 0
    while n_tot < args.max_seqs and any(len(s) > 0 for s in dsets.values()):
        for dname in dsets:
            if len(dsets[dname]) == 0 or n_tot >= args.max_seqs:
                continue
            sfo = dsets[dname].pop()
            sfo['dataset'] = dname
            sfo['weight'] = min(1., args.max_seqs_per_sample / float(sfo['n_seqs']))
            chosen.append(sfo)
            n_tot += sfo['weight'] * sfo['n_seqs']
    return chosen

# ----------------------------------------------------------------------------------------
def pool_locus(locus, samples):
    default_glfo = glutils.read_glfo(args.germline_dir, locus)
    counts = {}  # counts[gene][position][nuke]
    gene_counts = {r : collections.Counter() for r in utils.regions}
    hists = {}
    n_skipped = collections.Counter()
    for sfo in samples:
        glfo = glutils.read_glfo(sfo['hdir'] + '/germline-sets', locus)
        ok_genes = set(g for r in utils.regions for g, seq in glfo['seqs'][r].items() if default_glfo['seqs'][r].get(g) == seq)  # default alleles with identical sequence, i.e. no inferred alleles
        for region in utils.regions:
            fn = '%s/%s_gene-probs.csv' % (sfo['hdir'], region)
            if not os.path.exists(fn):
                continue
            for line in csv.DictReader(open(fn)):
                if line[region + '_gene'] in ok_genes:
                    gene_counts[region][line[region + '_gene']] += sfo['weight'] * int(line['count'])
                else:
                    n_skipped['gene-probs'] += 1
        for fn in glob.glob(sfo['hdir'] + '/mute-freqs/*.csv'):
            gene = utils.unsanitize_name(os.path.basename(fn).replace('.csv', ''))
            if gene not in ok_genes:
                n_skipped['mute-freqs'] += 1
                continue
            gcounts = counts.setdefault(gene, {})
            for line in csv.DictReader(open(fn)):
                pcounts = gcounts.setdefault(int(line['position']), collections.Counter())
                for nuke in utils.nukes:
                    pcounts[nuke] += sfo['weight'] * int(line[nuke + '_obs'])
        for hfn in hist_fnames:
            if not os.path.exists('%s/%s' % (sfo['hdir'], hfn)):
                continue
            hist = Hist(fname='%s/%s' % (sfo['hdir'], hfn))
            hist.bin_contents = [c * sfo['weight'] * sfo['n_seqs'] for c in hist.bin_contents]  # they're written normalized, so weight each sample by its (capped) number of sequences
            if hfn not in hists:
                hists[hfn] = hist
            elif hists[hfn].low_edges != hist.low_edges:
                n_skipped['hist with different binning'] += 1
            else:
                hists[hfn].add(hist)

    odir = '%s/%s/hmm' % (args.outdir, locus)
    utils.prep_dir(odir + '/mute-freqs', wildlings=['*.csv'])
    for gene, gcounts in counts.items():
        gl_seq = default_glfo['seqs'][utils.get_region(gene)][gene]
        with open('%s/mute-freqs/%s.csv' % (odir, utils.sanitize_name(gene)), utils.csv_wmode()) as ofile:  # mutefreqer.write_single_gene() format, but only the columns partis reads (paramutils.read_mute_counts() and read_mute_freqs_with_weights()), to keep the files small: the per-base frequencies and their errors can be recomputed from the counts
            writer = csv.DictWriter(ofile, ['position', 'mute_freq', 'lo_err', 'hi_err'] + [n + '_obs' for n in utils.nukes])
            writer.writeheader()
            for pos in sorted(gcounts):
                gcounts[pos] = {n : int(round(c)) for n, c in gcounts[pos].items()}  # weighted sums, so round back to counts
                total = sum(gcounts[pos].values())
                if total == 0:
                    continue
                n_mutated = total - gcounts[pos].get(gl_seq[pos], 0)
                row = {'position' : pos, 'mute_freq' : n_mutated / float(total)}
                row['lo_err'], row['hi_err'] = fraction_uncertainty.err(n_mutated, total)[:2]
                row['mute_freq'], row['lo_err'], row['hi_err'] = ['%.4g' % v for v in (row['mute_freq'], row['lo_err'], row['hi_err'])]
                for nuke in utils.nukes:
                    row[nuke + '_obs'] = gcounts[pos].get(nuke, 0)
                writer.writerow(row)
    for region in utils.regions:
        with open('%s/%s_gene-probs.csv' % (odir, region), utils.csv_wmode()) as ofile:
            writer = csv.DictWriter(ofile, [region + '_gene', 'count'])
            writer.writeheader()
            for gene, count in gene_counts[region].most_common():
                writer.writerow({region + '_gene' : gene, 'count' : int(round(count))})
    for hfn, hist in hists.items():
        hist.normalize()
        hist.write('%s/%s' % (odir, hfn))
    glutils.write_glfo(odir + '/' + glutils.glfo_dir, default_glfo)
    with open('%s/%s/samples.csv' % (args.outdir, locus), utils.csv_wmode()) as sfile:
        writer = csv.DictWriter(sfile, ['dataset', 'sample', 'n_seqs', 'weight', 'mean_mute_freq', 'source'])
        writer.writeheader()
        for sfo in samples:
            writer.writerow({'dataset' : sfo['dataset'], 'sample' : os.path.relpath(sfo['hdir'], args.base_dir), 'n_seqs' : sfo['n_seqs'], 'weight' : '%.4g' % sfo['weight'], 'mean_mute_freq' : '%.4f' % sfo['mean_mfreq'], 'source' : sfo['source']})
    print('  %s: pooled %d samples (%d weighted seqs) from %s' % (locus, len(samples), sum(s['weight'] * s['n_seqs'] for s in samples), ', '.join('%s %d' % (d, n) for d, n in collections.Counter(s['dataset'] for s in samples).items())))
    print('       %d genes with mute freqs, skipped (not default alleles or different binning): %s' % (len(counts), dict(n_skipped)))

# ----------------------------------------------------------------------------------------
for locus in args.loci:
    samples = choose_samples(locus)
    pool_locus(locus, samples)
