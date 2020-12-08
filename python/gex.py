import sys
import time
import os
import scipy
import numpy
import collections
import csv
import operator
import math
import itertools

import utils

mdir = "%s/work/partis/datascripts/meta/goo-dengue-10x" % os.getenv('HOME')
fabio_fname = '%s/fabio-pb-markers.tsv' % mdir
waickfname = '%s/waickman-markers.csv' % mdir
barcodefname = 'barcodes.txt'
pcafname = 'pca.txt'
umapfname = 'umap.txt'
clusterfname = 'clusters.txt'

# ----------------------------------------------------------------------------------------
def markfname(iclust):
    return 'markers-cluster-%d.txt' % iclust  # NOTE R indexing, starts from 1

# ----------------------------------------------------------------------------------------
def install():
    rcmds = ['install.packages("BiocManager", repos="http://cran.rstudio.com/"))',
             'BiocManager::install(c("scRNAseq", "scater", "scran", "uwot", "DropletUtils"), dependencies=TRUE)']  # "TENxPBMCData"
    workdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER'))
    os.makedirs(workdir)
    utils.run_r(rcmds, workdir)
    os.rmdir(workdir)
# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
def loadcmd(lib):
    return 'library(%s, warn.conflicts=F, quietly=T)' % lib

# ----------------------------------------------------------------------------------------
# takes the "dot product" (if normalized, it's cos theta) of two groups of logfc expression values to see how similar they are (yeah this is probably kind of dumb, but it'll give an idea of how similar they are)
def gexdot(gvals1, gvals2=None, normalize=True, recursed=False, lbstr='', debug=False):  # they're both ordered dicts gene : logfc (sorted by decreasing logfc)
    if gvals2 is None:
        gvals2 = gvals1
    dprod = 0.
    common_genes = set(gvals1) & set(gvals2)
    for gene in common_genes:  # loop over genes that they have in common
        dprod += gvals1[gene] * gvals2[gene]
    if normalize:
        dprod /= gexdot(gvals1, normalize=False, recursed=True) + gexdot(gvals2, normalize=False, recursed=True)
    if debug and not recursed:
        if debug > 1:
            lm = max(len(g) for g in gvals1.keys() + gvals2.keys())
            def dstr(vl): return '  '.join('%s %-5.1f'%(utils.color('red' if g in common_genes else None, g, width=lm), v) for g, v in vl.items())
            print '      %s' % dstr(gvals1)
            print '      %s' % dstr(gvals2)
        if len(common_genes) == 0:
            pass  # print '              none in common'
        else:
            print '            %s%5.2f  %2d / (%2d | %2d): %s' % (lbstr, dprod, len(common_genes), len(gvals1), len(gvals2), ' '.join(common_genes))
    return dprod

# ----------------------------------------------------------------------------------------
# read info from previous papers (from fabio + adam waickman, atm)
def read_ref_data():
    fabfo = {vt : [] for vt in ['pb', 'naive']}
    with open(fabio_fname) as ffile:
        reader = csv.DictReader(ffile, delimiter='\t')
        for line in reader:
            pbval, nval = float(line['avg_plasma']), float(line['avg_naive'])
            if pbval == 0 or nval == 0:
                continue
            fabfo['pb'].append((line['GeneName'], math.log(pbval / nval, 2)))  # NOTE these don't really match up with his KS statistic particularly wlel
            fabfo['naive'].append((line['GeneName'], math.log(nval / pbval, 2)))
    # print 'fabio'
    for vtype in fabfo:
        fabfo[vtype] = collections.OrderedDict(sorted(fabfo[vtype], key=operator.itemgetter(1), reverse=True))  # definitely not always sorted
        # print vtype, fabfo[vtype]

    waickfo = {}
    with open(waickfname) as wfile:
        reader = csv.DictReader(wfile)
        for line in reader:
            if line['type'] not in waickfo:
                waickfo[line['type']] = []
            waickfo[line['type']].append((line['gene'], float(line['logfc'])))
    # print 'waick'
    for vtype in waickfo:
        waickfo[vtype] = collections.OrderedDict(sorted(waickfo[vtype], key=operator.itemgetter(1), reverse=True))  # should already be sorted
        # print vtype, waickfo[vtype]

    return fabfo, waickfo

# ----------------------------------------------------------------------------------------
def read_gex(outdir, min_dprod=0.01, debug=True):
    # barcodes
    barcode_vals = []
    with open('%s/%s' % (outdir, barcodefname)) as bfile:
        for il, line in enumerate(bfile):
            lstrs = line.strip().split()
            icount = int(lstrs.pop(0).strip('[]'))
            assert icount == len(barcode_vals) + 1  # <icount> is the R-style (1-based) index of the first element in this line
            barcode_vals += [s.strip('"') for s in lstrs]
    if debug:
        print '    read %d barcodes' % len(barcode_vals)

    # pca values
    rotation_vals = collections.OrderedDict()  # relationship between pca and gene names (map from gene name to list of pca components)
    with open('%s/%s' % (outdir, pcafname)) as pfile:
        pca_comps = None  # names for each pca component (like PC3)
        for il, line in enumerate(pfile):
            if il == 0:
                pca_comps = line.strip().split()
                for ipc, pc in enumerate(pca_comps):
                    assert pc[:2] == 'PC'
                    assert int(pc[2:]) == ipc + 1
                continue
            lstrs = line.strip().split()
            gene = lstrs.pop(0)
            assert len(lstrs) == len(pca_comps)
            rotation_vals[gene] = [float(vstr) for vstr in lstrs]
    if debug:
        print '      %d pca components for %d genes: %s' % (len(pca_comps), len(rotation_vals), ' '.join(rotation_vals))

    # umap values
    umap_vals = []  # list of (x, y) umap values for each cell
    with open('%s/%s' % (outdir, umapfname)) as ufile:
        for il, line in enumerate(ufile):
            lstrs = line.strip().split()
            if il == 0:
                assert lstrs == ['[,%d]'%i for i in [1, 2]]
            else:
                icount = int(lstrs.pop(0).strip('[]').rstrip(','))
                assert icount == len(umap_vals) + 1
                umap_vals.append([float(v) for v in lstrs])
    if debug:
        print '      %d umap values' % len(umap_vals)
    assert len(umap_vals) == len(barcode_vals)

    # cluster assignments
    cluster_vals = []
    with open('%s/%s' % (outdir, clusterfname)) as cfile:
        for il, line in enumerate(cfile):
            lstrs = line.strip().split()
            if lstrs[0] != 'Levels:':
                icount = int(lstrs.pop(0).strip('[]'))
                assert icount == len(cluster_vals) + 1  # <icount> is the R-style (1-based) index of the first element in this line
                cluster_vals += [int(c) for c in lstrs]
            else:  # last line lists the clusters (not sure why they're called "levels"
                cluster_ints = [int(c) for c in lstrs[1:]]  # names of the clusters (1-based integer index)
                assert cluster_ints == list(range(min(cluster_ints), max(cluster_ints) + 1))
                assert set(cluster_ints) == set(cluster_vals)
    if debug:
        print '      %d values in %d clusters: %s' % (len(cluster_vals), len(cluster_ints), ' '.join(str(c) for c in cluster_ints))
    assert len(cluster_vals) == len(barcode_vals)

    # markers for each cluster
    cmarkers = {'%d-%d'%(c1, c2) : [] for c1, c2 in itertools.permutations(cluster_ints, 2)}  # reversing them (1-2 vs 2-1) the values are just the negative of each other, but you don't get all the same genes
    for cname in cluster_ints:
        other_clusters = [c for c in cluster_ints if c != cname]
        n_genes, n_columns = None, None
        with open('%s/%s' % (outdir, markfname(cname))) as cfile:
            for il, line in enumerate(cfile):
                lstrs = line.strip().split()
                if il == 0:  # intro line
                    assert line.find('DataFrame with') == 0
                    assert len(lstrs) == 7
                    n_genes = int(lstrs[2])
                    n_columns = int(lstrs[5])
                elif il == 1:  # column names
                    assert len(lstrs) == n_columns
                    assert lstrs.pop(0) == 'Top'
                    assert lstrs.pop(0) == 'p.value'
                    assert lstrs.pop(0) == 'FDR'
                    assert lstrs.pop(0) == 'summary.logFC'  # as currently configured below, this is the min log fc over the subsequence columns/cluster
                    assert lstrs == ['logFC.%d'%i for i in other_clusters]  # should be a column for each pairwise comparison with another cluster
                elif il == 2:  # column type
                    pass  # eh fuck it
                else:
                    assert len(lstrs) == n_columns + 1  # they don't count the first (name) one
                    gene, top, pval, fdr, summary_logfc = lstrs[:5]
                    logfc_vals = {i : float(l) for i, l in zip(other_clusters, lstrs[5:])}
                    for c2 in logfc_vals:
                        cmarkers['%d-%d'%(cname, c2)].append((gene, logfc_vals[c2]))
    for ckey in cmarkers:
        cmarkers[ckey] = collections.OrderedDict(sorted(cmarkers[ckey], key=operator.itemgetter(1), reverse=True))

    # reference marker genes
    fabfo, waickfo = read_ref_data()

    print '  interpretation: "this cluster looks very <type>-like compared to <clusters>,  based on relative upregulation of <N genes>"'
    print '                                                          fractional'
    print '        type   clusters             N genes               similarity                                genes'
    for cname in cluster_ints:
        print '  %s' % utils.color('green', 'cluster %d' % cname)
        for vtype in waickfo:
            clprods = []
            for ic2, c2 in enumerate([c for c in cluster_ints if c != cname]):
                dprod = gexdot(waickfo[vtype], cmarkers['%d-%d'%(cname, c2)], lbstr='%8s %s '%((vtype+':') if ic2==0 else '', utils.color('blue', str(c2)))) #, debug=True)
                if dprod < min_dprod:
                    continue
                common_genes = set(waickfo[vtype]) & set(cmarkers['%d-%d'%(cname, c2)])
                clprods.append((c2, dprod, common_genes))
            clprods = sorted(clprods, key=operator.itemgetter(1), reverse=True)
            if debug and len(clprods) > 0:
                print '    %s  %-20s  %-20s  %-40s  %s' % (utils.color('purple', vtype, width=8),
                                                        utils.color('blue', ' '.join('%d'%c for c, _, _ in clprods), width=20, padside='right'),
                                                        ' '.join('%d'%len(gl) for _, _, gl in clprods),
                                                        ' '.join('%.2f'%d for _, d, _ in clprods),
                                                        ' '.join(set(g for _, _, gl in clprods for g in gl)),
                )

    # return barcode_vals, rotation_vals, umap_vals, cluster_vals

# ----------------------------------------------------------------------------------------
def run_gex(feature_matrix_fname, outdir, make_plots=True, max_pca_components=25, n_top_genes=50):
    rcmds = [loadcmd(l) for l in ['DropletUtils', 'scater', 'scran', 'pheatmap']]
    rcmds += [
        'options(width=1000)',
        'sce <- read10xCounts("%s")' % feature_matrix_fname,
        'rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)',
        # quality control
        'is.mito <- grepl("^MT-", rownames(sce))',  # figure out which genes are mitochondrial
        'qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))',
        'filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")',  # identifies + removes outliers (in several qc metrics)
        'sce <- sce[, !filtered$discard]',
        # normalization
        'sce <- logNormCounts(sce)',
        # feature selection
        # 'fabio.pb.genes <- read.csv("%s", sep="\t", header=T)$GeneName' % , fabio_fname # $name  # genes from fabio (200 most discriminatory between plasmablast + naive B cell):
        'waick.genes <- read.csv("%s", header=T)$gene' % waickfname,  # 10 most up'd genes for naive, memory, pb, and prepb (40 total)
        'genelist <- waick.genes',  # fabio.pb.genes
        'print(sprintf("  using %d genes: %s", length(genelist), paste(genelist, collapse=" ")))',
        'gene.bools <- rowData(sce)$Symbol %in% genelist',  # $ID
        # dimensionality reduction
        'set.seed(1)',
        'n.comp <- min(%d, as.integer(length(genelist)/2))' % max_pca_components,
        'print(sprintf("running pca with %d components", n.comp))',
        'sce <- runPCA(sce, ncomponents=n.comp, subset_row=gene.bools)',
        'sce <- runUMAP(sce, dimred="PCA", external_neighbors=TRUE)',  # uses pca results from previous step TODO test variety of N neighbors and min_dist values
        # clustering
        'g <- buildSNNGraph(sce, use.dimred="PCA")',
        'colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)',
        # write output files (more written below)
        'capture.output(colData(sce)$Barcode, file="%s/%s")' % (outdir, barcodefname),
        'capture.output(attr(reducedDims(sce)$PCA, "rotation"), file="%s/%s")' % (outdir, pcafname),  # pca to gene name rotation
        # (reducedDim(sce, "PCA")[,]  # a table of the pca values for each cell
        'capture.output(reducedDim(sce, "UMAP")[,], file="%s/%s")' % (outdir, umapfname),
        'capture.output(colLabels(sce), file="%s/%s")' % (outdir, clusterfname),
    ]
    if make_plots:
        rcmds += [
            ## pdf(sprintf("%s/clusters.pdf", outdir))
            'png("%s/clusters.png")' % outdir,
            'plotUMAP(sce, colour_by="label")',
            'dev.off()',
        ]
    # find marker genes
    rcmds += [
        'markers <- findMarkers(sce)',  # <markers>: list of data frames for each cluster NOTE this uses *all* the genes, and i can't figure out a way to tell it not to
        'n.genes <- %d' % n_top_genes,
        'print(sprintf("  top %d genes for each cluster (total size %d)", n.genes, length(sce$label)))',
        'for(ich in seq(length(markers))) {'  # look at genes that distinguish cluster ich from all other clusters
        '    print(sprintf("   cluster %2d  size %4d  frac %.2f", ich, sum(sce$label==ich), sum(sce$label==ich) / length(sce$label)))',
        '    interesting <- markers[[ich]]',
        '    capture.output(interesting[1:n.genes,], file=sprintf("%s/markers-cluster-%%d.txt", ich))' % outdir,
        '    best.set <- interesting[interesting$Top <= n.genes,]',  # look at the top N genes from each pairwise comparison
        '    logFCs <- getMarkerEffects(best.set)',
    ]
    if make_plots:
        rcmds += [
            '    png(sprintf("%s/heatmap-%%d.png", ich))' % outdir,
            '    pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))',
            '    dev.off()',
        ]
    rcmds += [
        '}',
    ]

    workdir = utils.choose_random_subdir('/tmp/%s' % os.getenv('USER'))
    os.makedirs(workdir)
    utils.run_r(rcmds, workdir)
    os.rmdir(workdir)
