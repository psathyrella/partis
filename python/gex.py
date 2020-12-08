import sys
import time
import os
import scipy
import numpy
import collections
import csv
import operator

import utils

mdir = "%s/work/partis/datascripts/meta/goo-dengue-10x" % os.getenv('HOME')
fabio_fname = '%s/plasmablast_markers.tsv' % mdir
waickfname = '%s/waickman-markers.csv' % mdir
barcodefname = 'barcodes.txt'
pcafname = 'pca.txt'
umapfname = 'umap.txt'
clusterfname = 'clusters.txt'

# ----------------------------------------------------------------------------------------
def gmarkfname(iclust):
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
# read info from previous papers (from adam waickman + fabio, atm)
def read_ref_data():
    waickfo = {}
    with open(waickfname) as wfile:
        reader = csv.DictReader(wfile)
        for line in reader:
            if line['type'] not in waickfo:
                waickfo[line['type']] = []
            waickfo[line['type']].append((line['gene'], float(line['logfc'])))
    for vtype in waickfo:
        waickfo[vtype] = collections.OrderedDict(sorted(waickfo[vtype], key=operator.itemgetter(1), reverse=True))  # should already be sorted

# ----------------------------------------------------------------------------------------
def read_gex(outdir, debug=True):
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
                cluster_ints = set(int(c) for c in lstrs[1:])
                assert cluster_ints == set(cluster_vals)
    if debug:
        print '      %d values in %d clusters: %s' % (len(cluster_vals), len(cluster_ints), ' '.join(str(c) for c in cluster_ints))
    assert len(cluster_vals) == len(barcode_vals)

    read_ref_data()

    # return barcode_vals, rotation_vals, umap_vals, cluster_vals

# ----------------------------------------------------------------------------------------
def run_gex(feature_matrix_fname, outdir, make_plots=True, max_pca_components=25, n_top_genes=10):
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
