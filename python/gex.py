from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
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

from . import utils
from io import open

mdir = "%s/work/partis/datascripts/meta/goo-dengue-10x/gex-markers" % os.getenv('HOME')
fabio_fname = '%s/fabio-pb-markers.tsv' % mdir
waickfname = '%s/waickman-markers.csv' % mdir
msigdbfname = '%s/msigdb-markers.csv' % mdir
barcodefname = 'barcodes.txt'
pcafname = 'pca.txt'
umapfname = 'umap.txt'
clusterfname = 'clusters.txt'
cell_type_fname = 'cell-types.csv'
cluster_vs_subtype_fname = 'clusters-vs-subtype.csv'

msdsets = [  # still need _UP or _DN tacked on at the end
    ('gc', 'GSE4142_GC_BCELL_VS_MEMORY_BCELL'),  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4517294/ and https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1413911/
    ('memory', 'GSE42724_MEMORY_BCELL_VS_PLASMABLAST'),  # https://pubmed.ncbi.nlm.nih.gov/23613519/
    ('naive', 'GSE4142_NAIVE_BCELL_VS_PLASMA_CELL'),
    ('naive', 'GSE4142_NAIVE_VS_GC_BCELL'),
    ('naive', 'GSE4142_NAIVE_VS_MEMORY_BCELL'),
    ('naive', 'GSE42724_NAIVE_BCELL_VS_PLASMABLAST'),
    ('naive', 'GSE42724_NAIVE_VS_MEMORY_BCELL'),
    ('plasma', 'GSE4142_PLASMA_CELL_VS_GC_BCELL'),
    ('plasma', 'GSE4142_PLASMA_CELL_VS_MEMORY_BCELL'),
]
def msigdb_sets(updown):
    assert updown in ['UP', 'DN']
    return [(c, '%s_%s'%(n, updown)) for c, n in msdsets]

# ----------------------------------------------------------------------------------------
def markfname(iclust):
    return 'markers-cluster-%d.csv' % iclust  # NOTE R indexing, starts from 1

# ----------------------------------------------------------------------------------------
def install():
    rcmds = ['install.packages("BiocManager", repos="http://cran.rstudio.com/")',
             'BiocManager::install(c("scRNAseq", "scater", "scran", "uwot", "DropletUtils", "GSEABase", "AUCell", "celldex", "SingleR"), dependencies=TRUE)']  # "TENxPBMCData"
    # maybe should add these, since the regular package installer tried to update them but couldn't?
    # 'beachmat', 'BiocNeighbors', 'BiocStyle', 'biomaRt', 'DelayedArray', 'DelayedMatrixStats', 'edgeR', 'gdsfmt', 'GenomeInfoDb', 'HDF5Array', 'IRanges', 'MatrixGenerics', 'preprocessCore', 'Rhdf5lib', 'S4Vectors', 'scuttle', 'sparseMatrixStats'
    utils.run_r(rcmds, 'auto')
# install()
# sys.exit()

# ----------------------------------------------------------------------------------------
def loadcmd(lib):
    return 'library(%s, warn.conflicts=F, quietly=T)' % lib

# ----------------------------------------------------------------------------------------
def rplotcmds(plotdir, plotname, pcmd, rowcol=None, hw=None, ftype='png'):
    rcmds = [
        '%s("%s/%s.%s")' % (ftype, plotdir, plotname, ftype),
        pcmd,
        'dev.off()',
    ]
    if rowcol is not None:  # pair of (row, column) values for layout() command
        rcmds.insert(1, 'layout(mat=matrix(c(%s), nrow=%d, ncol=%d, byrow=T))' % (', '.join(str(i) for i in range(1, rowcol[0]*rowcol[1] + 1)), rowcol[0], rowcol[1]))
    if hw is not None:  # pair of (width, height)
        rcmds[0] = rcmds[0].rstrip(')') + ', width=%d, height=%d)' % tuple(hw)
    return rcmds

# ----------------------------------------------------------------------------------------
def dimredcmds(outdir, glist_name, max_pca_components=25, n_top_genes=100):
    # feature selection
    rcmds = [
        'print(sprintf("  using %%d genes: %%s", length(%s), paste(%s, collapse=" ")))' % (glist_name, glist_name),
        'gene.bools <- rowData(sce)$Symbol %%in%% %s' % glist_name,  # $ID
        # dimensionality reduction
        'set.seed(1)',
        'n.comp <- min(%d, as.integer(length(%s)/2))' % (max_pca_components, glist_name),
        'print(sprintf("running pca with %d components", n.comp))',
        'sce <- runPCA(sce, ncomponents=n.comp, subset_row=gene.bools)',
        'sce <- runUMAP(sce, dimred="PCA", external_neighbors=TRUE)',  # uses pca results from previous step TODO test variety of N neighbors and min_dist values
        # clustering
        'g <- buildSNNGraph(sce, use.dimred="PCA")',  # guild graph
        'colLabels(sce) <- factor(igraph::cluster_louvain(g)$membership)',  # use graph to cluster, and add the resulting labels to <sce>
        'capture.output(attr(reducedDims(sce)$PCA, "rotation"), file="%s/%s")' % (outdir, pcafname),  # pca to gene name rotation
        # (reducedDim(sce, "PCA")[,]  # a table of the pca values for each cell
        'capture.output(reducedDim(sce, "UMAP")[,], file="%s/%s")' % (outdir, umapfname),  # umap pair for each cell
        'capture.output(colLabels(sce), file="%s/%s")' % (outdir, clusterfname),  # cluster label for each cell
    ]
    rcmds += rplotcmds(outdir, 'clusters', 'plotUMAP(sce, colour_by="label")')

    # find marker genes
    rcmds += [
        'markers <- findMarkers(sce)',  # <markers>: list of data frames for each cluster NOTE this uses *all* the genes, and i can't figure out a way to tell it not to
        'print(sprintf("  top %d genes for each cluster (total size %%d)", length(sce$label)))' % n_top_genes,
        'for(ich in seq(length(markers))) {'  # look at genes that distinguish cluster ich from all other clusters
        '    print(sprintf("   cluster %2d  size %4d  frac %.2f", ich, sum(sce$label==ich), sum(sce$label==ich) / length(sce$label)))',
        '    interesting <- markers[[ich]]',
        '    best.set <- interesting[interesting$Top <= %d,]' % n_top_genes,  # takes all genes that were in the top N for any pairwise comparison
        '    write.csv(best.set, sprintf("%s/markers-cluster-%%d.csv", ich))' % outdir,
        '    logFCs <- getMarkerEffects(best.set)',
    ]
    # if make_plots:  # these plots aren't really readable any more with n_top_genes more than 10 or so
    #     rcmds += [
    #         # rplotcmds(outdir, 'sprintf("%s/heatmap-%%d", ich)',  # arg, this won't work this way
    #         '    png(sprintf("%s/heatmap-%%d.png", ich))' % outdir,
    #         '    pheatmap(logFCs, breaks=seq(-5, 5, length.out=101))',
    #         '    dev.off()',
    #     ]
    rcmds += [
        '}',
    ]

    return rcmds

# ----------------------------------------------------------------------------------------
def run_msigdbr(outdir):  # download the sets and write to csvs
    # NOTE still had to sort|uniq|sort -t, -k2 this by hand (after removing first column with just line numbers)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    rcmds = [
        loadcmd('msigdbr'),
        'all_gene_sets <- msigdbr(species="Homo sapiens", category="C7")',
        'alldf <- data.frame()',
    ]
    for ctype, gsname in msigdb_sets('UP'):  # TODO should probably use the 'DN' ones in some way?
        print('    %8s %s' % (ctype, gsname))
        rcmds += [
            'glist <- all_gene_sets[all_gene_sets$gs_name=="%s",]$human_gene_symbol' % gsname,  # gives list of gene names
            'df <- data.frame(glist, tag="%s")' % ctype,
            'names(df)[names(df) == "glist"] <- "gene"',
            'names(df)[names(df) == "tag"] <- "type"',
            'alldf <- rbind(alldf, df)'
        ]
    rcmds += [
        'write.csv(alldf, "%s/msigdb-markers.csv")' % outdir,
    ]
    utils.run_r(rcmds, 'auto', dryrun=False)

# ----------------------------------------------------------------------------------------
def ctype_ann_cmds(outdir, clnames):  # cell type annotation (although we do some with celldex at the start as well)
    # using custom references
    rcmds = [
        'clabels <- %s[%s$type!="", ]' % (clnames, clnames),  # remove rows/genes with empty 'type'
        'ctypes <- clabels$type[!duplicated(clabels$type)]',  # just gets all values for the 'type' column
        'all.sets <- lapply(ctypes, function(x) { GeneSet(lapply(clabels, `[`, clabels$type==x)$gene, setName=x) })',
        'all.sets <- GeneSetCollection(all.sets)',
        'rankings <- AUCell_buildRankings(counts(sce), plotStats=FALSE, verbose=FALSE)',
        'cell.aucs <- AUCell_calcAUC(all.sets, rankings)',
        'results <- t(assay(cell.aucs))',
        'new.labels <- colnames(results)[max.col(results)]',
        'write.csv(cbind(barcode=colData(sce)$Barcode, results, new.labels), "%s/%s")' % (outdir, cell_type_fname),
        'tab <- table(new.labels, sce$label)',  # only if we have clusters
        'write.csv(tab, "%s/%s")' % (outdir, cluster_vs_subtype_fname),
    ]
    rcmds += rplotcmds(outdir, 'auc-thresholds', 'AUCell_exploreThresholds(cell.aucs, plotHist=TRUE, assign=TRUE)', rowcol=(2, 2), hw=(1500, 1500))  # this is verbose as all hell

    return rcmds

# ----------------------------------------------------------------------------------------
# takes the "dot product" (if normalized, it's cos theta) of two groups of logfc expression values to see how similar they are (yeah this is probably kind of dumb, but it'll give an idea of how similar they are)
def gexdot(gvals1, gvals2=None, normalize=True, recursed=False, return_gene_contributions=False, lbstr='', debug=False):  # they're both ordered dicts gene : logfc (sorted by decreasing logfc)
    if gvals2 is None:
        gvals2 = gvals1
    dprod = 0.
    common_genes = set(gvals1) & set(gvals2)
    gene_contribs = {}
    for gene in common_genes:  # loop over genes that they have in common
        gene_contribs[gene] = gvals1[gene] * gvals2[gene]
        dprod += gene_contribs[gene]
    if normalize:
        dprod /= gexdot(gvals1, normalize=False, recursed=True) + gexdot(gvals2, normalize=False, recursed=True)
    if debug and not recursed:
        if debug > 1:
            lm = max(len(g) for g in list(gvals1.keys()) + list(gvals2.keys()))
            def dstr(vl): return '  '.join('%s %-5.1f'%(utils.color('red' if g in common_genes else None, g, width=lm), v) for g, v in vl.items())
            print('      %s' % dstr(gvals1))
            print('      %s' % dstr(gvals2))
        if len(common_genes) == 0:
            pass  # print '              none in common'
        else:
            print('            %s%5.2f  %2d / (%2d | %2d): %s' % (lbstr, dprod, len(common_genes), len(gvals1), len(gvals2), ' '.join(common_genes)))
    if return_gene_contributions:
        return dprod, gene_contribs
    else:
        return dprod

# ----------------------------------------------------------------------------------------
# read info from previous papers (from fabio + adam waickman, atm)
def read_ref_data():
    fabfo = {vt : [] for vt in ['pb', 'naive']}
    with open(fabio_fname) as ffile:
        reader = csv.DictReader(ffile, delimiter=str('\t'))
        for line in reader:
            pbval, nval = float(line['avg_plasma']), float(line['avg_naive'])
            if pbval == 0 or nval == 0:
                continue
            fabfo['pb'].append((line['gene'], math.log(pbval / nval, 2)))  # NOTE these don't really match up with his KS statistic particularly well
            fabfo['naive'].append((line['gene'], math.log(nval / pbval, 2)))
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
def read_gex(outdir, min_dprod=0.001, debug=True):
    # barcodes
    barcode_vals = []
    with open('%s/%s' % (outdir, barcodefname)) as bfile:
        for il, line in enumerate(bfile):
            lstrs = line.strip().split()
            icount = int(lstrs.pop(0).strip('[]'))
            assert icount == len(barcode_vals) + 1  # <icount> is the R-style (1-based) index of the first element in this line
            barcode_vals += [s.strip('"') for s in lstrs]
    if debug:
        print('    read %d barcodes' % len(barcode_vals))

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
        print('      %d pca components for %d genes: %s' % (len(pca_comps), len(rotation_vals), ' '.join(rotation_vals)))

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
        print('      %d umap values' % len(umap_vals))
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
        print('      %d values in %d clusters: %s' % (len(cluster_vals), len(cluster_ints), ' '.join(str(c) for c in cluster_ints)))
    assert len(cluster_vals) == len(barcode_vals)

    # markers for each cluster
    pairwise_cmarkers = {'%d-%d'%(c1, c2) : [] for c1, c2 in itertools.permutations(cluster_ints, 2)}  # reversing them (1-2 vs 2-1) the values are just the negative of each other if they're both there, but you don't get all the same genes
    summary_cmarkers = {'%d-summary'%c : [] for c in cluster_ints}
    for cname in cluster_ints:
        other_clusters = [c for c in cluster_ints if c != cname]
        with open('%s/%s' % (outdir, markfname(cname))) as cfile:
            reader = csv.DictReader(cfile)
            assert list(reader.fieldnames)[:5] == ['', 'Top', 'p.value', 'FDR', 'summary.logFC']  # summary.logFC is the log-fold change from the comparison with the lowest p-value (not necessarily the min/max log fold change)
            assert list(reader.fieldnames)[5:] == ['logFC.%d'%i for i in other_clusters]  # should be a column for each pairwise comparison with another cluster
            for il, line in enumerate(reader):
                gene = line['']
                logfc_vals = {i : float(line['logFC.%d'%i]) for i in other_clusters}
                summary_cmarkers['%d-summary'%cname].append((gene, float(line['summary.logFC'])))
                for c2 in logfc_vals:
                    pairwise_cmarkers['%d-%d'%(cname, c2)].append((gene, logfc_vals[c2]))
    for ckey in pairwise_cmarkers:
        pairwise_cmarkers[ckey] = collections.OrderedDict(sorted(pairwise_cmarkers[ckey], key=operator.itemgetter(1), reverse=True))
    for ckey in summary_cmarkers:
        summary_cmarkers[ckey] = collections.OrderedDict(sorted(summary_cmarkers[ckey], key=operator.itemgetter(1), reverse=True))

    # reference marker genes
    fabfo, waickfo = read_ref_data()

    print('  interpretation: "this cluster is much more <type>-like than <clusters>, based on relative upregulation of <N genes>"')
    print('        type    any (N genes)   vs. single clusters                                                     gene contributions (sum over clusters)')
    for cname in cluster_ints:
        print('  %s' % utils.color('green', 'cluster %d' % cname))
        for vtype in waickfo:
            clprods = []
            all_contribs = {}
            for ic2, c2 in enumerate([c for c in cluster_ints if c != cname]):
                dprod, gene_contribs = gexdot(waickfo[vtype], pairwise_cmarkers['%d-%d'%(cname, c2)], return_gene_contributions=True, lbstr='%8s %s '%((vtype+':') if ic2==0 else '', utils.color('blue', str(c2)))) #, debug=True)
                if dprod < min_dprod:
                    continue
                clprods.append({'c2' : c2, 'dprod' : dprod, 'gene_contribs' : gene_contribs})
                for tg, contr in gene_contribs.items():
                    if tg not in all_contribs:
                        all_contribs[tg] = 0.
                    all_contribs[tg] += gene_contribs[tg]
            clprods = sorted(clprods, key=lambda x: x['dprod'], reverse=True)
            anydprod, anygcontribs = gexdot(waickfo[vtype], summary_cmarkers['%d-summary'%cname], return_gene_contributions=True)  # lbstr=XXX
            sumclprod = {'dprod' : anydprod, 'gene_contribs' : anygcontribs}
            if debug and len(clprods) > 0:
                def dcol(d):
                    if d['dprod'] > 0.1:
                        return 'red'
                    elif d['dprod'] > 0.01:
                        return 'yellow'
                    else:
                        return None
                def dpstr(d): return utils.color(dcol(d), '%.3f'%d['dprod'])
                def cstr(d): return utils.color('blue', '%d' % d['c2'])
                tmpstr = '  '.join('%s %s' % (cstr(d), dpstr(d)) for d in clprods)
                anystr = ''
                if sumclprod['dprod'] > min_dprod:
                    anystr = '%s (%2d)' % (dpstr(sumclprod), len(sumclprod['gene_contribs']))
                print('      %s  %-s    %-s  %s' % (utils.color('purple', vtype, width=8),
                                                   # utils.color('blue', ' '.join('%d'%d['c2'] for d in clprods), width=20, padside='right'),
                                                   anystr + ' ' * (12 - utils.len_excluding_colors(anystr)),
                                                   tmpstr + ' ' * (70 - utils.len_excluding_colors(tmpstr)),
                                                   '  '.join('%s %.1f'%(g.lower(), c)for g, c in sorted(list(all_contribs.items()), key=operator.itemgetter(1), reverse=True)),
                ))

    # return barcode_vals, rotation_vals, umap_vals, cluster_vals

# ----------------------------------------------------------------------------------------
def get_init_cmds(feature_matrix_path):
    rcmds = [loadcmd(l) for l in ['DropletUtils', 'scater', 'scran', 'pheatmap', 'celldex', 'SingleR', 'GSEABase', 'AUCell']]  # i think some of these are only required for stuff in the next fcn, but i'm not checking
    rcmds += [
        'options(width=1000)',
        'sce <- read10xCounts("%s")' % feature_matrix_path,
        'rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)',
        # quality control
        'is.mito <- grepl("^MT-", rownames(sce))',  # figure out which genes are mitochondrial
        'qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))',
        'filtered <- quickPerCellQC(qcstats, percent_subsets="subsets_Mito_percent")',  # identifies + removes outliers (in several qc metrics)
        'sce <- sce[, !filtered$discard]',

        # # get reference labels from celldex (so we can remove HSCs) NOTE turning this off since it doesn't really change anything
        # 'ref <- celldex::BlueprintEncodeData()',  # get reference labels from cache or download
        # 'pred <- SingleR(test=sce, ref=ref, labels=ref$label.main)',  # assign labels to our cells (more SingleR detail here: https://ltla.github.io/SingleRBook)
        # 'table(pred$labels)',
        # 'sce <- sce[, pred$labels=="B-cells"]',  # discard non-b-cells
        # 'pred <- pred[pred$labels=="B-cells", ]',

    ]
    return rcmds

# ----------------------------------------------------------------------------------------
def write_gene_counts(feature_matrix_path, outdir, gene_name):  # would be nice to handle more than one gene, but i'd have to spend more time googling/chatgpt'ing R stuff
    workdir = '%s/work' % outdir
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    rcmds = get_init_cmds(feature_matrix_path)
    rcmds += ['sce <- logNormCounts(sce)']  # normalization
    rcmds += ['gene_index <- which(rownames(assay(sce)) == \"%s\")' % gene_name,
              'gene_counts <- assay(sce)[gene_index, ]',
              'output_data <- data.frame(Cell = colData(sce)[2], %s_counts = gene_counts)' % gene_name,
              'write.csv(output_data, file = \"%s/%s-counts.csv\", row.names = FALSE)' % (outdir, gene_name),
              ]
    utils.run_r(rcmds, workdir, logfname='%s/out'%outdir, dryrun=False)

# ----------------------------------------------------------------------------------------
def run_gex(feature_matrix_path, mname, outdir, make_plots=True):
    allowed_mnames = ['hvg', 'fabio', 'waick', 'msigdb']
    workdir = '%s/work' % outdir
    if not os.path.exists(workdir):
        os.makedirs(workdir)
    rcmds = get_init_cmds(feature_matrix_path)
    rcmds += ['capture.output(colData(sce)$Barcode, file="%s/%s")' % (outdir, barcodefname),
              # normalization
              'sce <- logNormCounts(sce)',
              ]
    if mname == 'hvg':  # hvg (hypervariable genes): this is the dumb/default/no prior info choice, where you use the ~700 most variable genes in our sample
        rcmds += [
            'dec <- modelGeneVar(sce)',
            'hvg <- getTopHVGs(dec, prop=0.1)',  # 0.1 gives ~700 most variable genes (if you combine these with the fabio/waick, these totally dominate everything, presumably because there's so many)
        ]
    elif mname == 'fabio':  # he gaves us the 200 most discriminative genes between pb and naive; here i use only the ones up'd in pb
        rcmds += [
            'fabio.markers <- read.csv("%s", sep="\t", header=T)' % fabio_fname, # $name  # genes from fabio (200 most up- or down-regulated in plasmablast as compared to naive B cells)
        ]
    elif mname == 'waick':  # he gave us the 10 most upregulated genes in his samples for each of naive, memory, pb, and prepb, so 40 total
        rcmds += [
            'waick.markers <- read.csv("%s", header=T)' % waickfname,  # 10 most up'd genes for naive, memory, pb, and prepb (40 total). Not sure if it's with respeect to each other, or other cells, or what
        ]
    elif mname == 'msigdb':  # msigdb: I searched through the msigdb "C7" immune cell sets for anything with plasma{blast,cell} and picked the sets that seemed most relevant, ended up with ~1000 genes from these two papers (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4517294/ and https://pubmed.ncbi.nlm.nih.gov/23613519/)
        rcmds += [
            'msigdb.markers <- read.csv("%s", header=T)' % msigdbfname,  # see <msdsets> above -- I just searched through the G7 sets for ones with plasma{blast,cell} and took the nearby ones
        ]
    else:
        raise Exception('mname must be among %s (but got %s)' % (' '.join(allowed_mnames), mname))
        # 'all_genes <- c(fabio.markers$gene, waick.markers$gene, hvg)',  # don't do this, the hvgs overwhelm everything

    mname_markers = mname
    if mname != 'hvg':
        mname_markers += '.markers$gene'
    rcmds += dimredcmds(outdir, mname_markers)
    # reference labels from celldex
    rcmds += [
        'ref <- celldex::BlueprintEncodeData()',  # get reference labels from cache or download
        'pred <- SingleR(test=sce, ref=ref, labels=ref$label.main)',  # assign labels to our cells (more SingleR detail here: https://ltla.github.io/SingleRBook)
        'table(pred$labels)',
    ]
    rcmds += rplotcmds(outdir, 'celldex-label-heatmap', 'plotScoreHeatmap(pred)')
    # only if we have clusters:
    rcmds += ['tab <- table(Assigned=pred$pruned.labels, Cluster=colLabels(sce))',]  # table (and then heatmap) comparing these new labels to our existing clusters
    rcmds += rplotcmds(outdir, 'celldex-label-vs-cluster-heatmap', 'pheatmap(log2(tab+10), color=colorRampPalette(c("white", "blue"))(101))')  # this will crash if you've filtered to only B cells
    if mname != 'hvg':  # doesn't @#%$!$ing work any more (failes with "unable to find an inherited method for function AUCell_buildRankings for signature DelayedMatrix")
        rcmds += ctype_ann_cmds(outdir, mname_markers.replace('$gene', ''))

    utils.run_r(rcmds, workdir, logfname='%s/out'%outdir, dryrun=False)
