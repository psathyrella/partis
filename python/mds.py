import os
import string

import utils

# ----------------------------------------------------------------------------------------
def read_kmeans_clusterfile(clusterfname, seqfos, debug=False):

    # holy crap the need for this function [and its consequent form] make me angry

    all_uids = set([sfo['name'] for sfo in seqfos])
    partition = []
    with open(clusterfname) as clusterfile:
        lines = [l.strip() for l in clusterfile.readlines()]
        iline = -1
        while iline < len(lines) - 1:
            iline += 1

            clidline = lines[iline]
            if debug:
                print 'clid  ', clidline
            if clidline[0] != '$' or int(clidline.lstrip('$').strip('`')) != len(partition) + 1:
                raise Exception('couldn\'t convert %s to the expected cluster id %d' % (clidline, len(partition) + 1))
            partition.append([])

            while True:
                if iline + 2 >= len(lines):
                    break
                iline += 1
                uidline = lines[iline]
                if debug:
                    print 'uid   ', uidline
                iline += 1
                floatline = lines[iline]  # some info about the kmean cluster quality i think? don't care a.t.m.
                if debug:
                    print 'float ', floatline

                uids = set([u for u in uidline.split()])
                if len(uids - all_uids) > 0:
                    raise Exception('read unexpected uid[s] \'%s\' from %s' % (' '.join(uids - all_uids), clusterfname))
                all_uids -= uids
                partition[-1] += list(uids)

                floats = [float(istr) for istr in floatline.split()]
                if len(floats) != len(uids):
                    raise Exception('uid line %d and floats line %d have different lengths:\n  %s\n  %s' % (len(uids), len(floats), uidline, floatline))

                if lines[iline + 1] == '':
                    iline += 1
                    break

            # if emptyline != '':
            #     raise Exception('expected empty line but got \'%s\'' % emptyline)

    if len(all_uids) > 0:
        raise Exception('didn\'t read %d expected queries from %s (%s)' % (len(all_uids), clusterfname, ' '.join(all_uids)))

    os.remove(clusterfname)
    return partition


# ----------------------------------------------------------------------------------------
def kmeans_cluster(n_clusters, seqfos, all_qr_seqs, base_workdir, seed, reco_info=None, n_components=None, max_iterations=1000, max_runs=10, debug=False):
    # NOTE duplication in plotting fcn
    workdir = base_workdir + '/mds'
    msafname = workdir + '/msa.fa'
    clusterfname = workdir + '/clusters.txt'
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    utils.align_many_seqs(seqfos, outfname=msafname)

    # build the R cmd file
    cmdlines = [
        'require(bios2mds, quietly=TRUE)',
        'set.seed(%d)' % seed,
        'human <- import.fasta("%s")' % msafname,
        'active <- mat.dif(human, human)',  # mat.dif or mat.dis?
    ]
    if n_components is not None:  # TODO now that you moved mmds to here, you need to move the plotting stuff into here (or, really, get it to output the PCA components and stop doing the plotting in R)
        cmdlines += ['mmds_active <- mmds(active, pc=%d)' % n_components]
    else:
        raise Exception('need to implement')
    cmdlines += [  # TODO iter.max (max iteratoins and nb.run (max runs), also maybe add method (default "euclidean")
        'kmeans.run1 <- kmeans.run(mmds_active$coord, nb.clus=%d, iter.max=%d, nb.run=%d)' % (n_clusters, max_iterations, max_runs),
        # 'kmeans.run1$clusters',
        # 'kmeans.run1$elements',
        'options(width=10000)',
        'capture.output(kmeans.run1$clusters, file="%s")' % clusterfname,

        # sil.score(mat, nb.clus = c(2:13), nb.run = 100, iter.max = 1000,  # run for every possible number of clusters (?)
        #               method = "euclidean")
        # random.msa  # builds a random [...]
    ]

    utils.run_r(cmdlines, workdir, print_time='kmeans')
    partition = read_kmeans_clusterfile(clusterfname, seqfos)

    clusterfos = []
    for cluster in partition:
        cfo = {
            'seqfos' : [{'name' : uid, 'seq' : all_qr_seqs[uid]} for uid in cluster],
            # 'centroid'  # placeholder to remind you that vsearch clustering adds this, but I think it isn't subsequently used
        }
        cfo['cons_seq'] = utils.cons_seq(0.1, unaligned_seqfos=cfo['seqfos'])
        clusterfos.append(cfo)

        # debug = True
        # if debug and reco_info is not None:
        #     print len(cluster)
        #     for uid in cluster:
        #         print '    %s' % utils.color_gene(reco_info[uid][XXX region + '_gene'])  # need to pass in <region> if I want to uncomment this

    os.remove(msafname)
    os.rmdir(workdir)

    return clusterfos

# ----------------------------------------------------------------------------------------
colors = ['red', 'blue', 'forestgreen', 'grey', 'orange', 'green', 'skyblue4', 'maroon', 'salmon', 'chocolate4', 'magenta']

# ----------------------------------------------------------------------------------------
def write_reco_info_group_colors(reco_info, region, seqfos, group_csv_fname, untranslate):
    all_genes = set([reco_info[untranslate(seqfo['name'])][region + '_gene'] for seqfo in seqfos])
    if len(all_genes) > len(colors):
        print '%s more genes %d than colors %d' % (color('yellow', 'warning'), len(all_genes), len(colors))
    all_gene_list = list(all_genes)
    gene_colors = {all_gene_list[ig] : colors[ig % len(colors)] for ig in range(len(all_gene_list))}
    with open(group_csv_fname, 'w') as groupfile:
        for iseq in range(len(seqfos)):
            seqfo = seqfos[iseq]
            gene = reco_info[untranslate(seqfo['name'])][region + '_gene']
            if len(all_genes) == 1 and iseq == 0:  # R code crashes if there's only one group
                gene += '-dummy'
            groupfile.write('"%s","%s","%s"\n' % (seqfo['name'], gene, gene_colors.get(gene, 'black')))

# ----------------------------------------------------------------------------------------
def write_kmeans_group_colors(clusterfos, seqfos, group_csv_fname, untranslate):  # <seqfos> correspond to the family (<plotname>), <clusterfos> is everybody
    idstrs = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
    uid_to_cluster_id_map = {sfo['name'] : idstrs[iclust] for iclust in range(len(clusterfos)) for sfo in clusterfos[iclust]['seqfos']}
    all_cluster_ids = set(uid_to_cluster_id_map.values()) #range(len(clusterfos)))
    if len(all_cluster_ids) > len(colors):
        print '%s more clusters %d than colors %d' % (color('yellow', 'warning'), len(all_cluster_ids), len(colors))
    all_cluster_id_list = list(all_cluster_ids)
    cluster_id_colors = {all_cluster_id_list[ig] : colors[ig % len(colors)] for ig in range(len(all_cluster_id_list))}
    with open(group_csv_fname, 'w') as groupfile:
        for iseq in range(len(seqfos)):
            seqfo = seqfos[iseq]
            cluster_id = uid_to_cluster_id_map[untranslate(seqfo['name'])]
            if len(all_cluster_ids) == 1 and iseq == 0:  # R code crashes if there's only one group
                cluster_id += '-dummy'
            groupfile.write('"%s","%s","%s"\n' % (seqfo['name'], cluster_id, cluster_id_colors.get(cluster_id, 'black')))

# ----------------------------------------------------------------------------------------
def plot(plotname, seqfos, clusterfos, n_components, plotdir, base_workdir, seed, reco_info=None, region=None, debug=False):
    if not os.path.exists(plotdir):
        os.makedirs(plotdir)

    # R does some horrible truncation or some bullshit when it reads the group csv
    chmap = ['0123456789', 'abcdefghij']
    translations = string.maketrans(*chmap)
    reverse_translations = string.maketrans(*reversed(chmap))
    def translate(name):
        return name.translate(translations)
    def untranslate(trans_name):
        return trans_name.translate(reverse_translations)
    # this somewhat wastefully makes a whole new <seqfos>, but it's better than modifying the original one, and it should only happen when we're plotting
    seqfos = [{'name' : translate(sfo['name']), 'seq' : sfo['seq']} for sfo in seqfos]

    # NOTE duplication in clustering fcn
    workdir = base_workdir + '/mds'
    msafname = workdir + '/msa.fa'
    group_csv_fname = workdir + '/groups.csv'
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    utils.align_many_seqs(seqfos, outfname=msafname)

# ----------------------------------------------------------------------------------------
    # write_reco_info_group_colors(reco_info, region, seqfos, group_csv_fname, untranslate)
    write_kmeans_group_colors(clusterfos, seqfos, group_csv_fname, untranslate)
# ----------------------------------------------------------------------------------------

    cmdlines = [
        'require(bios2mds, quietly=TRUE)',
        'set.seed(%d)' % seed,
        'human <- import.fasta("%s")' % msafname,

        # mat.dif or mat.dis?
        'active <- mat.dif(human, human)',
        'mmds_active <- mmds(active, pc=%d, group.file=%s)' % (n_components, 'NULL' if reco_info is None else '"' + group_csv_fname + '"'),

        'scree.plot(mmds_active$eigen.perc, lab=TRUE, title="%s", pdf.file="%s/scree.pdf")' % (plotname, plotdir),
        'mmds.2D.plot(mmds_active, title="%s", outfile.name="%s/mmds-2d", outfile.type="pdf")' % (plotname, plotdir),

        # random.msa  # builds a random [...]
    ]

    utils.run_r(cmdlines, workdir, print_time='mds plot')

    os.remove(msafname)
    if reco_info is not None:
        os.remove(group_csv_fname)
    os.rmdir(workdir)
