#!/usr/bin/env python3
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import time
import yaml
import itertools
import glob
import argparse
import copy
import random
import os
import tempfile
import subprocess
import sys
import colored_traceback.always
from collections import OrderedDict
from io import open
import ete3

# ----------------------------------------------------------------------------------------
def pairkey(name1, name2):
    return '-&-'.join(sorted([name1, name2]))

# ----------------------------------------------------------------------------------------
scolors = {
    'novel' : '#ffc300',  # 'Gold'
    'data' : 'LightSteelBlue',
    'pale-green' : '#85ad98',
    'pale-blue' : '#94a3d1',
    'tigger-default' : '#d77c7c', #'#c32222',  # red
    'igdiscover' : '#85ad98', #'#29a614',  # green
    'partis' : '#94a3d1', #'#2455ed',  # blue
}

listfaces = [
    'red',
    'blue',
    'green',
]
used_colors, used_faces = {}, {}
simu_colors = OrderedDict((
    ('ok', 'DarkSeaGreen'),
    ('missing', '#d77c7c'),
    ('spurious', '#a44949'),
))

# ----------------------------------------------------------------------------------------
def set_colors(gl_sets, ref_label=None, mix_primary_colors=False):
    listcolors = [plotting.getgrey('medium') for _ in range(10)]
    if ref_label is not None:  # simulation
        for status, color in simu_colors.items():
            scolors[status] = color
        return

    names = sorted(gl_sets.keys())

    if len(names) == 1:  # single-sample data
        scolors[names[0]] = scolors['data']
        return

    assert len(names) in [2, 3]

    if len(names) == 2:
        scolors['all'] = plotting.getgrey('light')
    else:
        scolors['all'] = plotting.getgrey('white')

    for name in names:
        if name not in scolors:
            scolors[name] = listcolors[names.index(name) % len(listcolors)]
            facestr = listfaces[names.index(name) % len(listfaces)]
            used_colors[name] = scolors[name]
            used_faces[name] = facestr
    for name1, name2 in itertools.combinations(names, 2):
        if len(gl_sets) == 2:
            shade = 'white'
        else:
            shade = 'medium' if len(names) == 2 else 'light-medium'
        scolors[pairkey(name1, name2)] = plotting.getgrey(shade)

# ----------------------------------------------------------------------------------------
def get_cmdfos(cmdstr, workdir, outfname):
    return [{'cmd_str' : cmdstr,
             'workdir' : workdir,
             'outfname' : outfname}]

# ----------------------------------------------------------------------------------------
def make_tree(all_genes, workdir, use_cache=False):
    aligned_fname = workdir + '/all-aligned.fa'
    raxml_label = 'xxx'
    raxml_output_fnames = ['%s/RAxML_%s.%s' % (workdir, fn, raxml_label) for fn in ['parsimonyTree', 'log', 'result', 'info', 'bestTree']]
    treefname = [fn for fn in raxml_output_fnames if 'result' in fn][0]
    if use_cache:  # don't re-run muxcle & raxml, just use the previous run's output tree file
        return treefname
    utils.prep_dir(workdir, wildlings=['*.' + raxml_label, os.path.basename(aligned_fname), 'out', 'err', os.path.basename(aligned_fname) + '.reduced'])

    # write and align an .fa with all alleles from any gl set
    start = time.time()
    with tempfile.NamedTemporaryFile(mode='w') as tmpfile:
        for name, seq in all_genes.items():
            tmpfile.write('>%s\n%s\n' % (name, seq))
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        cmdstr = '%s -in %s -out %s' % (args.muscle_path, tmpfile.name, aligned_fname)
        if args.debug:
            print('    %s %s' % (utils.color('red', 'run'), cmdstr))
        utils.run_cmds(get_cmdfos(cmdstr, workdir, aligned_fname), ignore_stderr=True)

    # get a tree for the aligned .fa
    cmdstr = '%s -mGTRCAT -n%s -s%s -p1 -w%s' % (args.raxml_path, raxml_label, aligned_fname, workdir)
    if args.debug:
        print('    %s %s' % (utils.color('red', 'run'), cmdstr))
    utils.run_cmds(get_cmdfos(cmdstr, workdir, treefname), ignore_stderr=True)
    print('    raxml time: %.1f' % (time.time() - start))

    os.remove(aligned_fname)  # rm muscle output
    for fn in [f for f in raxml_output_fnames if f != treefname]:  # rm all the raxml outputs except what the one file we really want
        os.remove(fn)

    return treefname

# ----------------------------------------------------------------------------------------
def getstatus(gene_categories, node, ref_label=None, debug=False):
    gene = node.name

    if not node.is_leaf():
        return 'internal'
    cats = [cat for cat, genes in gene_categories.items() if gene in genes]
    if len(cats) == 0:
        raise Exception('[probably need to bust plot cache/rewrite tree file] couldn\'t find a category for %s among:\n %s' % (node.name, '\n  '.join(['%s:\n     %s' % (k, ' '.join(gene_categories[k])) for k in gene_categories])))
    elif len(cats) > 1:
        raise Exception('wtf?')
    if debug:
        print('%-50s   %s' % (gene, cats[0]))
    return cats[0]

# ----------------------------------------------------------------------------------------
def print_results(gene_categories, gl_sets, ref_label=None):
    pwidth = str(max([len(n) for n in gene_categories]))
    for name, genes in gene_categories.items():
        if name not in scolors:
            raise Exception('status \'%s\' not in scolors' % name)
        if name == 'ok':
            genestr = ''
        else:
            genestr = ' '.join([utils.color_gene(g) for g in genes])
        print(('    %-' + pwidth + 's') % name, end=' ')
        # print '%20s' % scolors[name],
        if name in gl_sets:
            print('  total %2d' % len(gl_sets[name]), end=' ')
        else:
            print('          ', end=' ')
        only_str = 'only' if ref_label is None else ''
        if len(genes) == 0:
            print('  %s %s' % (only_str, utils.color('blue', 'none')))
        else:
            print('  %s %2d    %s' % (only_str, len(genes), genestr))

# ----------------------------------------------------------------------------------------
def write_results(outdir, gene_categories, gl_sets):
    with open(outdir + '/results.yaml', 'w') as yamlfile:
        yamlfo = {gcat : list(genes) for gcat, genes in gene_categories.items()}
        yaml.dump(yamlfo, yamlfile, width=150)

# ----------------------------------------------------------------------------------------
def get_gene_sets(glsfnames, glslabels, ref_label=None, classification_fcn=None, debug=False):
    # debug = True
    glfos = {}
    for label, fname in zip(glslabels, glsfnames):
        if os.path.isdir(fname):
            raise Exception('directory passed instead of germline file name: %s' % fname)
        if os.path.basename(os.path.dirname(fname)) != args.locus:
            raise Exception('unexpected germline directory structure (should have locus \'%s\' at end): %s' % (args.locus, fname))
        gldir = os.path.dirname(fname).replace('/' + args.locus, '')
        glfos[label] = glutils.read_glfo(gldir, args.locus)

    if args.region != 'v':
        print('  not synchronizing gl sets for %s' % args.region)
    if args.region == 'v':  # don't want to deal with d and j synchronization yet
        # synchronize to somebody -- either simulation (<ref_label>) or the first one
        if ref_label is not None:
            sync_label = ref_label
        elif 'partis' in glslabels:
            sync_label = 'partis'
        else:
            sync_label = glslabels[0]
        for label in [l for l in glslabels if l != sync_label]:
            if debug:
                print('  synchronizing %s names to match %s' % (label, sync_label))
            glutils.synchronize_glfos(ref_glfo=glfos[sync_label], new_glfo=glfos[label], region=args.region, ref_label=sync_label, debug=debug)

    gl_sets = {label : {g : seq for g, seq in glfos[label]['seqs'][args.region].items()} for label in glfos}
    all_genes = {g : s for gls in gl_sets.values() for g, s in gls.items()}

    if classification_fcn is not None:
        all_primary_versions = set([classification_fcn(g) for g in all_genes])
        gl_sets = {pv : {label : {g : gl_sets[label][g] for g in gl_sets[label] if classification_fcn(g) == pv} for label in gl_sets} for pv in all_primary_versions}
        all_genes = {pv : {g : s for g, s in all_genes.items() if classification_fcn(g) == pv} for pv in all_primary_versions}

    if len(gl_sets) > 3:
        raise Exception('not implemented')
    gcats = OrderedDict()

    gcats['all'] = set()
    if len(gl_sets) > 2:
        gcats['all'] = set(all_genes)  # gcats['all'] is genes that are in *every* gl set, whereas <all_genes> is genes that're in *any* of 'em (i know, i know...)
        for genes in gl_sets.values():
            gcats['all'] &= set(genes)

    for name, genes in gl_sets.items():
        gcats[name] = set(genes) - gcats['all']

    for ds_1, ds_2 in itertools.combinations(gl_sets, 2):
        gcats[ds_1] -= set(gl_sets[ds_2])
        gcats[ds_2] -= set(gl_sets[ds_1])
        gcats[pairkey(ds_1, ds_2)] = set(gl_sets[ds_2]) & set(gl_sets[ds_1]) - gcats['all']

    if ref_label is not None:
        assert len(gl_sets) == 2
        assert ref_label in gl_sets
        inf_label = [ds for ds in gl_sets if ds != ref_label][0]
        gcats['missing'] = gcats[ref_label]
        gcats['spurious'] = gcats[inf_label]
        gcats['ok'] = gcats[pairkey(ref_label, inf_label)]
        del gcats[ref_label]
        del gcats[inf_label]
        del gcats[pairkey(ref_label, inf_label)]

    any_check = []
    for key, genes in gcats.items():
        any_check += genes
    if sorted(any_check) != sorted(all_genes):
        raise Exception('you done messed up')

    if len(gl_sets) <= 2:
        del gcats['all']

    return all_genes, gl_sets, gcats

# ----------------------------------------------------------------------------------------
def set_node_style(node, status, n_gl_sets, ref_label=None):
    if status != 'internal':
        if status not in scolors:
            raise Exception('status \'%s\' not in scolors' % status)
        node.img_style['bgcolor'] = scolors[status]
        if status not in used_colors:
            used_colors[status] = scolors[status]

        if glutils.is_novel(node.name):
            node.add_face(ete3.CircleFace(args.novel_dot_size, scolors['novel']), column=1) #, position='float') # if args.leaf_names else 'branch')

    # linewidth = 2
    # node.img_style['hz_line_width'] = linewidth
    # node.img_style['vt_line_width'] = linewidth

    names = status.split('-&-')
    if node.is_leaf():
        if args.pie_chart_faces and len(names) > 1:
            pcf = ete3.PieChartFace(percents=[100./len(names) for _ in range(len(names))], width=args.leafheight, height=args.leafheight, colors=[scolors[n] for n in names], line_color=None)
            # pcf = ete3.StackedBarFace(percents=[100./len(names) for _ in range(len(names))], width=30, height=50, colors=[scolors[n] for n in names], line_color=None)
            node.add_face(pcf, column=0, position='aligned')
        elif len(names) == 1 and names[0] in used_faces:
            node.add_face(ete3.RectFace(width=5, height=args.leafheight, bgcolor=used_faces[names[0]], fgcolor=None), column=0, position='aligned')
        elif n_gl_sets > 2:
            rectnames = [n for n in names if n in used_faces]
            node.add_face(ete3.StackedBarFace(percents=[100./len(names) for _ in range(len(rectnames))], width=5 * len(rectnames), height=args.leafheight, colors=[used_faces[rn] for rn in rectnames], line_color=None), column=0, position='aligned')
        else:  # every leaf has to have a face, so that every leaf takes up the same vertical space
            node.add_face(ete3.RectFace(width=1, height=args.leafheight, bgcolor=None, fgcolor=None), column=0, position='aligned')

# ----------------------------------------------------------------------------------------
def get_entirety_of_gene_family(root, family):
    # return set([leaf.name for leaf in node.get_tree_root() if utils.gene_family(leaf.name) == gene_family])
    return set([leaf.name for leaf in root if utils.gene_family(leaf.name) == family])

# ----------------------------------------------------------------------------------------
def set_distance_to_zero(node, debug=False):
    if node.is_root():
        return True
    if node.is_leaf():
        if len(get_entirety_of_gene_family(node.get_tree_root(), utils.gene_family(node.name))) == 1:  # if this is the *only* one from this family
            if debug:
                print('  family %s is of length 1 %s (set to zero)' % (utils.gene_family(node.name), node.name))
            return True
        else:
            return False

    descendents = set([leaf.name for leaf in node])
    gene_families = set([utils.gene_family(d) for d in descendents])
    if debug:
        print('  %s' % ' '.join([utils.shorten_gene_name(d) for d in descendents]))
        print('      %s' % ' '.join(gene_families))
    if len(gene_families) == 0:
        raise Exception('zero length gene family set from %s' % ' '.join([leaf.name for leaf in node]))
    if len(gene_families) > 1:
        return True

    gene_family = list(gene_families)[0]
    entirety_of_gene_family = get_entirety_of_gene_family(node.get_tree_root(), gene_family)
    if debug:
        if len(entirety_of_gene_family - descendents) > 0:
            print('    missing %d / %d of family' % (len(entirety_of_gene_family - descendents), len(entirety_of_gene_family)))
        elif len(descendents - entirety_of_gene_family) > 0:
            raise Exception('wtf should\'ve already returned')
        else:
            print('    setting to zero')
    return descendents == entirety_of_gene_family

# ----------------------------------------------------------------------------------------
def write_legend(plotdir):
    def get_leg_name(status):
        if args.legends is not None and status in args.glslabels:
            lname = args.legends[args.glslabels.index(status)]
        elif status == 'both':
            if len(args.glsfnames) == 2:
                lname = 'both'
            elif len(args.glsfnames) == 3:
                lname = 'two'
            else:
                raise Exception('can\'t make a legend when --glsfnames is length %d' % len(args.glsfnames))
        elif status == 'all':
            if len(args.glsfnames) == 2:
                lname = 'both'
            elif len(args.glsfnames) == 3:
                lname = 'all three'
            else:
                raise Exception('can\'t make a legend when --glsfnames is length %d' % len(args.glsfnames))
        else:
            lname = status
        return lname
    def add_stuff(status, leg_name, color):
        legfo[leg_name] = color
        if status in used_faces:
            facefo[leg_name] = used_faces[status]

    legfo, facefo = {}, {}
    if args.ref_label is not None:
        for status, color in simu_colors.items():
            add_stuff(status, status, color)
    else:
        added_two_method_color = False
        for status, color in used_colors.items():
            if '-&-' in status:
                for substatus in status.split('-&-'):  # arg, have to handle cases where the single one isn't in there
                    if get_leg_name(substatus) not in legfo:
                        add_stuff(substatus, get_leg_name(substatus), scolors[substatus])
                if not added_two_method_color:
                    leg_name = get_leg_name('both')
                    added_two_method_color = True
                else:
                    continue
            else:
                leg_name = get_leg_name(status)

            add_stuff(status, leg_name, color)

    # figure out the order we want 'em in
    lnames = sorted(legfo.keys())
    for status in ['both', 'all']:
        if get_leg_name(status) in lnames:
            lnames.remove(get_leg_name(status))
            lnames.append(get_leg_name(status))

    etree = ete3.ClusterTree() #'(a);')
    tstyle = ete3.TreeStyle()
    tstyle.show_scale = False
    # tstyle.show_leaf_name = False
    # for node in etree.traverse():
    #     print node.name
    #     node.add_face(ete3.CircleFace(args.novel_dot_size, scolors['novel']), column=1) #, position='float') # if args.leaf_names else 'branch')

    dummy_column = 0
    pic_column = 1
    text_column = 2
    leg_title_height = 1.5 * args.leafheight  # if args.legend_title is not None else 0.75 * args.leafheight

    for icol in range(text_column + 1):  # add a top border
        tstyle.title.add_face(ete3.RectFace(0.9*args.leafheight, 0.9*args.leafheight, fgcolor=None, bgcolor=None), column=icol)

    tstyle.title.add_face(ete3.TextFace(' ', fsize=leg_title_height), column=dummy_column)  # adds a left border

    if args.legend_title is not None:
        tstyle.title.add_face(ete3.TextFace('', fsize=leg_title_height), column=pic_column)  # keeps the first legend entry from getting added on this line
        tstyle.title.add_face(ete3.TextFace(args.legend_title, fsize=leg_title_height, fgcolor='black', bold=True), column=text_column)  # add an empty title so there's some white space at the top, even with no actual title text

    for leg_name in lnames:
        color = legfo[leg_name]
        size_factor = 2.
        if leg_name in facefo:
            tstyle.title.add_face(ete3.StackedBarFace([80., 20.], width=size_factor*args.leafheight, height=size_factor*args.leafheight, colors=[color, facefo[leg_name]], line_color='black'), column=pic_column)  # looks like maybe they reversed fg/bg kwarg names
        else:
            tstyle.title.add_face(ete3.RectFace(size_factor*args.leafheight, size_factor*args.leafheight, fgcolor='black', bgcolor=color), column=pic_column)  # looks like maybe they reversed fg/bg kwarg names
        tstyle.title.add_face(ete3.TextFace(' ' + leg_name, fsize=args.leafheight, fgcolor='black'), column=text_column)

    tstyle.title.add_face(ete3.CircleFace(1.5*args.novel_dot_size, scolors['novel']), column=pic_column)
    tstyle.title.add_face(ete3.TextFace('novel allele', fsize=args.leafheight), column=text_column)  # keeps the first legend entry from getting added on this line

    etree.render(plotdir + '/legend.svg', tree_style=tstyle)

# ----------------------------------------------------------------------------------------
def draw_tree(plotdir, plotname, treestr, gl_sets, all_genes, gene_categories, ref_label=None, arc_start=None, arc_span=None):
    etree = ete3.ClusterTree(treestr)
    node_names = set()  # make sure we get out all the genes we put in
    for node in etree.traverse():
        if set_distance_to_zero(node):
            node.dist = 0. if ref_label is not None else 1e-9  # data crashes sometimes with float division by zero if you set it to 0., but simulation sometimes gets screwed up for some other reason (that I don't understand) if it's 1e-9
        # node.dist = 1.
        status = getstatus(gene_categories, node, ref_label=ref_label)
        set_node_style(node, status, len(gl_sets), ref_label=ref_label)
        if node.is_leaf():
            node_names.add(node.name)
    if len(set(all_genes) - node_names) > 0:
        raise Exception('missing genes from final tree: %s' % ' '.join(node_names))

    if args.param_dirs is not None:
        countfo = OrderedDict()
        for label, pdir in zip(args.glslabels, args.param_dirs):  # it would be cleaner to do this somewhere else
            if pdir == 'None':  # not the best way to do this
                continue
            countfo[label] = utils.read_overall_gene_probs(pdir, normalize=True)[args.region]
        for node in etree.traverse():
            node.countstr = '%s' % ' '.join([('%.2f' % (100 * cfo[node.name])) if node.name in cfo else '-' for cfo in countfo.values()])

    if ref_label is None:  # have to do it in a separate loop so it doesn't screw up the distance setting
        for node in [n for n in etree.traverse() if n.is_leaf()]:  # yeah I'm sure there's a fcn for that
            node.name = utils.shorten_gene_name(node.name)

    tstyle = ete3.TreeStyle()
    tstyle.show_scale = False

    if len(args.glslabels) > 1:
        write_legend(plotdir)
    if args.title is not None:
        fsize = 13
        tstyle.title.add_face(ete3.TextFace(args.title, fsize=fsize, bold=True), column=0)
        if args.title_color is not None:
            # tstyle.title.add_face(ete3.CircleFace(fsize, scolors[args.title]), column=1)
            tcol = scolors[args.title_color] if args.title_color in scolors else args.title_color
            rect_width = 3 if len(args.title) < 12 else 2
            tstyle.title.add_face(ete3.RectFace(width=rect_width*fsize, height=fsize, bgcolor=tcol, fgcolor=None), column=1)
    suffix = '.svg'
    imagefname = plotdir + '/' + plotname + suffix
    print('      %s' % imagefname)
    etree.render(utils.insert_before_suffix('-leaf-names', imagefname), tree_style=tstyle)
    tstyle.show_leaf_name = False
    etree.render(imagefname, tree_style=tstyle)

    # NOTE all the node names are screwed up after this, so you'll have to fix them if you add another step
    if args.param_dirs is not None:
        for node in etree.traverse():
            node.name = node.countstr
        tstyle.show_leaf_name = True
        etree.render(utils.insert_before_suffix('-gene-counts', imagefname), tree_style=tstyle)

# ----------------------------------------------------------------------------------------
def plot_trees(args, plotdir, plotname, glsfnames, glslabels):
    all_genes, gl_sets, gene_categories = get_gene_sets(glsfnames, glslabels, ref_label=args.ref_label)
    set_colors(gl_sets, ref_label=args.ref_label)
    print_results(gene_categories, gl_sets, ref_label=args.ref_label)
    write_results(plotdir, gene_categories, gl_sets)
    if args.only_print:
        return

    treefname = make_tree(all_genes, plotdir + '/workdir', use_cache=args.use_cache)
    with open(treefname) as treefile:
        treestr = treefile.read().strip()

    draw_tree(plotdir, plotname, treestr, gl_sets, all_genes, gene_categories, ref_label=args.ref_label)

# ----------------------------------------------------------------------------------------
example_str = '\n    '.join(['example usage (note that this example as it is will be a) really slow, since the files are the full imgt set, with ~250 genes, and b) not very interesting, since the two .fasta files are the same):',
                             './bin/plot-gl-set-trees.py --glsfnames data/germlines/human/igh/ighv.fasta:data/germlines/human/igh/ighv.fasta --glslabels foo:bar --locus igh'])
parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, epilog=example_str)
parser.add_argument('--plotdir', default=os.getcwd() + '/gl-set-tree-plots')
parser.add_argument('--plotname', default='test')
parser.add_argument('--glsfnames', required=True, help='colon-separated list of germline ig fasta file names')
parser.add_argument('--glslabels', required=True, help='colon-separated list of labels corresponding to --glsfnames')
parser.add_argument('--param-dirs', help='parameter dirs for each gls fname, for getting counts for each gene')
parser.add_argument('--locus', required=True, choices=['igh', 'igk', 'igl'])
parser.add_argument('--legends', help='colon-separated list of legend labels')
parser.add_argument('--legend-title')
parser.add_argument('--pie-chart-faces', action='store_true')
parser.add_argument('--use-cache', action='store_true', help='use existing raxml output from a previous run (crashes if it isn\'t there)')
parser.add_argument('--only-print', action='store_true', help='just print the summary, without making any plots')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--title')
parser.add_argument('--title-color')
parser.add_argument('--region', default='v')
parser.add_argument('--partis-dir', default=os.getcwd(), help='path to main partis install dir')
parser.add_argument('--muscle-path', default='./packages/muscle/muscle3.8.31_i86linux64')
parser.add_argument('--raxml-path', default='./packages/standard-RAxML/raxmlHPC-SSE3')  # laptop: '-AVX'
parser.add_argument('--ref-label', help='label (in --glslabels) corresponding to simulation/truth')

args = parser.parse_args()

sys.path.insert(1, args.partis_dir) # + '/python')
try:
    import partis.utils as utils
    import partis.glutils as glutils
    import partis.plotting as plotting
except ImportError as e:
    print(e)
    raise Exception('couldn\'t import from main partis dir \'%s\' (set with --partis-dir)' % args.partis_dir)

args.glsfnames = utils.get_arg_list(args.glsfnames)
args.glslabels = utils.get_arg_list(args.glslabels)
args.param_dirs = utils.get_arg_list(args.param_dirs)
args.legends = utils.get_arg_list(args.legends)
if not os.path.exists(args.muscle_path):
    raise Exception('muscle binary %s doesn\'t exist (set with --muscle-path)' % args.muscle_path)
if not os.path.exists(args.raxml_path):
    raise Exception('raxml binary %s doesn\'t exist (set with --raxml-path)' % args.raxml_path)
if not os.path.exists(args.plotdir):
    os.makedirs(args.plotdir)
args.leafheight = 10  #20 if args.leaf_names else 10  # arg, kinda messy
args.novel_dot_size = 2.5

assert len(args.glslabels) == len(set(args.glslabels))  # no duplicates

plot_trees(args, args.plotdir, args.plotname, args.glsfnames, args.glslabels)
