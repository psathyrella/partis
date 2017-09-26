#!/usr/bin/env python
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
import itertools
import glob
import argparse
import copy
import random
import os
import tempfile
import subprocess
import sys
import ete3
import colored_traceback.always
from collections import OrderedDict

sys.path.insert(1, './python')
import utils
import glutils
sys.path.insert(1, './datascripts')
import heads

# custom interactive faces: http://etetoolkit.org/docs/latest/tutorial/tutorial_drawing.html#id34

methodnames = ['tigger-default', 'partis', 'igdiscover']

# ----------------------------------------------------------------------------------------
def methstr(meth):
    return 'tigger' if meth == 'tigger-default' else meth

# ----------------------------------------------------------------------------------------
def tpstr(timepoint):
    if timepoint[0] == 'W':
        return 'week %d' % int(timepoint[1:])
    elif timepoint[0] == 'M':
        return 'month %d' % int(timepoint[1:])
    elif 'dpi' in timepoint:
        return timepoint.replace('dpi', ' dpi')
    elif timepoint[0] in ['p', 'm'] and timepoint[-1] in ['h', 'd']:
        plusminusstr = timepoint[0].replace('p', '+').replace('m', '-')
        number = int(timepoint[1:-1])
        unitstr = timepoint[-1].replace('h', 'hour').replace('d', 'day')
        return '%s%d %s%s' % (plusminusstr, number, unitstr, utils.plural(number))
    else:
        raise Exception('not sure what to do with %s' % timepoint)

# ----------------------------------------------------------------------------------------
def get_title(nstr, other_nstr=None):
    if nstr in methodnames:
        return_str = methstr(nstr)
    elif args.metafo is not None and nstr in args.metafo:
        return_str = args.metafo[nstr]['subject']
        extras = []
        if args.metafo[nstr]['isotype'] != '':
            extras.append('Ig%s' % args.metafo[nstr]['isotype'].upper())
        if args.metafo[nstr]['timepoint'] != '':
            extras.append('%s' %tpstr(args.metafo[nstr]['timepoint']))
        if len(extras) > 0:
            return_str += '  (%s)' % '  '.join(extras)
    else:
        return_str = nstr

    return return_str

# ----------------------------------------------------------------------------------------
def getgrey(gtype='medium'):
    if gtype == 'medium':
        return '#929292'
    elif gtype == 'light-medium':
        return '#cdcdcd'
    elif gtype == 'light':
        return '#d3d3d3'
    elif gtype == 'white':
        return '#ffffff'
    else:
        assert False

# ----------------------------------------------------------------------------------------
def pairkey(name1, name2):
    return '-&-'.join(sorted([name1, name2]))

# ----------------------------------------------------------------------------------------
scolors = {
    'ok' : 'DarkSeaGreen',
    'missing' : '#d77c7c',
    'spurious' : '#a44949',
    'data' : 'LightSteelBlue',
    'all' : getgrey('light'),
    'pale-green' : '#85ad98',
    'pale-blue' : '#94a3d1',
    'tigger-default' : '#d77c7c', #'#c32222',  # red
    'igdiscover' : '#85ad98', #'#29a614',  # green
    'partis' : '#94a3d1', #'#2455ed',  # blue
}

listcolors = [
    getgrey('medium'),
    getgrey('medium'),
    # hopefully only really using the first two a.t.m.
    '#d5aaf4',  # purple
    '#87cfcf',  #
    scolors['spurious'],
    '#2455ed',  # blue
    '#29a614',  # green
    '#c32222',  # red
]

listfaces = [
    'red',
    'blue',
]
used_faces = {}

# ----------------------------------------------------------------------------------------
def set_colors(gl_sets, ref_label=None, mix_primary_colors=False):
    if ref_label is not None:  # simulation
        return

    names = sorted(gl_sets.keys())

    if len(names) == 1:  # single-sample data
        scolors[names[0]] = scolors['data']
        return

    for name in names:
        if name not in scolors:
            scolors[name] = listcolors[names.index(name) % len(listcolors)]
            facestr = listfaces[names.index(name) % len(listfaces)]
            used_faces[name] = facestr
    for name1, name2 in itertools.combinations(names, 2):
        if len(names) == 3:
            pair_color = getgrey('medium')
        else:
            pair_color = getgrey('light')
        # if name1 in methodnames or name2 in methodnames:
        #     if len(names) == 3:
        #         pair_color = getgrey('medium')
        #     else:
        #         pair_color = getgrey('light')
        # else:
        #     pair_color = getgrey('white')
        scolors[pairkey(name1, name2)] = pair_color

# ----------------------------------------------------------------------------------------
def get_cmdfos(cmdstr, workdir, outfname):
    return [{'cmd_str' : cmdstr,
             'workdir' : workdir,
             'outfname' : outfname}]

# ----------------------------------------------------------------------------------------
def shorten_name(name):
    if name[:2] != 'IG':
        raise Exception('bad node name %s' % name)

    pv, sv, allele = utils.split_gene(name)
    if glutils.is_novel(name):
        _, template_name, mutstrs = glutils.split_inferred_allele_name(name)
        if mutstrs is None:
            allele = '%s (+...)' % (utils.allele(template_name))
        else:
            allele = '%s (+%d snp%s)' % (utils.allele(template_name), len(mutstrs), utils.plural(len(mutstrs)))
    if sv is not None:
        return '%s-%s*%s' % (pv, sv, allele)
    else:
        return '%s*%s' % (pv, allele)

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
    with tempfile.NamedTemporaryFile() as tmpfile:
        for name, seq in all_genes.items():
            tmpfile.write('>%s\n%s\n' % (name, seq))
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        cmdstr = '%s -in %s -out %s' % (args.muscle_path, tmpfile.name, aligned_fname)
        if args.debug:
            print '    %s %s' % (utils.color('red', 'run'), cmdstr)
        utils.run_cmds(get_cmdfos(cmdstr, workdir, aligned_fname), ignore_stderr=True)

    # get a tree for the aligned .fa
    cmdstr = '%s -mGTRCAT -n%s -s%s -p1 -w%s' % (args.raxml_path, raxml_label, aligned_fname, workdir)
    if args.debug:
        print '    %s %s' % (utils.color('red', 'run'), cmdstr)
    utils.run_cmds(get_cmdfos(cmdstr, workdir, treefname), ignore_stderr=True)

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
        raise Exception('couldn\'t find a category for %s (among %s)' % (node.name, gene_categories))
    elif len(cats) > 1:
        raise Exception('wtf?')
    if debug:
        print '%-50s   %s' % (gene, cats[0])
    return cats[0]

# ----------------------------------------------------------------------------------------
def print_results(gene_categories, ref_label=None):
    pwidth = str(max([len(n) for n in gene_categories]))
    for name, genes in gene_categories.items():
        if name not in scolors:
            raise Exception('status \'%s\' not in scolors' % name)
        if name == 'ok':
            genestr = ''
        else:
            genestr = ' '.join([utils.color_gene(g) for g in genes])
        print ('    %-' + pwidth + 's') % name,
        print '%20s' % scolors[name],
        if len(genes) == 0:
            print ' %s' % utils.color('blue', 'none')
        else:
            print ' %2d    %s' % (len(genes), genestr)

# ----------------------------------------------------------------------------------------
def get_gene_sets(glsfnames, glslabels, ref_label=None, classification_fcn=None, debug=False):
    glfos = {}
    for label, fname in zip(glslabels, glsfnames):
        gldir = os.path.dirname(fname).replace('/' + args.locus, '')
        glfos[label] = glutils.read_glfo(gldir, args.locus)  # this is gonna fail for tigger since you only have the .fa

    # synchronize to somebody -- either simulation (<ref_label>) or the first one
    if ref_label is not None:
        sync_label = ref_label
    else:
        sync_label = glslabels[0]
    for label in [l for l in glslabels if l != sync_label]:
        if debug:
            print '    syncronizing %s names to match %s' % (label, sync_label)
        glutils.synchronize_glfos(ref_glfo=glfos[sync_label], new_glfo=glfos[label], region=args.region, debug=debug)

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
def set_node_style(node, status, n_gl_sets, used_colors, ref_label=None):
    if status != 'internal':
        if status not in scolors:
            raise Exception('status \'%s\' not in scolors' % status)
        node.img_style['bgcolor'] = scolors[status]
        if status not in used_colors:
            used_colors[status] = scolors[status]

        if glutils.is_novel(node.name):
            node.add_face(ete3.CircleFace(2.5, 'Gold'), column=1) #, position='float') # if args.leaf_names else 'branch')

    # linewidth = 2
    # node.img_style['hz_line_width'] = linewidth
    # node.img_style['vt_line_width'] = linewidth

    stlist = status.split('-&-')
    if len(stlist) > 0:
        if len(set(stlist) & set(methodnames)) == len(stlist):  # method vs method
            names = status.split('-&-')
            pcf = ete3.PieChartFace(percents=[100./len(names) for _ in range(len(names))], width=20, height=20, colors=[scolors[n] for n in names], line_color=None)
            # pcf = ete3.StackedBarFace(percents=[100./len(names) for _ in range(len(names))], width=30, height=50, colors=[scolors[n] for n in names], line_color=None)
            node.add_face(pcf, column=0, position='aligned')
    if len(stlist) == 1 and status in used_faces:
        node.add_face(ete3.RectFace(width=5, height=20, bgcolor=used_faces[status], fgcolor=None), column=0, position='aligned')

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
                print '  family %s is of length 1 %s (set to zero)' % (utils.gene_family(node.name), node.name)
            return True
        else:
            return False

    descendents = set([leaf.name for leaf in node])
    gene_families = set([utils.gene_family(d) for d in descendents])
    if debug:
        print '  %s' % ' '.join([shorten_name(d) for d in descendents])
        print '      %s' % ' '.join(gene_families)
    if len(gene_families) == 0:
        raise Exception('zero length gene family set from %s' % ' '.join([leaf.name for leaf in node]))
    if len(gene_families) > 1:
        return True

    gene_family = list(gene_families)[0]
    entirety_of_gene_family = get_entirety_of_gene_family(node.get_tree_root(), gene_family)
    if debug:
        if len(entirety_of_gene_family - descendents) > 0:
            print '    missing %d / %d of family' % (len(entirety_of_gene_family - descendents), len(entirety_of_gene_family))
        elif len(descendents - entirety_of_gene_family) > 0:
            raise Exception('wtf should\'ve already returned')
        else:
            print '    setting to zero'
    return descendents == entirety_of_gene_family

# # ----------------------------------------------------------------------------------------
# def add_legend(tstyle, used_colors):
#     for status, color in used_colors.items():
#         if status == 'all' or '-&-' in status:
#             continue
#         tstyle.title.add_face(ete3.RectFace(20, 20, color, color, label=status), column=0)

# ----------------------------------------------------------------------------------------
def write_legend(used_colors, plotdir):
    def add_stuff(status, leg_name, color):
        legfo[leg_name] = color
        if status in used_faces:
            facefo[leg_name] = used_faces[status]

    added_two_method_color = False
    legfo, facefo = OrderedDict(), OrderedDict()
    for status, color in used_colors.items():
        if '-&-' in status:
            for substatus in status.split('-&-'):  # arg, have to handle cases where the single one isn't in there
                if get_title(substatus) not in legfo:
                    add_stuff(substatus, get_title(substatus), scolors[substatus])
            if not added_two_method_color:
                leg_name = 'both'
                added_two_method_color = True
            else:
                continue
        else:
            leg_name = get_title(status)

        add_stuff(status, leg_name, color)

    # reorder some of 'em
    for leg_name in ['both', 'all']:
        if leg_name in legfo:
            tmpcolor = legfo[leg_name]
            del legfo[leg_name]
            legfo[leg_name] = tmpcolor

    etree = ete3.ClusterTree()
    tstyle = ete3.TreeStyle()
    tstyle.show_scale = False
    for leg_name, color in legfo.items():
        if leg_name in facefo:
            tstyle.title.add_face(ete3.StackedBarFace([85., 15.], width=20, height=20, colors=[color, facefo[leg_name]], line_color='black'), column=0)  # looks like maybe they reversed fg/bg kwarg names
        else:
            tstyle.title.add_face(ete3.RectFace(20, 20, fgcolor='black', bgcolor=color), column=0)  # looks like maybe they reversed fg/bg kwarg names
        tstyle.title.add_face(ete3.TextFace(' ' + leg_name, fgcolor='black'), column=1)

    etree.render(plotdir + '/legend.svg', h=300, tree_style=tstyle)  # w=500, h=500, 

# ----------------------------------------------------------------------------------------
def draw_tree(plotdir, plotname, treestr, gl_sets, all_genes, gene_categories, ref_label=None, arc_start=None, arc_span=None):
    etree = ete3.ClusterTree(treestr)
    node_names = set()  # make sure we get out all the genes we put in
    used_colors = {}
    for node in etree.traverse():
        if set_distance_to_zero(node):
            node.dist = 0. if ref_label is not None else 1e-9  # data crashes sometimes with float division by zero if you set it to 0., but simulation sometimes gets screwed up for some other reason (that I don't understand) if it's 1e-9
        # node.dist = 1.
        status = getstatus(gene_categories, node, ref_label=ref_label)
        set_node_style(node, status, len(gl_sets), used_colors, ref_label=ref_label)
        if node.is_leaf():
            node_names.add(node.name)
    if len(set(all_genes) - node_names) > 0:
        raise Exception('missing genes from final tree: %s' % ' '.join(node_names))

    if ref_label is None:  # have to do it in a separate loop so it doesn't screw up the distance setting
        for node in [n for n in etree.traverse() if n.is_leaf()]:  # yeah I'm sure there's a fcn for that
            node.name = shorten_name(node.name)

    tstyle = ete3.TreeStyle()
    tstyle.show_scale = False
    # if not args.leaf_names:
    #     tstyle.show_leaf_name = False

    # tstyle.mode = 'c'
    # if arc_start is not None:
    #     tstyle.arc_start = arc_start
    # if arc_span is not None:
    #     tstyle.arc_span = arc_span

    # if args.legend:
    #     add_legend(tstyle, used_colors)
    write_legend(used_colors, plotdir)
    if args.title is not None:
        title = get_title(args.title)
        fsize = 13
        tstyle.title.add_face(ete3.TextFace(title, fsize=fsize, bold=True), column=0)
        if args.title in scolors:
            tstyle.title.add_face(ete3.CircleFace(fsize, scolors[args.title]), column=1)
    suffix = '.svg'
    imagefname = plotdir + '/' + plotname + suffix
    print '      %s' % imagefname
    etree.render(imagefname, tree_style=tstyle) # h=750, 

# ----------------------------------------------------------------------------------------
def plot_trees(args, plotdir, plotname, glsfnames, glslabels):
    all_genes, gl_sets, gene_categories = get_gene_sets(glsfnames, glslabels, ref_label=args.ref_label)
    set_colors(gl_sets, ref_label=args.ref_label)
    print_results(gene_categories, ref_label=args.ref_label)

    treefname = make_tree(all_genes, plotdir + '/workdir', use_cache=args.use_cache)
    with open(treefname) as treefile:
        treestr = treefile.read().strip()

    draw_tree(plotdir, plotname, treestr, gl_sets, all_genes, gene_categories, ref_label=args.ref_label)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--plotdir', required=True)
parser.add_argument('--plotname', required=True)
parser.add_argument('--glsfnames', required=True)
parser.add_argument('--glslabels', required=True)
parser.add_argument('--locus', required=True)
parser.add_argument('--legend', action='store_true')
parser.add_argument('--leaf-names', action='store_true')
parser.add_argument('--use-cache', action='store_true', help='just print results and remake the plots, without remaking the tree (which is the slow part)')
parser.add_argument('--debug', action='store_true')
parser.add_argument('--title')
parser.add_argument('--study')
parser.add_argument('--region', default='v')
parser.add_argument('--muscle-path', default='./packages/muscle/muscle3.8.31_i86linux64')
parser.add_argument('--raxml-path', default=glob.glob('./packages/standard-RAxML/raxmlHPC-*')[0])
parser.add_argument('--ref-label')  # label corresponding to simulation

args = parser.parse_args()
args.metafo = None
if args.study is not None:
    args.metafo = heads.read_metadata(args.study)

args.glsfnames = utils.get_arg_list(args.glsfnames)
args.glslabels = utils.get_arg_list(args.glslabels)
if not os.path.exists(args.muscle_path):
    raise Exception('muscle path %s does not exist' % args.muscle_path)
if not os.path.exists(args.raxml_path):
    raise Exception('raxml path %s does not exist' % args.raxml_path)

assert len(args.glslabels) == len(set(args.glslabels))  # no duplicates

plot_trees(args, args.plotdir, args.plotname, args.glsfnames, args.glslabels)
