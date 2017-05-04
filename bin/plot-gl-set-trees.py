#!/usr/bin/env python
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
import glob
import argparse
import copy
import random
import os
import tempfile
import subprocess
import sys
import ete3

sys.path.insert(1, './python')
import utils

def get_cmdfos(cmdstr, workdir, outfname):
    return [{'cmd_str' : cmdstr,
             'workdir' : workdir,
             'outfname' : outfname}]

# ----------------------------------------------------------------------------------------
def make_tree(all_genes, workdir, glsfnames, glslabels, use_cache=False):
    aligned_fname = workdir + '/all-aligned.fa'
    raxml_label = 'xxx'
    raxml_output_fnames = ['%s/RAxML_%s.%s' % (workdir, fn, raxml_label) for fn in ['parsimonyTree', 'log', 'result', 'info', 'bestTree']]
    treefname = [fn for fn in raxml_output_fnames if 'result' in fn][0]
    if use_cache:  # don't re-run muxcle & raxml, just use the previous run's output tree file
        return treefname
    utils.prep_dir(workdir, wildlings=['*.' + raxml_label, os.path.basename(aligned_fname), 'out', 'err'])

    # write and align an .fa with all alleles from any gl set
    with tempfile.NamedTemporaryFile() as tmpfile:
        for name, seq in all_genes.items():
            tmpfile.write('>%s\n%s\n' % (name, seq))
        tmpfile.flush()  # BEWARE if you forget this you are fucked
        cmdstr = '%s -in %s -out %s' % (args.muscle_path, tmpfile.name, aligned_fname)
        print '    %s %s' % (utils.color('red', 'run'), cmdstr)
        utils.run_cmds(get_cmdfos(cmdstr, workdir, aligned_fname), ignore_stderr=True)

    # get a tree for the aligned .fa
    cmdstr = '%s -mGTRCAT -n%s -s%s -p1 -w%s' % (args.raxml_path, raxml_label, aligned_fname, workdir)
    print '    %s %s' % (utils.color('red', 'run'), cmdstr)
    utils.run_cmds(get_cmdfos(cmdstr, workdir, treefname), ignore_stderr=True)

    os.remove(aligned_fname)  # rm muscle output
    for fn in [f for f in raxml_output_fnames if f != treefname]:  # rm all the raxml outputs except what the one file we really want
        os.remove(fn)

    return treefname

# ----------------------------------------------------------------------------------------
def getstatus(gl_sets, node):
    gene = node.name
    if not node.is_leaf():
        return 'internal'
    elif gene in gl_sets['sim'] and gene in gl_sets['inf']:
        return 'ok'
    elif gene in gl_sets['sim']:
        return 'missing'
    elif gene in gl_sets['inf']:
        return 'spurious'
    else:
        assert False

# ----------------------------------------------------------------------------------------
def print_results(gl_sets):
    tmpfo = {'missing' : set(gl_sets['sim']) - set(gl_sets['inf']),
             'spurious' : set(gl_sets['inf']) - set(gl_sets['sim']),
             'ok' : set(gl_sets['inf']) & set(gl_sets['sim'])}
    for name, genes in tmpfo.items():
        print '    %9s %2d: %s' % (name, len(genes), ' '.join([utils.color_gene(g) for g in genes]))

# ----------------------------------------------------------------------------------------
def plot_gls_gen_tree(args, plotdir, plotname, glsfnames, glslabels, leg_title=None, title=None):
    assert len(glslabels) == len(set(glslabels))  # no duplicates
    assert glslabels == ['sim', 'inf']  # otherwise stuff needs to be updated

    # read the input germline sets
    all_genes, gl_sets = {}, {}
    for label, fname in zip(glslabels, glsfnames):
        gl_sets[label] = {gfo['name'] : gfo['seq'] for gfo in utils.read_fastx(fname)}
        for name, seq in gl_sets[label].items():
            if name not in all_genes:
                all_genes[name] = seq

    print_results(gl_sets)

    workdir = plotdir + '/workdir'
    treefname = make_tree(all_genes, workdir, glsfnames, glslabels, use_cache=args.use_cache)
    with open(treefname) as treefile:
        treestr = treefile.read().strip()
    # treestr = "(A:0.7,B:0.7):0.3;"

    scolors = {'ok' : 'DarkSeaGreen', 'missing' : 'IndianRed', 'spurious' : 'IndianRed'}
    faces = {'missing'  : ete3.CircleFace(10, 'white'), 
             'spurious' : ete3.CircleFace(10, 'black')}

    # ----------------------------------------------------------------------------------------
    def set_node_style(node, status):
        linewidth = 4
        if not status == 'internal':
            node.img_style['bgcolor'] = scolors[status]
        if '_' in node.name:  # tigger names
            raise Exception('unhandled gene name %s' % node.name)
        if '+' in node.name:
            linewidth = 8
            node.img_style['hz_line_color'] = 'Gold'
            node.img_style['vt_line_color'] = 'Gold'
        node.img_style['hz_line_width'] = linewidth
        node.img_style['vt_line_width'] = linewidth

    etree = ete3.ClusterTree(treestr)
    node_names = set()  # make sure we get out all the genes we put in
    for node in etree.traverse():
        node.dist = 1
        status = getstatus(gl_sets, node)
        set_node_style(node, status)
        if node.is_leaf():
            if status in faces:
                node.add_face(copy.deepcopy(faces[status]), column=0)
            node_names.add(node.name)
    if len(set(all_genes) - node_names) > 0:
        raise Exception('missing genes from final tree: %s' % ' '.join(node_names))

    tstyle = ete3.TreeStyle()
    tstyle.show_leaf_name = False
    tstyle.mode = 'c'
    tstyle.show_scale = False
    etree.render(plotdir + '/' + plotname + '.svg', h=750, tree_style=tstyle)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--plotdir', required=True)
parser.add_argument('--plotname', required=True)
parser.add_argument('--glsfnames', required=True)
parser.add_argument('--glslabels', required=True)
parser.add_argument('--use-cache', action='store_true')
parser.add_argument('--title')
parser.add_argument('--muscle-path', default='./packages/muscle/muscle3.8.31_i86linux64')
parser.add_argument('--raxml-path', default=glob.glob('./packages/standard-RAxML/raxmlHPC-*')[0])

args = parser.parse_args()
args.glsfnames = utils.get_arg_list(args.glsfnames)
args.glslabels = utils.get_arg_list(args.glslabels)
if not os.path.exists(args.muscle_path):
    raise Exception('muscle path %s does not exist' % args.muscle_path)
if not os.path.exists(args.raxml_path):
    raise Exception('raxml path %s does not exist' % args.raxml_path)

plot_gls_gen_tree(args, args.plotdir, args.plotname, args.glsfnames, args.glslabels, title=args.title)
