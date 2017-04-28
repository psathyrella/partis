#!/usr/bin/env python
# has to be its own script, since ete3 requires its own god damn python version, installed in a separated directory
import argparse
import copy
import random
import os
import tempfile
import subprocess
import sys
sys.path.insert(1, './python')
import ete3

import utils

# ----------------------------------------------------------------------------------------
def plot_gls_gen_tree(plotdir, plotname, glsfnames, glslabels, leg_title=None, title=None):
    assert len(glslabels) == len(set(glslabels))  # no duplicates
    if 'ete3' not in sys.modules:
        import ete3
    ete3 = sys.modules['ete3']

    workdir = '/tmp/' + os.getenv('USER') + '/gls-trees/' + str(random.randint(0, 999999))
    musclepath = './packages/muscle/muscle3.8.31_i86linux64'
    if not os.path.exists(musclepath):
        raise Exception('muscle path %s does not exist' % musclepath)
    sys.exit()
    aligned_fname = workdir + '/all-aligned.fa'
    raxml_label = 'xxx'
    raxml_output_fnames = ['%s/RAxML_%s.%s' % (workdir, fn, raxml_label) for fn in ['parsimonyTree', 'log', 'result', 'info', 'bestTree']]
    tree_fname = [fn for fn in raxml_output_fnames if 'result' in fn][0]
    utils.prep_dir(workdir)

    # read the input germline sets
    all_genes, gl_sets = {}, {}
    for label, fname in zip(glslabels, glsfnames):
        gl_sets[label] = {gfo['name'] : gfo['seq'] for gfo in utils.read_fastx(fname)}
        for name, seq in gl_sets[label].items():
            if name not in all_genes:
                all_genes[name] = seq

    # write and align an .fa with all alleles from any gl set
    with tempfile.NamedTemporaryFile() as tmpfile:
        for name, seq in all_genes.items():
            tmpfile.write('>%s\n%s\n' % (name, seq))
        cmd_str = '%s -in %s -out %s' % (musclepath, tmpfile.name, aligned_fname)
        utils.simplerun(cmd_str)

    # get a tree for the aligned .fa
    cmd_str = './packages/standard-RAxML/raxmlHPC-AVX -mGTRCAT -n%s -s%s -p1 -w %s' % (raxml_label, aligned_fname, workdir)
    utils.simplerun(cmd_str)

    with open(tree_fname) as treefile:
        treestr = treefile.read().strip()

    t = ete3.ClusterTree(treestr)
    ts = ete3.TreeStyle()
    # treestr = "(A:0.7,B:0.7):0.3;"

    scolors = {'ok' : 'DarkSeaGreen', 'missing' : 'IndianRed', 'spurious' : 'IndianRed'}
    faces = {'missing'  : ete3.CircleFace(10, 'white'), 
             'spurious' : ete3.CircleFace(10, 'black')}

    def getstatus(gene):
        if gene in gl_sets['sim'] and gene in gl_sets['inf']:
            return 'ok'
        elif gene in gl_sets['sim']:
            return 'missing'
        elif gene in gl_sets['inf']:
            return 'spurious'
        else:
            return 'xxx'

    def get_nst(status):
        nst = ete3.NodeStyle()
        nst['bgcolor'] = scolors[status]
        return nst

    for n in t.traverse():
        n.dist = 1
        n.img_style['hz_line_width'] = 2
        n.img_style['vt_line_width'] = 2
        if n.is_leaf():
            status = getstatus(n.name)
            n.set_style(get_nst(status))
            if status in faces:
                n.add_face(copy.deepcopy(faces[status]), column=0)

    ts.show_leaf_name = False
    ts.mode = 'c'
    ts.show_scale = False
    outfname = 'out.png' #'../papers/germline-set-generation/prefigs/gls-gen/2.svg'
    t.render(outfname, h=750, tree_style=ts)  # , layout='heatmap'

    sys.exit()
    os.remove(aligned_fname)
    for fn in raxml_output_fnames:
        os.remove(fn)
    os.rmdir(workdir)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--plotdir', required=True)
parser.add_argument('--plotname', required=True)
parser.add_argument('--glsfnames', required=True)
parser.add_argument('--glslabels', required=True)
parser.add_argument('--title')

args = parser.parse_args()
args.glsfnames = utils.get_arg_list(args.glsfnames)
args.glslabels = utils.get_arg_list(args.glslabels)

plot_gls_gen_tree(args.plotdir, args.plotname, args.glsfnames, args.glslabels, title=args.title)
