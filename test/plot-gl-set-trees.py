#!/usr/bin/env python
import copy
import random
import os
import tempfile
import subprocess
import sys
sys.path.insert(1, './python')
import ete3

import utils

gls_labels = ['sim', 'inf']
workdir = '/tmp/' + os.getenv('USER') + '/gls-trees/' + str(random.randint(0, 999999))
aligned_fname = workdir + '/all-aligned.fa'
raxml_label = 'test'
raxml_output_fnames = ['%s/RAxML_%s.%s' % (workdir, fn, raxml_label) for fn in ['parsimonyTree', 'log', 'result', 'info', 'bestTree']]
tree_fname = [fn for fn in raxml_output_fnames if 'result' in fn][0]
utils.prep_dir(workdir)

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print '%s %s' % (utils.color('red', 'run'), cmd_str)
    sys.stdout.flush()
    subprocess.check_call(cmd_str.split())

# read the input germline sets
all_genes, gl_sets = {}, {}
for label in gls_labels:
    gfn = label + '-fiddles-2.fa'
    gl_sets[label] = {gfo['name'] : gfo['seq'] for gfo in utils.read_fastx(gfn)}
    for name, seq in gl_sets[label].items():
        if name not in all_genes:
            all_genes[name] = seq

# write and align an .fa with all genes from any gl set
with tempfile.NamedTemporaryFile() as tmpfile:
    for name, seq in all_genes.items():
        tmpfile.write('>%s\n%s\n' % (name, seq))
    cmd_str = './packages/muscle/muscle3.8.31_i86linux64 -in %s -out %s' % (tmpfile.name, aligned_fname)
    run(cmd_str)

# get a tree for the aligned .fa
cmd_str = './packages/standard-RAxML/raxmlHPC-AVX -mGTRCAT -n%s -s%s -p1 -w %s' % (raxml_label, aligned_fname, workdir)
run(cmd_str)

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
