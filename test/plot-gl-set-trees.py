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
workdir = '/tmp/' + os.getenv('USER') + '/gls-trees/230626' # + str(random.randint(0, 999999))
aligned_fname = workdir + '/all-aligned.fa'
raxml_label = 'test'
raxml_output_fnames = ['%s/RAxML_%s.%s' % (workdir, fn, raxml_label) for fn in ['parsimonyTree', 'log', 'result', 'info', 'bestTree']]
tree_fname = [fn for fn in raxml_output_fnames if 'result' in fn][0]
utils.prep_dir(workdir)

# ----------------------------------------------------------------------------------------
def run(cmd_str):
    print '%s %s' % (utils.color('red', 'run'), cmd_str)
    sys.stdout.flush()
    # subprocess.check_call(cmd_str.split())

# read the input germline sets
all_genes, gl_sets = {}, {}
for label in gls_labels:
    gfn = label + '.fa'
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

# make the text matrix for ete
sorted_all_genes = sorted(all_genes)
matrix_strs = ['#Names\t' + '\t'.join(gls_labels)]
# present = {n : [g in gl_sets[n] for g in sorted_all_genes] for n in gls_labels}
# for igene in range(len(sorted_all_genes)):
#     line = sorted_all_genes[igene] + '\t' + '\t'.join([str(int(present[label][igene])) for label in gls_labels])
#     matrix_strs.append(line)
def color_fcn(gene):
    if gene in gl_sets['sim'] and gene in gl_sets['inf']:
        return '0'  # ok
    elif gene in gl_sets['sim']:
        return '1'  # missing
    elif gene in gl_sets['inf']:
        return '2'  # spurious
    else:
        assert False
for igene in range(len(sorted_all_genes)):
    line = sorted_all_genes[igene] + '\t' + '\t'.join([color_fcn(sorted_all_genes[igene]) for label in gls_labels])
    matrix_strs.append(line)


matrix = '\n'.join(matrix_strs)
# matrix = """
# #Names\tcol1\tcol2
# A\t1\t1
# B\t0\t0
# """
# print matrix
# partis_dir = os.path.dirname(os.path.realpath(__file__)).replace('/test', '')
# sys.path.insert(1, partis_dir + '/packages/baltic')
# import treegenerator
# print treegenerator.get_ascii_tree(treestr.replace('*', '_s').replace('+', '_p_'))

t = ete3.ClusterTree(treestr) #, text_array=matrix)
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
    if n.is_leaf():
        status = getstatus(n.name)
        n.set_style(get_nst(status))
        if status in faces:
            print n.name, status
            n.add_face(copy.deepcopy(faces[status]), column=0)

ts.show_leaf_name = False
ts.mode = 'c'
ts.show_scale = False
t.render('out.png', h=750, tree_style=ts)  # , layout='heatmap'

sys.exit()
os.remove(aligned_fname)
for fn in raxml_output_fnames:
    os.remove(fn)
os.rmdir(workdir)
