#!/usr/bin/env python
import random
import os
import tempfile
import subprocess
import sys
sys.path.insert(1, './python')
from ete3 import ClusterTree

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
present = {n : [g in gl_sets[n] for g in sorted_all_genes] for n in gls_labels}
matrix_strs = ['#Names\t' + '\t'.join(gls_labels)]
for igene in range(len(sorted_all_genes)):
    line = sorted_all_genes[igene] + '\t' + '\t'.join([str(int(present[label][igene])) for label in gls_labels])
    matrix_strs.append(line)

matrix = '\n'.join(matrix_strs)
# treestr = "(A:0.7,B:0.7):0.3;"
# matrix = """
# #Names\tcol1\tcol2
# A\t1\t1
# B\t0\t0
# """
print matrix

t = ClusterTree(treestr, text_array=matrix)
t.show("heatmap")

sys.exit()
os.remove(aligned_fname)
for fn in raxml_output_fnames:
    os.remove(fn)
os.rmdir(workdir)
