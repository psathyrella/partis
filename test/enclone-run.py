#!/usr/bin/env python3
from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy
import random
import glob
import time
import colored_traceback.always
import argparse
import subprocess
import sys
import os
import json
from io import open

import utils
import glutils
import paircluster
from clusterpath import ClusterPath

ig_or_tr = 'ig'

# ----------------------------------------------------------------------------------------
def wkdir():
    return '%s/work' % args.outdir

# ----------------------------------------------------------------------------------------
def encifn():
    if args.test_data:
        return '%s/packages/enclone/downloads/datasets/1279049/outs/multi/vdj_b/all_contig_annotations.json' % utils.get_partis_dir()
    else:
        return '%s/enclone/all_contig_annotations.json' % wkdir() #/%s, idlstr(locus))

# ----------------------------------------------------------------------------------------
def translation_fname():
    return '%s/translations.json' % wkdir()

# ----------------------------------------------------------------------------------------
def simfn(locus, lpair=None):
    return paircluster.paired_fn(args.simdir, locus, lpair=lpair, suffix='.yaml')

# ----------------------------------------------------------------------------------------
def getofn(locus):
    return paircluster.paired_fn(args.outdir, locus, actstr='partition', suffix='.yaml')

# ----------------------------------------------------------------------------------------
def gloci():  # if no sw file, we probably only made one light locus (or i guess all input is missing, oh well)
    return [l for l in utils.sub_loci(ig_or_tr) if args.simdir is None or os.path.exists(simfn(l))]

# ----------------------------------------------------------------------------------------
def encofn():
    return '%s/enclone/out.csv' % wkdir()

dummy_jdatas = [
    {
        "barcode": "AAACCTGAGACACGAC-1",
        "cdr3": "CAAWDDSLNVVVF",
        "cdr3_seq": "TGTGCAGCATGGGATGACAGCCTGAATGTTGTGGTATTC",
        "cdr3_start": 361,
        "cdr3_stop": 400,
        "contig_name": "AAACCTGAGACACGAC-1_contig_1",
        "filtered": True,
        "fraction_of_reads_for_this_barcode_provided_as_input_to_assembly": 1.0,
        "high_confidence": True,
        "is_asm_cell": True,
        "is_cell": True,
        "is_gex_cell": True,
        "productive": True,
        "quals": "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
        "read_count": 4054,
        "sequence": "GCTGCGGGTAGAGAAGACAGGACTCAGGACAATCTCCAGCATGGCCTGGTCCCCTCTCTTCCTCACCCTCATCACTCACTGTGCAGGGTCCTGGGCCCAGTCTGTGCTGACTCAGCCACCCTCGGTGTCTGAAGCCCCCAGGCAGAGGGTCACCATCTCCTGTTCTGGAAGCAGCTCCAACATCGGAAATAATGCTGTAAACTGGTACCAGCAGCTCCCAGGAAAGGCTCCCAAACTCCTCATCTATTATGATGATCTGCTGCCCTCAGGGGTCTCTGACCGATTCTCTGGCTCCAAGTCTGGCACCTCAGCCTCCCTGGCCATCAGTGGGCTCCAGTCTGAGGATGAGGCTGATTATTACTGTGCAGCATGGGATGACAGCCTGAATGTTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAGGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCCTCTGAGGAGCTTCAAGCCAACAAGGCCACACTGGTGTGTCTCATAAGTGACTTCTACCCGGGAGCCGTGACAGTGGCCTGGAAGGCAGATAGCAGCCCCGTCAAGGCGGGAGTGGAGACCACCACACCCTCCAAACAAAGCAACAACAAGTACGCGGCCAGCAGCTA",
        "umi_count": 39
    },
    {
        "barcode": "AAACCTGAGACACGAC-1",
        "cdr3": "CARAGAVAASKGYYYYYYGMDVW",
        "cdr3_seq": "TGTGCGAGAGCAGGAGCAGTGGCTGCCTCCAAGGGTTATTACTACTACTACTACGGTATGGACGTCTGG",
        "cdr3_start": 406,
        "cdr3_stop": 475,
        "contig_name": "AAACCTGAGACACGAC-1_contig_2",
        "filtered": True,
        "fraction_of_reads_for_this_barcode_provided_as_input_to_assembly": 1.0,
        "high_confidence": True,
        "is_asm_cell": True,
        "is_cell": True,
        "is_gex_cell": True,
        "productive": True,
        "quals": "]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]",
        "read_count": 4981,
        "sequence": "GGGAGCATCACCCAGCAACCACATCTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTGGAGGATCCTCTTCTTGGTGGCAGCAGCCACAGGAGCCCACTCCCAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCTTCTGGATACACCTTCACCGGCTACTATATGCACTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGATGGATCAACCCTAACAGTGGTGGCACAAACTATGCACAGAAGTTTCAGGGCAGGGTCACCATGACCAGGGACACGTCCATCAGCACAGCCTACATGGAGCTGAGCAGGCTGAGATCTGACGACACGGCCGTGTATTACTGTGCGAGAGCAGGAGCAGTGGCTGCCTCCAAGGGTTATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAGGGAGTGCATCCGCCCCAACCCTTTTCCCCCTCGTCTCCTGTGAGAATTCCCCGTCGGATACGAGCAGCGTG",
        "umi_count": 38
    },
]

leader_seqs = { 'h' : 'GGGAGCATCACCCAGCAACCACATCTGTCCTCTAGAGAATCCCCTGAGAGCTCCGTTCCTCACCATGGACTGGACCTGGAGGATCCTCTTCTTGGTGGCAGCAGCCACAGGAGCCCACTCC',
                'l' : 'GCTGCGGGTAGAGAAGACAGGACTCAGGACAATCTCCAGCATGGCCTGGTCCCCTCTCTTCCTCACCCTCATCACTCACTGTGCAGGGTCCTGGGCC'}
base_jdt = {
        "filtered": True,
        "fraction_of_reads_for_this_barcode_provided_as_input_to_assembly": 1.0,
        "high_confidence": True,
        "is_asm_cell": True,
        "is_cell": True,
        "is_gex_cell": True,
        "productive": True,
        "read_count": 4981,
        "umi_count": 38,
    }

# ----------------------------------------------------------------------------------------
def ctgid(tch):
    assert tch in 'hl'
    return 'contig_%d' % ('hl'.index(tch) + 1)

# ----------------------------------------------------------------------------------------
def get_jdata(sfo, barstr, translations, tch=None):
    jfo = {k : v for k, v in base_jdt.items()}
    jfo['barcode'] = '%s-1' % barstr
    if tch is None:  # data
        jfo['contig_name'] = sfo['name']
        jfo['sequence'] = sfo['seq']
        if barstr not in translations:
            translations[barstr] = []
        translations[barstr].append(sfo['name'])
    else:
        jfo['contig_name'] = '%s_%s' % (jfo['barcode'], ctgid(tch))
        jfo['sequence'] = leader_seqs[tch] + sfo['%s_seq'%tch]
        antn = sfo['%s_antn'%tch]
        # jfo['cdr3_seq'] = utils.get_cdr3_seq(antn, antn['unique_ids'].index(sfo['%s_id'%tch]))
        # jfo['cdr3'] = utils.ltranslate(jfo['cdr3_seq'])
        # jfo['cdr3_start'] = antn['codon_positions']['v'] + 3
        # jfo['cdr3_stop'] = antn['codon_positions']['j'] + 3
        translations[jfo['contig_name']] = sfo['%s_id'%tch]
    jfo['quals'] = len(jfo['sequence']) * ']'
    return jfo

# ----------------------------------------------------------------------------------------
def convert_input():
    ofn = encifn()
    if utils.output_exists(args, ofn, debug=False):
        print('    converted input already there %s' % ofn)
        return
    jdatas, translations = [], {}
    if args.simdir is not None:
        lp_infos = paircluster.read_lpair_output_files(utils.locus_pairs[ig_or_tr], simfn, dont_add_implicit_info=True, debug=True)
        _, antn_lists, _ = paircluster.concat_heavy_chain(utils.locus_pairs[ig_or_tr], lp_infos, dont_deep_copy=True)
        outfos = paircluster.find_seq_pairs(antn_lists)
        potential_names, used_names = None, None
        for sfo in outfos:
            barstr, potential_names, used_names = utils.choose_new_uid(potential_names, used_names, initial_length=16, n_initial_names=len(outfos), available_chars='ACGT', repeat_chars=True, dont_extend=True)
            for tch in 'hl':
                jdatas.append(get_jdata(sfo, barstr, translations, tch=tch))
    elif args.datafname is not None:
        seqfos = utils.read_fastx(args.datafname)
        for sfo in seqfos:
            did = utils.get_droplet_id(sfo['name'], '-_', [0])
            jdatas.append(get_jdata(sfo, did, translations))
    else:
        assert False
    # jdatas = dummy_jdatas
    utils.mkdir(ofn, isfile=True)
    raise Exception('see comment in next line')
    utils.jsdump(ofn, jdatas)  # need to pass , indent=0 but don't want to add atm
    flines = []
    with open(ofn) as jfile:
        for line in jfile:
            if line.strip() in ['{', '[', ']']:
                flines.append([line.rstrip()])
            else:
                flines[-1].append(line.rstrip())
    with open(ofn, 'w') as jfile:
        for fl in flines:
            jfile.write('%s\n' % ' '.join(fl))
    utils.jsdump(translation_fname(), translations)

# ----------------------------------------------------------------------------------------
def run_enclone():
    ofn, cmdfos, n_already_there, n_total = None, [], 0, 1
    ofn = encofn()
    if utils.output_exists(args, ofn, debug=False): # and not args.dry:  # , offset=8):
        print('    enclone output already there %s' % ofn)
        return
    # description of output fields here https://10xgenomics.github.io/enclone/pages/auto/help.parseable.html
    out_columns = ['barcode', 'group_id', 'datasets_cell']  # still don't understand what these other ones ('clonotype_id', 'exact_subclonotype_id') are, but they're not the index of the family
    cmd = '%s/packages/enclone/enclone BCR=%s BUILT_IN SUMMARY NOPRINT POUT=%s PCELL PCOLS=%s' % (utils.get_partis_dir(), os.path.dirname(encifn()), encofn(), ','.join(out_columns))  #  FASTA=_output/enclone-test/out.fa
    cmdfos += [{
        'cmd_str' : cmd,
        'outfname' : ofn,
        'logdir' : wkdir(),
        'workdir' : wkdir(),
    }]
    start = time.time()
    utils.run_scan_cmds(args, cmdfos, 'enclone.log', n_total, n_already_there, ofn, dbstr='enclone run')
    print('enclone time: %.1f' % (time.time() - start))

# ----------------------------------------------------------------------------------------
def convert_output():
    if utils.all_outputs_exist(args, [getofn(l) for l in gloci()], debug=False):
        print('    converted ouputs already there %s' % ', '.join(getofn(l) for l in gloci()))
        return
    airr_fn = '%s/airr-partition.tsv' % wkdir()
    utils.simplerun('sed -e \"s/group_id/clone_id/\" -e \"s/barcode/sequence_id/\" %s >%s' % (encofn(), airr_fn), shell=True)
    utils.simplerun('sed -i \"s/,/\t/g\" %s' % airr_fn, shell=True)
    tptfn = '%s/tmp-partitions.yaml' % wkdir()
    utils.simplerun('%s/bin/parse-output.py %s %s --airr-input' % (utils.get_partis_dir(), airr_fn, tptfn))
    _, _, barcode_cpath = utils.read_output(tptfn)
    with open(translation_fname()) as tfile:
        translations = json.load(tfile)
    if args.simdir is None:  # data
        seqfos = utils.read_fastx(args.datafname)
        input_uids = set(s['name'] for s in seqfos)
        for cluster in barcode_cpath.best():
            for barstr in cluster:
                for uid in translations[barstr.replace('-1', '')]:
                    input_uids.remove(uid)
        if len(input_uids) > 0:
            print('  %s missing %d / %d input uids: %s' % (utils.wrnstr(), len(input_uids), len(seqfos), ' '.join(input_uids)))
        else:
            print('  found all input uids')
    else:
        cpaths = {l : ClusterPath(partition=[]) for l in gloci()}
        for bclust in barcode_cpath.best():
            ch_clusts = [[] for _ in 'hl']
            loci = {c : None for c in 'hl'}
            for bcode in bclust:
                for tch, cclust in zip('hl', ch_clusts):
                    encid = '%s_%s' % (bcode, ctgid(tch))
                    uid = translations[encid]
                    cclust.append(uid)
                    if loci[tch] is None:
                        loci[tch] = uid.split('-')[-1]
                        assert loci[tch] in gloci()  # will need to update to run on data that doesn't have the locus in the uid
                        cpaths[loci[tch]].partitions[0].append(cclust)
        for locus in gloci():
            _, _, true_cpath = utils.read_output(simfn(locus), skip_annotations=True)
            plines = cpaths[locus].get_partition_lines(true_partition=true_cpath.best(), calc_missing_values='best', fail_frac=0.10)
            print('    writing partition to %s' % getofn(locus))
            utils.write_annotations(getofn(locus), {}, [], utils.annotation_headers, partition_lines=plines)

# ----------------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--simdir')  # set either --simdir (partis simulation dir) or --datafname (10x fasta data file)
parser.add_argument('--datafname')
parser.add_argument('--outdir', required=True)
parser.add_argument('--test-data', action='store_true', help='run on enclone test data')
parser.add_argument('--overwrite', action='store_true')
parser.add_argument('--dry', action='store_true')
parser.add_argument('--n-max-procs', type=int, help='NOT USED')
parser.add_argument('--n-procs', type=int, help='NOT USED')
args = parser.parse_args()
assert [args.simdir, args.datafname].count(None) == 1

if not args.test_data:
    convert_input()
run_enclone()
convert_output()
