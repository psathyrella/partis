from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import platform
import resource
import psutil
import numpy
import tempfile
import string
import time
import sys
import os
import random
import itertools
import ast
import math
import glob
from collections import Counter
from collections import OrderedDict
from collections import defaultdict
import csv
import subprocess
import multiprocessing
import copy
import traceback
import json
import types
import collections
import operator
import yaml
import six
import hashlib
from pathlib import Path
from io import open
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

existing_output_actions = []  # set in bin/partis (so it's next to the action dict)

all_ptn_plot_cfg = ['shm-vs-size', 'cluster-bubble', 'mut-bubble', 'diversity', 'sizes', 'trees', 'subtree-purity', 'mds', 'laplacian-spectra', 'sfs', 'vrc01-muts', 'no-slug']
default_ptn_plot_cfg = ['shm-vs-size', 'diversity', 'cluster-bubble', 'sizes', 'trees']

dummy_str = 'x-dummy-x'

def csv_wmode(mode='w'):
    if sys.version_info.major < 3:  # 2.7 csv module doesn't support unicode, this is the hackey fix
        return mode + 'b'
    else:
        return mode

# ----------------------------------------------------------------------------------------
def get_partis_dir():
    return str(Path(__file__).resolve().parent.parent)

# ----------------------------------------------------------------------------------------
def fsdir():
    fsdir = '/fh/fast/matsen_e'
    if not os.path.exists(fsdir):
        fsdir = '/tmp'
    if os.getenv('USER') is not None:
        fsdir += '/' + os.getenv('USER')
    return fsdir

# ----------------------------------------------------------------------------------------
def choose_random_subdir(dirname, make_dir=False):
    subname = str(random.randint(0, 999999))
    while os.path.exists(dirname + '/' + subname):
        subname = str(random.randint(0, 999999))
    if make_dir:
        prep_dir(dirname + '/' + subname)
    return dirname + '/' + subname

# ----------------------------------------------------------------------------------------
def timeprinter(fcn):
    def wrapper(*args, **kwargs):
        start = time.time()
        # print fcn.__name__,
        fcn(*args, **kwargs)
        print('    %s: (%.1f sec)' % (fcn.__name__, time.time()-start))
    return wrapper

# ----------------------------------------------------------------------------------------
# get string to print with a specific width
# NOTE there are many many many places where i could use this
# see it's short for "width format"...
def wfmt(val, width, fmt='s', jfmt=''):  # return string with <val> formatted to <width> and format <fmt> (the necessary parentheses really make stuff hard to read without this fcn)
    return ('%' + jfmt + str(width) + fmt) % val

# ----------------------------------------------------------------------------------------
parameter_type_choices = ('multi-hmm', 'hmm', 'sw')  # NOTE this order determines default priority, i.e. if not set on the command line we choose the first one in this order that exists
default_parameter_type = 'hmm'  # not 'default' in the sense of we always use it if user doesn't set something, but default in terms of we want to set something if none of them exist (especially when caching parameters)
def get_parameter_type(args, paramdir):
    if paramdir is None:
        return None
    ptype = default_parameter_type
    if args.parameter_type is not None:  # if it was set on the command line, use that
        ptype = args.parameter_type
    else:  # otherwise take the first one that actually exists
        for tpt in parameter_type_choices:
            if os.path.exists('%s/%s' % (paramdir, tpt)):
                ptype = tpt
                break
    return ptype
def parameter_type_subdir(args, paramdir):  # NOTE this is similar to a bit in getpdir() in the paired locus stuff in bin/partis.py
    if paramdir is None:
        return None
    return '%s/%s' % (paramdir, get_parameter_type(args, paramdir))

# ----------------------------------------------------------------------------------------
# putting these up here so glutils import doesn't fail... I think I should be able to do it another way, though
regions = ['v', 'd', 'j']
constant_regions = ['c', 'm', 'g', 'a', 'd', 'e']  # NOTE d is in here, which is stupid but necessary, so use is_constant_gene()
chains = ['h', 'l']
loci = collections.OrderedDict((
    ('igh', 'vdj'),
    ('igk', 'vj'),
    ('igl', 'vj'),
    ('tra', 'vj'),
    ('trb', 'vdj'),
    ('trg', 'vj'),
    ('trd', 'vdj'),
))
isotypes = ['m', 'g', 'k', 'l']
locus_pairs = {'ig' : [['igh', 'igk'], ['igh', 'igl']],
               'tr' : [['trb', 'tra'], ['trd', 'trg']]}
# ----------------------------------------------------------------------------------------
def locstr(l):
    cdict = {'igh' : 'blue', 'igk' : 'purple', 'igl' : 'green'}
    return color(cdict.get(l, 'red'), l.replace('ig', ''))

def lpstr(lp):
    return '%s+%s' % (locstr(lp[0]), locstr(lp[1]))

def sub_loci(ig_or_tr):  # ok i probably should have just split <loci> by ig/tr, but too late now
    return [l for l in loci if ig_or_tr in l]

def heavy_locus(ig_or_tr):
    return get_single_entry([l for l in sub_loci(ig_or_tr) if has_d_gene(l)])

def light_loci(ig_or_tr):
    return [l for l in sub_loci(ig_or_tr) if not has_d_gene(l)]

def getlpair(ltmp):  # return heavy/light locus pair for light chain locus <ltmp>
    if has_d_gene(ltmp): raise Exception('only makes sense for light chain, but got %s' % ltmp)
    return get_single_entry([lp for lplist in locus_pairs.values() for lp in lplist if ltmp in lp])

def getregions(locus):  # for clarity, don't use the <loci> dictionary directly to access its .values()
    return list(loci[locus])  # doesn't really need to be a list, but it's more clearly analagous to regions & co that way

def has_d_gene(locus):  # for clarity, don't use the <loci> dictionary directly to access its .values()
    return 'd' in loci[locus]

def samechain(l1, l2):  # true if l1 and l2 are either both heavy, or both light
    if has_d_gene(l1) and has_d_gene(l2):
        return True
    elif has_d_gene(l1) or has_d_gene(l2):
        return False
    else:
        return True

def get_boundaries(locus):  # NOTE almost everything still uses the various static boundaries variables, rather than calling this. It may or may not be more sensible to switch to this eventually
    rlist = getregions(locus)
    return [rlist[i] + rlist[i+1] for i in range(len(rlist)-1)]

def region_pairs(locus):
    return [{'left' : bound[0], 'right' : bound[1]} for bound in get_boundaries(locus)]

#----------------------------------------------------------------------------------------
# NOTE I also have an eps defined in hmmwriter. Simplicity is the hobgoblin of... no, wait, that's just plain ol' stupid to have two <eps>s defined
eps = 1.0e-10  # if things that should be 1.0 are this close to 1.0, blithely keep on keepin on. kinda arbitrary, but works for the moment
def is_normed(probs, this_eps=eps, total=1.):
    if hasattr(probs, 'keys'):  # if it's a dict, call yourself with a list of the dict's values
        return is_normed([val for val in probs.values()], this_eps=this_eps)
    elif hasattr(probs, '__iter__'):  # if it's a list call yourself with their sum
        return is_normed(sum(probs), this_eps=this_eps)
    else:  # and if it's a float actually do what you're supposed to do
        return math.fabs(probs - total) < this_eps

# ----------------------------------------------------------------------------------------
# return copy of vlist with floats normalized to 1 (I kind of thought there was already a fcn that did this? maybe not)
def normalize(vlist):
    tot = sum(vlist)
    return [v / float(tot) for v in vlist]

# ----------------------------------------------------------------------------------------
def pass_fcn(val):  # dummy function for conversions (see beloww)
    return val

# ----------------------------------------------------------------------------------------
# NOTE returns *all* lines as list by default, i.e. will read entire file at once before returning
def csvlines(fn, n_max_queries=None, delimiter=','):  # for search: csv_lines()
    with open(fn) as cfile:
        reader = csv.DictReader(filter(lambda l: l[0]!='#', cfile), delimiter=delimiter)
        if n_max_queries is None:
            return list(reader)
        else:  # ick, this seems overly verbose
            rlist = []
            for tline in reader:
                rlist.append(tline)
                if len(rlist) >= n_max_queries:
                    break
            return rlist

# ----------------------------------------------------------------------------------------
# make lists from args that are passed as strings of colon-separated values
def get_arg_list(arg, intify=False, intify_with_ranges=False, floatify=False, boolify=False, translation=None, list_of_lists=False, key_val_pairs=False, key_list=None, choices=None, forbid_duplicates=False, default_vals=None):
    if arg is None:
        return None

    convert_fcn = pass_fcn
    if intify:
        convert_fcn = int
    elif intify_with_ranges:  # allow both plain integers and ranges (specied with a dash), e.g. 0:3-6:50 --> 0:3:4:5:50
        def iwr_fcn(vstr):
            if '-' in vstr:
                istart, istop = [int(v) for v in vstr.split('-')]
                return list(range(istart, istop))  # isn't really right, since we still need to flatten this sublist
            else:
                return [int(vstr)]  # make single values list of length one for easier flattening below
        convert_fcn = iwr_fcn
    elif floatify:
        convert_fcn = float
    elif boolify:
        convert_fcn = bool

    arglist = arg.strip().split(':')  # to allow ids with minus signs, you can add a space (if you don't use --name=val), which you then have to strip() off
    if key_list:  # like key_val_pairs, but the arg just has a list, and here we supply the keys (obviously assumes order is the same)
        assert not key_val_pairs
        if len(key_list) != len(arglist):
            raise Exception('key_list %s different length to arglist %s: %d vs %d' % (key_list, arglist, len(key_list), len(arglist)))
        arglist = ['%s,%s' % (k, v) for k, v in zip(key_list, arglist)]
        key_val_pairs = True
    if list_of_lists or key_val_pairs:
        arglist = [substr.split(',') for substr in arglist]
        if list_of_lists:
            arglist = [[convert_fcn(p) for p in sublist] for sublist in arglist]
    else:
        arglist = [convert_fcn(x) for x in arglist]

    if intify_with_ranges:
        arglist = [v for vl in arglist for v in vl]

    if translation is not None:
        for ia in range(len(arglist)):
            if arglist[ia] in translation:
                arglist[ia] = translation[arglist[ia]]

    if key_val_pairs:
        if any(len(l) != 2 for l in arglist):
            raise Exception('need arg list element of length 2 to convert to key/val pair but got: %s' % str([l for l in arglist if len(l) != 2]))
        vdict = {}
        for tk, tval in arglist:
            if default_vals is not None and tval == '':  # if <default_vals> is set, empty strings get replaced with value from <default_vals>
                tval = default_vals[tk]
            if tk in vdict:
                vdict[tk] += ':'+convert_fcn(tval)
            else:
                vdict[tk] = convert_fcn(tval)
        arglist = vdict

    if choices is not None:  # note that if <key_val_pairs> is set, this (and <forbid_duplicates) is just checking the keys, not the vals
        for arg in arglist:
            if arg not in choices:
                raise Exception('unexpected argument \'%s\' (choices: %s)' % (str(arg), [str(c) for c in choices]))

    if forbid_duplicates:
        if any(arglist.count(a) > 1 for a in arglist):
            raise Exception('duplicate values for arg: %s' % arg)

    return arglist

# ----------------------------------------------------------------------------------------
def add_to_scan_cmd(args, vname, vstr, cmd, replacefo=None):
    if vstr == 'None':  # we want 'None' to show up in the scan dir structure, but it's nicer to just not set the command line arg (assuming the script we're running has None as the default!)
        return cmd
    if hasattr(args, 'bool_args') and vname in args.bool_args:  # in cf-tree-metrics this was handled for e.g. context dependence in bcr-phylo-run
        if vstr == '0':
            return cmd
        elif vstr == '1':
            vstr = ''  # adds extra space, darn it
        else:
            assert False
    if replacefo is not None and vname in replacefo:
        if replacefo[vname][vstr] is None:
            return cmd
        vname, vstr = replacefo[vname][vstr]
    return cmd + ' --%s%s' % (vname, '' if vstr=='' else ' '+vstr)

# ----------------------------------------------------------------------------------------
def svoutdir(args, varnames, valstrs, svtype):
    assert len(varnames) == len(valstrs)
    outdir = [args.base_outdir, args.label]
    if hasattr(args, 'version'):
        outdir.append(args.version)
    for vn, vstr in zip(varnames, valstrs):
        if vn not in args.scan_vars[svtype]:
            continue
        if len(varnames)==1 and vstr is None:  # used force_val when getting var info
            continue
        outdir.append('%s-%s' % (vn, vstr))
    return '/'.join(outdir)

# ----------------------------------------------------------------------------------------
def get_scanvar_arg_lists(args):
    # ----------------------------------------------------------------------------------------
    def set_arg_list(aname):  # NOTE we *don't* want to set intify or floatify here since the dir structure stuff is too hard if we don't have strings; conversions happen only for plotting axes
        attr_name = aname.replace('-', '_') + '_list'
        if not hasattr(args, attr_name):
            raise Exception('args does not have attribute \'%s\', i.e. you probably need to add an arg called --%s-list' % (attr_name, aname))
        arglist = get_arg_list(getattr(args, attr_name), list_of_lists=aname in args.str_list_vars, #, intify=aname in args.svartypes['int'], floatify=aname in args.svartypes['float'],
                               forbid_duplicates=args.zip_vars is None or aname not in args.zip_vars)
        setattr(args, attr_name, arglist)  # if we're zipping the var, we have to allow duplicates, but then check for them again after we've done combos in get_var_info()
    # ----------------------------------------------------------------------------------------
    for aname in set([v for vlist in args.scan_vars.values() for v in vlist]):
        if aname == 'seed':
            continue
        set_arg_list(aname)

# ----------------------------------------------------------------------------------------
def sargval(args, sv):  # ick this name sucks
    def dkey(sv):
        return sv.replace('-', '_') + '_list'
    if sv == 'seed':
        riter = list(range(args.n_replicates)) if args.iseeds is None else args.iseeds
        return [args.random_seed + i for i in riter]
    else:
        return args.__dict__[dkey(sv)]

# ----------------------------------------------------------------------------------------
# return value in <vlist> corresponding to <vname> using index of <vname> in <varnames>
def vlval(args, vlist, varnames, vname):  # ok this name also sucks, but they're doing complicated things while also needing really short names...
    if vname in varnames:
        return vlist[varnames.index(vname)]
    else:
        assert len(sargval(args, vname))  # um, I think?
        return sargval(args, vname)[0]

# ----------------------------------------------------------------------------------------
# handles the complicated machinations of scanning over a variety of (mostly simulation) variables, taking as input the command line args, and returning several lists e.g. with these separated into those with/without multiple values
# <args.str_list_vars>: vars that are lists of lists of strs, rather than simple lists of strs
def get_var_info(args, scan_vars, action=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def handle_var(svar, val_lists, valstrs, force_val=False):  # force_val: this is the last var, and we haven't yet gotten one with more than one value, so we need to put one in or there won't be any, and things will crash later on
        # convert_fcn = (lambda vlist: ':'.join(str(v) for v in vlist)) if svar in args.str_list_vars else str  # old version, leaving here cause i'm chicken that i messed something up
        convert_fcn = str
        if svar in args.str_list_vars:
            def convert_fcn(vlist):
                subfcn = str
                if hasattr(args, 'recurse_replace_vars') and svar in args.recurse_replace_vars:
                    subfcn = lambda s: s.replace('M', ',').replace('L', ':')  # if you need , or : in subprocess args, put in M and L and they'll get replaced here
                return ':'.join(subfcn(v) for v in vlist)
        sargv = sargval(args, svar)
        if sargv is None:  # no default value, and it wasn't set on the command line
            pass
        elif len(sargv) > 1 or (svar == 'seed' and args.iseeds is not None):  # if --iseeds is set, then we know there must be more than one replicate, but/and we also know the fcn will only be returning one of 'em
            varnames.append(svar)
            val_lists = [vlist + [sv] for vlist in val_lists for sv in sargv]
            valstrs = [vlist + [convert_fcn(sv)] for vlist in valstrs for sv in sargv]
            if debug:
                print('    %30s  %4s    %s    %s' % (svar, 'y' if svar in args.str_list_vars else '-', val_lists, valstrs))
        else:
            if force_val:
                if debug:
                    print('    force_val: setting %s val str to None' % svar)
                varnames.append(svar)  # this adds it to the list of <varnames> (so the calling script has something to loop over) but sets the val str to None (so it knows not to e.g. include it in the output dir name)
                val_lists = [[None]]
                valstrs = [[None]]
            if hasattr(args, 'bool_args') and svar in args.bool_args:  # duplicates code in add_to_scan_cmd()
                assert sargv[0] in ['0', '1']
                if sargv[0] == '1':
                    base_args.append('--%s'%svar)
            else:
                base_args.append('--%s %s' % (svar, convert_fcn(sargv[0])))
        return val_lists, valstrs

    # ----------------------------------------------------------------------------------------
    base_args = []  # list of command line args/values that only have one value (but support the possibility of more than one), i.e. that don't need to be looped over and don't need a subdir in the output path (e.g. ['--carry-cap 1000', '--obs-times 100:150:200:250:300', '--metric-for-target-distance aa', '--selection-strength 1.0', '--leaf-sampling-scheme uniform-random', '--target-count 1'])
    varnames = []  # names of variables that currently have more than one value that we want to scan over (e.g. ['n-sim-seqs-per-gen', 'seed'])
    val_lists, valstrs = [[]], [[]]  # these are both lists in which each entry is a combination of each of the vars in <varnames> that we want to run; in the former they're all lists, whereas in the latter they're strings, i.e. colon-separated lists (e.g. [[[10, 20], 0], [[10, 20], 1], [[10, 20], 2], [[50, 100], 0], [[50, 100], 1], [[50, 100], 2]] and [['10:20', '0'], ['10:20', '1'], ['10:20', '2'], ['50:100', '0'], ['50:100', '1'], ['50:100', '2']])
    if debug:
        print('                             name    list   val_lists              valstrs')
    non_null_vars = [v for v in scan_vars if sargval(args, v) is not None]
    for iv, svar in enumerate(non_null_vars):
        val_lists, valstrs = handle_var(svar, val_lists, valstrs, force_val=iv==len(non_null_vars)-1 and len(varnames)==0)  # val_lists and valstrs get updated each time through

    if args.zip_vars is not None:
        tzvars = [v for v in varnames if v in args.zip_vars]  # want order to come from <varnames> so it matches below
        if any(v not in varnames for v in args.zip_vars):
            print('  %s zip var[s] %s not in varnames (%s) for %s action' % (wrnstr(), ', '.join(v for v in args.zip_vars if v not in varnames), ', '.join(varnames), action if action is not None else 'this'))
        zipvals = list(zip(*[sargval(args, zv) for zv in tzvars]))
        assert len(set(len(vlist) for vlist in zipvals)) == 1  # doesn't make sense unless you provide a corresponding value for each
        zval_lists, zvalstrs = [], []  # new ones, only containing zipped values
        for vlist, vstrlist in zip(val_lists, valstrs):
            zvals = tuple([vlval(args, vlist, varnames, zv) for zv in tzvars])  # values for this combo of the vars we want to zip
            if zvals in zipvals and vlist not in zval_lists:  # second clause is to avoid duplicates (duh), which we get because when we're zipping vars we have to allow duplicate vals in each zip'd vars arg list, and then (above) we make combos including all those duplicate combos
                zval_lists.append(vlist)
                zvalstrs.append(vstrlist)
        val_lists = zval_lists
        valstrs = zvalstrs
        print('    %s: zipped values for %s:  lists %s    strs %s' % (color('blue', action), ' '.join(tzvars), val_lists, valstrs))

    if any(valstrs.count(vstrs) > 1 for vstrs in valstrs):
        raise Exception('duplicate combinations for %s' % ' '.join(':'.join(vstr) for vstr in valstrs if valstrs.count(vstr) > 1))

    if debug:
        print('    base args: %s' % ' '.join(base_args))
    if len(varnames) == 0:  # shouldn't need this any more, but don't want to remove it yet
        print('  %s zero length varnames in get_var_info() (probably because no variables have multiple values)' % wrnstr())
    return base_args, varnames, val_lists, valstrs

# ----------------------------------------------------------------------------------------
def run_scan_cmds(args, cmdfos, logfname, n_total, n_already_there, single_ofn, example_existing_ofn=None, n_missing_input=0, single_ifn=None, shell=False, dbstr=None):
    if n_missing_input > 0:
        print('      %2d / %d %s missing input%s' % (n_missing_input, n_total, '' if dbstr is None else dbstr+' ', '' if single_ifn is None else ', e.g. %s'%single_ifn))
    if n_already_there > 0:
        print('      %2d / %d %sjobs skipped (outputs exist, e.g. %s)' % (n_already_there, n_total, '' if dbstr is None else dbstr+' ', single_ofn if single_ofn is not None else example_existing_ofn))
    if len(cmdfos) > 0:
        print('      %s %d %sjobs' % ('--dry: would start' if args.dry else 'starting', len(cmdfos), '' if dbstr is None else dbstr+' '))
        if args.dry:
            if hasattr(args, 'print_all') and args.print_all:
                print('  would run:')
                for cfo in cmdfos:
                    print('    %s' % cfo['cmd_str'])
            else:
                print('  first command: %s' % cmdfos[0]['cmd_str'])
        else:
            run_cmds(cmdfos, debug='write:%s'%logfname, n_max_procs=args.n_max_procs, allow_failure=True, shell=shell)

# ----------------------------------------------------------------------------------------
def add_scanvar_args(parser, script_base, all_perf_metrics, default_plot_metric='partition'):
    parser.add_argument('--base-outdir', default='%s/partis/%s'%(os.getenv('fs', default=os.getenv('HOME')), script_base))
    parser.add_argument('--version', default='v0')
    parser.add_argument('--label', default='test')
    parser.add_argument('--dry', action='store_true')
    parser.add_argument('--print-all', action='store_true')
    parser.add_argument('--paired-loci', action='store_true', help='same as partis arg')
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--test', action='store_true', help='don\'t parallelize \'plot\' action')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--simu-extra-args', help='extra args to pass for simulation step')
    parser.add_argument('--inference-extra-args', help='extra args to pass to inference steps')
    parser.add_argument('--plot-metrics', default=default_plot_metric, help='these can be either methods (paired-loci) or metrics (tree metrics), but they\'re called metrics in scanplot so that\'s what they have to be called everywhere')
    parser.add_argument('--perf-metrics', default=':'.join(all_perf_metrics), help='performance metrics (i.e. usually y axis) that we want to plot vs the scan vars. Only necessary if --plot-metrics are actually methods (as in cf-paired-loci.py).')
    parser.add_argument('--zip-vars', help='colon-separated list of variables for which to pair up values sequentially, rather than doing all combinations')
    parser.add_argument('--final-plot-xvar', help='variable to put on the x axis of the final comparison plots')
    parser.add_argument('--legend-var', help='non-default "component" variable (e.g. obs-frac) to use to label different lines in the legend')
    parser.add_argument('--x-legend-var', help='derived variable with which to label the x axis (e.g. mfreq [shm %] when --final-plot-x-var is scratch-mute-freq)')
    parser.add_argument('--pvks-to-plot', help='only plot these line/legend values when combining plots')
    parser.add_argument('--use-val-cfgs', action='store_true', help='use plotting.val_cfgs dict (we can\'t always use it)')
    parser.add_argument('--plot-metric-extra-strs', help='extra strs for each metric in --plot-metrics (i.e. corresponding to what --extra-plotstr was set to during get-tree-metrics for that metric)')
    parser.add_argument('--dont-plot-extra-strs', action='store_true', help='while we still use the strings in --plot-metric-extra-strs to find the right dir to get the plot info from, we don\'t actually put the str in the plot (i.e. final plot versions where we don\'t want to see which dtr version it is)')
    parser.add_argument('--combo-extra-str', help='extra label for combine-plots action i.e. write to combined-%s/ subdir instead of combined/')
    parser.add_argument('--make-hist-plots', action='store_true')
    parser.add_argument('--empty-bin-range', help='remove empty bins only outside this range (will kick a warning if this ignores any non-empty bins)')
    parser.add_argument('--workdir')  # default set below
    parser.add_argument('--dataset-in-list', help='list of samples from --data-in-cfg')
    parser.add_argument('--data-in-cfg', help='instead of running simulation to get input for subsequent actions, use input files specified in this yaml cfg (the actual samples to use are specified with --dataset-in-list)')
    parser.add_argument('--n-replicates', default=1, type=int)
    parser.add_argument('--iseeds', help='if set, only run these replicate indices (i.e. these corresponds to the increment *above* the random seed)')
    parser.add_argument('--n-max-procs', type=int, help='Max number of *child* procs (see --n-sub-procs). Default (None) results in no limit.')
    parser.add_argument('--n-sub-procs', type=int, default=1, help='Max number of *grandchild* procs (see --n-max-procs)')
    parser.add_argument('--random-seed', default=0, type=int, help='note that if --n-replicates is greater than 1, this is only the random seed of the first replicate')
    parser.add_argument('--single-light-locus')

# ----------------------------------------------------------------------------------------
def process_scanvar_args(args, after_actions, plot_actions, all_perf_metrics, data_actions=None):
    for act in after_actions + plot_actions:  # actions that happen after simu need to have all the simu scan vars included in them
        if data_actions is not None and act in data_actions:
            continue
        if act not in args.scan_vars:
            args.scan_vars[act] = []
        args.scan_vars[act] = args.scan_vars['simu'] + args.scan_vars[act]

    args.actions = get_arg_list(args.actions, choices=['simu', 'plot', 'combine-plots', 'data'] + after_actions + plot_actions)
    args.plot_metrics = get_arg_list(args.plot_metrics)
    args.zip_vars = get_arg_list(args.zip_vars)
    args.plot_metric_extra_strs = get_arg_list(args.plot_metric_extra_strs)
    if args.plot_metric_extra_strs is None:
        args.plot_metric_extra_strs = ['' for _ in args.plot_metrics]
    if len(args.plot_metrics) != len(args.plot_metric_extra_strs):
        raise Exception('--plot-metrics %d not same length as --plot-metric-extra-strs %d' % (len(args.plot_metrics), len(args.plot_metric_extra_strs)))
    args.pvks_to_plot = get_arg_list(args.pvks_to_plot)
    args.perf_metrics = get_arg_list(args.perf_metrics, choices=all_perf_metrics)
    args.iseeds = get_arg_list(args.iseeds, intify=True)
    args.empty_bin_range = get_arg_list(args.empty_bin_range, floatify=True)

    if args.dataset_in_list is not None and args.data_in_cfg is None:
        raise Exception('have to set --data-in-cfg if --dataset-in-list is set')
    if args.data_in_cfg is not None:
        assert args.n_replicates == 1
        with open(args.data_in_cfg) as dfile:
            args.data_in_cfg = yaml.load(dfile, Loader=yaml.CLoader)
        if 'simu' in args.actions and not args.dry:
            print('  %s turning on --dry, since the links still get made, and not sure what happens if it isn\'t set' % wrnstr())
            args.dry = True

    get_scanvar_arg_lists(args)
    if args.final_plot_xvar is None:  # set default value based on scan vars
        base_args, varnames, _, valstrs = get_var_info(args, args.scan_vars['simu'], action='simu')
        svars = [v for v in varnames if v != 'seed']
        args.final_plot_xvar = svars[0] if len(svars) > 0 else 'seed'  # if we're not scanning over any vars, i'm not sure what we should use
    if args.dataset_in_list is not None and len(args.dataset_in_list) < 2:
        raise Exception('at least for now you need to specify at least two samples')

# ----------------------------------------------------------------------------------------
def make_scanvar_data_links(args, varnames, vstrs, action, outfn, ofn_fcn, odir_fcn):
    if varnames != ['dataset-in']:  # at least for now
        raise Exception('at least for now, if --data-in-cfg/--dataset-in-list are set, \'dataset-in\' must be a scan var (and be the only one), but got %s' % (varnames))
    sample = vstrs[0]
    if sample not in args.data_in_cfg:
        raise Exception('sample \'%s\' not found in meta info with choices: %s' % (sample, list(args.data_in_cfg.keys())))
    inpath = args.data_in_cfg[sample]['infname']
    if inpath[0] != '/':
        inpath = os.path.realpath(inpath)
    if args.paired_loci and not os.path.isdir(inpath) and os.path.exists(inpath):  # if the data isn't paired (but we're running paired loci), we need to make the paired loci structure UPDATE maybe don't need this any more?
        mkdir(inpath, isfile=True)
        link_name = ofn_fcn(args, varnames, vstrs, action, locus=args.data_in_cfg[sample]['locus'])
        assert args.data_in_cfg[sample]['locus'] == 'igh'  # would need to update things if not
        write_empty_annotations(ofn_fcn(args, varnames, vstrs, action, locus='igk'), 'igk')  # make a fake light file (ick)
    else:  # NOTE haven't run this, it may need testing
        link_name = '%s%s' % (odir_fcn(args, varnames, vstrs, 'simu'), '' if args.paired_loci else '/simu.yaml')
    if not os.path.exists(outfn):
        makelink(os.path.dirname(link_name), inpath, link_name, debug=True)  # , dryrun=True

# ----------------------------------------------------------------------------------------
def add_lists(list_a, list_b):  # add two lists together, except if one is None treat it as if it was zero length (allows to maintain the convention that command line arg variables are None if unset, while still keeping things succinct)
    if list_b is None:
        return copy.deepcopy(list_a)
    elif list_a is None:
        return copy.deepcopy(list_b)
    else:
        return list_a + list_b

# ----------------------------------------------------------------------------------------
def get_single_entry(tl, errmsg=None):  # adding this very late, so there's a lot of places it could be used
    if len(tl) != 1:
        raise Exception('length must be 1 in get_single_entry(), but got %d (%s)%s' % (len(tl), tl, '' if errmsg is None else ': %s'%errmsg))
    return tl[0]

# ----------------------------------------------------------------------------------------
# igh lengths (from code above):
#        mean   min  max
# v     296.7   293  305 (from germline genes)
# d      25.7    11   37
# j      53.6    48   62
# cdr3   52.2    23   96  (from test/ output)
# all   378.2   348  420
typical_lengths = {
    'v' : 297,
    'd' : 26,
    'j' : 54,
    'cdr3' : 52,
    'all' : 378,
}
# ----------------------------------------------------------------------------------------
# values used when simulating from scratch (mutation stuff is controlled by command line args, e.g. --scratch-mute-freq)
scratch_mean_erosion_lengths = {'igh' : {'v_3p' : 1.3, 'd_5p' : 5.6, 'd_3p' : 4.7, 'j_5p' : 5.1},
                                'igk' : {'v_3p' : 2.8, 'd_5p' : 1.0, 'd_3p' : 0.0, 'j_5p' : 1.5},
                                'igl' : {'v_3p' : 2.8, 'd_5p' : 1.0, 'd_3p' : 0.0, 'j_5p' : 1.8}}
scratch_mean_insertion_lengths = {l : {'vd' : 5.9 if has_d_gene(l) else 0.,
                                       'dj' : 5.1 if has_d_gene(l) else 1.6}
                                  for l in loci}

real_erosions = ['v_3p', 'd_5p', 'd_3p', 'j_5p']
# NOTE since we now handle v_5p and j_3p deletions by padding with Ns, the hmm does *not* allow actual v_5p and j_3p deletions.
# This means that while we write parameters for v_5p and j_3p deletions to the parameter dir, these are *not* used in making the
# hmm yamels -- which is what we want, because we want to be able to read in short data reads but make full-length simulation.
effective_erosions = ['v_5p', 'j_3p']
all_erosions = real_erosions + effective_erosions
boundaries = ['vd', 'dj']  # NOTE needs to be updated to be consistent with erosion names
effective_boundaries = ['fv', 'jf']
all_boundaries = boundaries + effective_boundaries
nukes = ['A', 'C', 'G', 'T']
ambig_base = 'N'  # this is the only ambiguous base that we allow/use internally
all_ambiguous_bases = list('NRYKMSWBDHV')  # but the other ones are sometimes handled when parsing input, often just by turning them into Ns
ambiguous_amino_acids = ['X', ]
alphabet = set(nukes + [ambig_base])  # NOTE not the greatest naming distinction, but note difference to <expected_characters>
gap_chars = ['.', '-']
expected_characters = set(nukes + [ambig_base] + gap_chars)  # NOTE not the greatest naming distinction, but note difference to <alphabet>
conserved_codons = {l : {'v' : 'cyst',
                         'j' : 'tryp' if l == 'igh' else 'phen'}  # e.g. heavy chain has tryp, light chain has phen
                     for l in loci}
def get_all_codons():  # i'm only make these two fcns rather than globals since they use fcns that aren't defined til way down below
    return [''.join(c) for c in itertools.product('ACGT', repeat=3)]
def get_all_amino_acids(no_stop=False, include_ambig=False):  # i'm not sure why i didn't include ambig by default originally, but i don't want to change it now
    all_aas = set(ltranslate(c) for c in get_all_codons())  # note: includes stop codons (*)
    if no_stop:
        all_aas.remove('*')
    if include_ambig:
        all_aas |= set(ambiguous_amino_acids)
    return all_aas
def check_nuc_alphabet(seq, only_warn=False):
    check_alphabet(seq, aa=False)
def check_aa_alphabet(seq, only_warn=False):
    check_alphabet(seq, aa=True)
def check_alphabet(seq, aa=False, only_warn=False):
    if aa:
        alph = get_all_amino_acids(include_ambig=True) | set(gap_chars)
    else:
        alph = set(alphabet) | set(gap_chars)
    if len(set(seq) - set(alph)) > 0:
        fstr = 'unexpected %s chars[s]: %s (expected %s)' % ('amino acid' if aa else 'nuc', ' '.join(set(seq) - set(alph)), ' '.join(alph))
        if only_warn:
            print('  %s %s' % (color('yellow', 'warning'), fstr))
        else:
            raise Exception(fstr)

def cdn(glfo, region):  # returns None for d
    return conserved_codons[glfo['locus']].get(region, None)
def cdn_positions(glfo, region):
    if cdn(glfo, region) is None:
        return None
    return glfo[cdn(glfo, region) + '-positions']
def cdn_pos(glfo, region, gene):
    if cdn(glfo, region) is None:
        return None
    return cdn_positions(glfo, region)[gene]
def gseq(glfo, gene):  # adding this fcn very late in the game, i could really stand to use it in a lot more places (grep for "glfo\['seqs'\]\[.*\]\[")
    return glfo['seqs'][get_region(gene)][gene]
def remove_gaps(seq):
    return ''.join([c for c in seq if c not in gap_chars])
def gap_len(seq):  # NOTE see two gap-counting fcns below (_pos_in_alignment())
    return len([c for c in seq if c in gap_chars])
def non_gap_len(seq):  # NOTE see two gap-counting fcns below (_pos_in_alignment())
    return len(seq) - gap_len(seq)
def ambig_frac(seq, aa=False):
    assert not aa  # just to remind you that this fcn doesn't handle aa
    # ambig_seq = ''.join(filter(all_ambiguous_bases.__contains__, seq))
    # return float(len(ambig_seq)) / len(seq)
    return float(seq.count(ambig_base)) / len(seq)

# NOTE these fcns probably shouldn't exist -- what should really happen is we keep track of the regional/insertion boundaries for the aa seq as well... but that isn't really feasible, so for the moment we have this hack
# ----------------------------------------------------------------------------------------
def n_variable_ambig_nuc(line, nuc_seq):  # number of ambiguous bases in the variable region (start of v to end of j)
    check_nuc_alphabet(nuc_seq)
    istart, istop = line['regional_bounds']['v'][0], line['regional_bounds']['j'][1]
    return nuc_seq[istart : istop].count(ambig_base)
# ----------------------------------------------------------------------------------------
def n_variable_ambig_aa(line, aa_seq, nuc_seq):  # omfg don't ask why it's better to have two different fcns
    # NOTE this duplicates [is reverse of?] code in pad_seq_for_translation() (but not really sure how to clean up/combine them)
    check_aa_alphabet(aa_seq)
    istart = int(math.ceil(len(line['fv_insertion']) / 3.))  # add_seqs_aa() pads partial fv insertion codons, which will often [always?] be ambigous, so we don't want to include them here
    istop = len(aa_seq) - int(math.ceil(len(line['jf_insertion']) / 3.))  # on the right we also want to chop off any partial codons (since they'll usually [always?] be ambiguous
    if len(pad_seq_for_translation(line, nuc_seq)) % 3 != 0:  # furthermore, the 3' end of j seems to at least often not be a complete codon, which'll mean it's usually ambiguous, so we also need to get rid of it
        istop -=1
    return aa_seq[istart : istop].count(ambiguous_amino_acids[0])

# ----------------------------------------------------------------------------------------
def reverse_complement_warning():
    return '%s maybe need to take reverse complement (partis only searches in forward direction) or set --locus (default is igh). Both of these can be fixed using bin/split-loci.py. Or, maybe it\'s fine and your sequences just have super high mutation, which will be handled just fine in smith-waterman and the hmm.' % color('red', 'note:')

codon_table = {
    'cyst' : ['TGT', 'TGC'],
    'tryp' : ['TGG', ],
    'phen' : ['TTT', 'TTC'],
    'stop' : ['TAG', 'TAA', 'TGA']
}
# Infrastrucure to allow hashing all the columns together into a dict key.
# Uses a tuple with the variables that are used to index selection frequencies
# NOTE fv and jf insertions are *effective* (not real) insertions between v or j and the framework. They allow query sequences that extend beyond the v or j regions
index_columns = tuple(['v_gene', 'd_gene', 'j_gene', 'v_5p_del', 'v_3p_del', 'd_5p_del', 'd_3p_del', 'j_5p_del', 'j_3p_del', 'fv_insertion', 'vd_insertion', 'dj_insertion', 'jf_insertion'])

index_keys = {}
for i in range(len(index_columns)):  # dict so we can access them by name instead of by index number
    index_keys[index_columns[i]] = i

# don't think this has been used in a long time
# # ----------------------------------------------------------------------------------------
# def get_codon(fname):
#     codon = fname.split('-')[0]
#     if codon not in [c for locus in loci for c in conserved_codons[locus].values()]:
#         raise Exception('couldn\'t get codon from file name %s' % fname)
#     return codon

# ----------------------------------------------------------------------------------------
# Info specifying which parameters are assumed to correlate with which others. Taken from mutual
# information plot in bcellap repo

# key is parameter of interest, and associated list gives the parameters (other than itself) which are necessary to predict it
# TODO am I even using these any more? I think I might just be using the all-probs file
column_dependencies = {}
column_dependencies['v_gene'] = [] # NOTE v choice actually depends on everything... but not super strongly, so a.t.m. I ignore it
column_dependencies['v_5p_del'] = ['v_gene']
column_dependencies['v_3p_del'] = ['v_gene']
column_dependencies['d_gene'] = []
column_dependencies['d_5p_del'] = ['d_gene']  # NOTE according to the mebcell mutual information plot d_5p_del is also correlated to d_3p_del (but we have no way to model that a.t.m. in the hmm)
column_dependencies['d_3p_del'] = ['d_gene']  # NOTE according to the mebcell mutual information plot d_3p_del is also correlated to d_5p_del (but we have no way to model that a.t.m. in the hmm)
column_dependencies['j_gene'] = []
column_dependencies['j_5p_del'] = ['j_gene']  # NOTE mebcell plot showed this correlation as being small, but I'm adding it here for (a possibly foolish) consistency
column_dependencies['j_3p_del'] = ['j_gene']
column_dependencies['fv_insertion'] = []
column_dependencies['vd_insertion'] = ['d_gene']
column_dependencies['dj_insertion'] = ['j_gene']
column_dependencies['jf_insertion'] = []
# column_dependencies['vd_insertion_content'] = []
# column_dependencies['dj_insertion_content'] = []

# tuples with the column and its dependencies mashed together
# (first entry is the column of interest, and it depends upon the following entries)
column_dependency_tuples = []
for column, deps in column_dependencies.items():
    tmp_list = [column]
    tmp_list.extend(deps)
    column_dependency_tuples.append(tuple(tmp_list))

# pairs of parameters for which you're allowed to specify correlations in simulation (note that the behavior is highly dependent on the order in which we choose to simulate parameter choice in scratch rearrangement)
available_simu_correlations = [('v_gene', 'd_gene'), ('d_gene', 'j_gene'), ('v_gene', 'j_gene'),
                               ('v_gene', 'v_3p_del'),
                               ('d_gene', 'd_5p_del'), ('d_gene', 'd_3p_del'),
                               ('j_gene', 'j_5p_del'),
                               ('v_gene', 'vd_insertion'),
                               ('d_gene', 'vd_insertion'),
                               ('d_gene', 'dj_insertion'),
                               ('j_gene', 'dj_insertion'),
]
paired_available_simu_correlations = available_simu_correlations + [(p, p) for p in set(q for qpair in available_simu_correlations for q in qpair)]  # all the single-chain ones, plus each parameter with itself
paired_available_simu_correlations = [ppair for ppair in paired_available_simu_correlations if ppair[1] not in ['d_gene', 'd_5p_del', 'd_3p_del', 'vd_insertion']]  # doesn't make sense to have d gene stuff since there's only one option for light chain (and these are all correlations on choosing the light chain)

# ----------------------------------------------------------------------------------------
adaptive_headers = {
    'seqs' : 'nucleotide',
    'v_gene' : 'vMaxResolved',
    'd_gene' : 'dMaxResolved',
    'j_gene' : 'jMaxResolved',
    'v_3p_del' : 'vDeletion',
    'd_5p_del' : 'd5Deletion',
    'd_3p_del' : 'd3Deletion',
    'j_5p_del' : 'jDeletion'
}

# ----------------------------------------------------------------------------------------
# make new-style (unicode) translation dict from two lists of strings (i.e. from input to old-style string.maketrans())
def make_unicode_translation(str1, str2):
    return {ord(a) : ord(b) for a, b in zip(str1, str2)}

forbidden_characters = set([':', ';', ','])  # strings that are not allowed in sequence ids
forbidden_character_translations = make_unicode_translation(':;,', 'csm')
warn_chars = set('()=')  # these may cause issues with phylo inference programs
ambig_translations = make_unicode_translation(all_ambiguous_bases, ambig_base * len(all_ambiguous_bases))

functional_columns = ['mutated_invariants', 'in_frames', 'stops']

# ----------------------------------------------------------------------------------------
def useful_bool(bool_str):
    if bool_str == 'True':
        return True
    elif bool_str == 'False':
        return False
    elif bool_str == '1':
        return True
    elif bool_str == '0':
        return False
    else:
        raise Exception('couldn\'t convert \'%s\' to bool' % bool_str)

# ----------------------------------------------------------------------------------------
def get_z_scores(vals):  # return list of <vals> normalized to units of standard deviations (i.e. z scores)
    mean, std = numpy.mean(vals), numpy.std(vals, ddof=1)
    return [(v - mean) / std for v in vals]

# ----------------------------------------------------------------------------------------
def get_str_float_pair_dict(strlist):
    def getpair(pairstr):
        pairlist = pairstr.split(':')
        assert len(pairlist) == 2
        return (pairlist[0], float(pairlist[1]))
    return OrderedDict(getpair(pairstr) for pairstr in strlist)

# ----------------------------------------------------------------------------------------
def get_list_of_str_list(strlist):
    if strlist == '':
        return []
    return [[] if substr == '' else substr.split(':') for substr in strlist]

# keep track of all the *@*@$!ing different keys that happen in the <line>/<hmminfo>/whatever dictionaries
linekeys = {}
# I think 'per_family' is pretty incomplete at this point, but I also think it isn't being used
linekeys['per_family'] = ['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds', 'reco_id', 'clone_id'] + \
                         ['invalid', 'tree', 'consensus_seq', 'consensus_seq_aa', 'naive_seq_aa', 'cons_dists_nuc', 'cons_dists_aa'] + \
                         [r + '_gene' for r in regions] + \
                         [e + '_del' for e in all_erosions] + \
                         [b + '_insertion' for b in all_boundaries] + \
                         [r + '_gl_seq' for r in regions] + \
                         [r + '_per_gene_support' for r in regions]
# NOTE some of the indel keys are just for writing to files, whereas 'indelfos' is for in-memory
# note that, as a list of gene matches, all_matches would in principle be per-family, except that it's sw-specific, and sw is all single-sequence (this is also true of fv/jf insertions)
linekeys['per_seq'] = ['seqs', 'unique_ids', 'mut_freqs', 'n_mutations', 'shm_aa', 'input_seqs', 'indel_reversed_seqs', 'cdr3_seqs', 'cdr3_seqs_aa', 'full_coding_input_seqs', 'padlefts', 'padrights', 'indelfos', 'duplicates',
                       'leader_seqs', 'c_gene_seqs',  'leaders', 'c_genes', # these are kind of replacing fv/jf insertions, and the latter probably should just be removed, since they're really more per-seq things, but i don't know if it'd really work, and it'd for sure be hard, so whatever (otoh, maybe fv/jf insertions are necessary for padding to same length in sw? not sure atm)
                       'has_shm_indels', 'qr_gap_seqs', 'gl_gap_seqs', 'loci', 'paired-uids', 'all_matches', 'seqs_aa', 'input_seqs_aa', 'cons_dists_nuc', 'cons_dists_aa', 'lambdas', 'nearest_target_indices',
                       'min_target_distances', 'vrc01-muts', 'mut_positions', 'subjects'] + \
                      [r + '_qr_seqs' for r in regions] + \
                      ['aligned_' + r + '_seqs' for r in regions] + \
                      functional_columns
linekeys['hmm'] = ['logprob', 'errors', 'tree-info', 'alternative-annotations'] + [r + '_per_gene_support' for r in regions]
linekeys['sw'] = ['k_v', 'k_d', 'all_matches', 'padlefts', 'padrights']
linekeys['simu'] = ['reco_id', 'affinities', 'relative_affinities', 'lambdas', 'tree', 'target_seqs', 'nearest_target_indices', 'min_target_distances', 'heavy-chain-correlation-info']

# keys that are required to specify a naive rearrangement event
minimal_linekeys = [r+'_gene' for r in regions] + [b+'_insertion' for b in boundaries] + [e+'_del' for e in real_erosions]
# keys that are added by add_implicit_info() (these are *not* all the non-minimal ones)
implicit_linekeys = set(['naive_seq', 'cdr3_length', 'codon_positions', 'lengths', 'regional_bounds', 'invalid', 'indel_reversed_seqs'] + \
                        [r + '_gl_seq' for r in regions] + \
                        ['mut_freqs', 'n_mutations'] + functional_columns + [r + '_qr_seqs' for r in regions] + ['aligned_' + r + '_seqs' for r in regions])

extra_annotation_headers = [  # you can specify additional columns (that you want written to csv) on the command line from among these choices (in addition to <annotation_headers>)
    'cdr3_seqs',
    'full_coding_naive_seq',
    'full_coding_input_seqs',
    'linearham-info',
    'consensus_seq',
    'consensus_seq_aa',
    'naive_seq_aa',
    'seqs_aa',
    'cdr3_seqs_aa',
    'cons_dists_nuc',
    'cons_dists_aa',
    'tree',
] + list(implicit_linekeys)  # NOTE some of the ones in <implicit_linekeys> are already in <annotation_headers>

fake_paired_columns = ['is_fake_paired', 'h_seqs', 'l_seqs']

linekeys['extra'] = extra_annotation_headers
all_linekeys = set([k for cols in linekeys.values() for k in cols])
def add_per_seq_keys(line):  # NOTE should really combine this with add_input_meta_keys()
    pskeys = set(line) & set(linekeys['per_seq'])
    other_keys = [k for k in line if k not in pskeys and isinstance(line[k], list) and len(line[k])==len(line['unique_ids'])]  # custom meta info keys (yeah they should get added some other way, but it's hard to do it super reliably
    if len(other_keys) > 0:
        print('    %s adding missing per-seq key%s: %s' % (wrnstr(), plural(len(other_keys)), ', '.join(other_keys)))
        add_input_meta_keys(other_keys, are_line_keys=True)

# ----------------------------------------------------------------------------------------
forbidden_metafile_keys = ['name', 'seq']  # these break things if you add them
input_metafile_keys = {  # map between the key we want the user to put in the meta file, and the key we use in the regular <line> dicts (basically just pluralizing)
    'affinity' : 'affinities',  # should maybe add all of these to <annotation_headers>?
    'relative_affinity' : 'relative_affinities',
    'timepoint' : 'timepoints',
    'multiplicity' : 'multiplicities',
    'subject' : 'subjects',
    'constant-region' : 'constant-regions',
    'paired-uids' : 'paired-uids',
    'locus' : 'loci',  # i *think* it's worth having this, even though it duplicates info in the glfo, since it's nice to be able to pass in locus info for some input sequences (e.g. --queries-to-include-fname), but you should generally avoid using the 'loci' key in annotations unless you really need it, since it's annoying trying to make sure it's always there)
    'cell-type' : 'cell-types',
    'umis' : 'umis',
    'reads' : 'reads',
    'c_gene' : 'c_genes',
    'gc-round' : 'gc-rounds',
    'generation-time' : 'generation-times',
    # 'chosen' : 'chosens',
    # 'paired' : 'paireds',
}
linekeys['per_seq'] += list(input_metafile_keys.values())
def input_metafile_defaults(mkey):  # default values to use if the info isn't there (have to use a fcn rather than dict so e.g. [] doesn't give the same object each time)
    if mkey == 'multiplicities':
        return 1
    elif mkey == 'paired-uids':
        return []
    else:
        return None

reversed_input_metafile_keys = {v : k for k, v in input_metafile_keys.items()}

# ----------------------------------------------------------------------------------------
# add any keys in <meta_keys> that aren't in input_metafile_keys
# NOTE should really combine this with add_per_seq_keys()
def add_input_meta_keys(meta_keys, are_line_keys=False):  # NOTE I'm adding this late, and not completely sure that it's ok to modify these things on the fly like this (but i think the ability to add arbitrary keys is super important)
    extra_keys = set(meta_keys) - set(input_metafile_keys) - set(forbidden_metafile_keys) - set(linekeys['per_seq'])
    # should really also check that it's a list with the right len here (like in add_per_seq_keys), but i should combine them anyway
    if len(extra_keys) == 0:
        return
    new_keys = []
    for ekey in extra_keys:
        if are_line_keys:
            if ekey[-1] != 's':
                print('  %s last char of input meta key \'%s\' isn\'t \'s\'' % (color('yellow', 'warning'), ekey))
            ekey = ekey[:-1]  # if we could, we'd use reversed_input_metafile_keys
        input_metafile_keys[ekey] = ekey + 's'
        reversed_input_metafile_keys[ekey+'s'] = ekey
        linekeys['per_seq'].append(ekey+'s')
        sw_cache_headers.append(ekey+'s')
        annotation_headers.append(ekey+'s')
        new_keys.append(ekey)
    print('  note: added %d new input meta key%s to allowed keys (add \'s\'/plural to access it in the final annotations): %s (%s)' % (len(new_keys), plural(len(new_keys)), ' '.join(sorted(new_keys)), ' '.join(k+'s' for k in sorted(new_keys))))

# ----------------------------------------------------------------------------------------
def meta_emph_arg_process(args):
    if args.meta_info_to_emphasize is not None:
        add_input_meta_keys(set(args.meta_info_to_emphasize), are_line_keys=True)  # these would get added when we read the meta info file, if we read it, but we don't e.g. when reading existing output
        if args.meta_emph_formats is not None:  # have to make it so calling len() on this *and* on the actual values from annotations are comparable (ick, this is ugly, but... that's plotting for you)
            for mkey, mval in args.meta_info_to_emphasize.items():
                if args.meta_emph_formats.get(mkey) == 'len':
                    args.meta_info_to_emphasize[mkey] = [None for _ in range(int(mval))]
        assert len(args.meta_info_to_emphasize) == 1  # should at some point let there be more than one key
    if args.meta_info_key_to_color is not None:
        add_input_meta_keys([args.meta_info_key_to_color], are_line_keys=True)

# ----------------------------------------------------------------------------------------
def meta_emph_str(key, val, formats=None, legend_groups=None):  # ick
    kstr, use_len = key.rstrip('s'), False  # (rstrip removes plural, and yes will probably break something at some point)
    if formats is not None and key in formats:
        if formats[key] == 'len':
            use_len = True
        # elif formats[key] == 'translate':  # ugh, something like this would be better
        elif '=' in formats[key]:
            old_val, new_val = formats[key].split('=')
            if val == old_val:
                val = new_val
        else:
            kstr = formats[key].replace('@', ' ')
    if use_len:
        return 'None' if val is None else '%d' % len(val)
    elif isinstance(val, float):
        return '%.2f' % val
    elif isinstance(val, bool) or val in ['True', 'False']:
        return kstr if val else 'nope'  # not sure what to do if it's False, but probably you'll only call it asking for True?
    else:
        return str(val) if legend_groups is None else legend_groups.get(val, val)

# ----------------------------------------------------------------------------------------
def meta_info_equal(key, val1, val2, formats=None):  # also ick (don't have another way to convert the value from the command line arg to the proper thing [the values coming from yaml file will be properly converted)
    def cfcn():
        if formats is not None and formats.get(key) == 'len':
            return len
        for tname in [bool, float, int]:
            if any(isinstance(v, tname) for v in [val1, val2]):  # if either value is float/int convert both values to float/int (for bool use str)
                if tname == bool:
                    bstrs = ['True', 'False']
                    if any(str(v) not in bstrs for v in [val1, val2]):
                        raise Exception('bool string[s] %s not among %s' % ([v for v in [val1, val2] if str(v) not in bstrs], bstrs))
                    return str
                else:
                    return tname
        return lambda x: x
    # ----------------------------------------------------------------------------------------
    return cfcn()(val1) == cfcn()(val2)

# ----------------------------------------------------------------------------------------
# return dict of "meta emph value" : [fraction of cluster], e.g. for cell types could be {'pb' : 0.9, 'naive' : 0.1}
def get_meta_emph_fractions(mekey, all_emph_vals, cluster, antn, formats=None):
    def gfcn(x): return meta_emph_str(mekey, per_seq_val(antn, mekey, x, use_default=True), formats=formats)
    vgroups = group_seqs_by_value(cluster, gfcn, return_values=True)
    emph_fracs = {v : len(grp) / float(len(cluster)) for v, grp in vgroups}
    emph_fracs.update({v : 0. for v in all_emph_vals - set(emph_fracs)})  # need to include families with no seqs with a given value (otherwise all the lines go to 1 as cluster size goes to 1)
    return emph_fracs

# ----------------------------------------------------------------------------------------
special_indel_columns_for_output = ['qr_gap_seqs', 'gl_gap_seqs', 'indel_reversed_seqs']  # arg, ugliness (but for reasons...)  NOTE used to also include 'has_shm_indels' (also note that 'indel_reversed_seqs' is treated differently for some purposes than are the gap seq keys)
annotation_headers = ['unique_ids', 'invalid', 'v_gene', 'd_gene', 'j_gene', 'cdr3_length', 'mut_freqs', 'n_mutations', 'input_seqs', 'indel_reversed_seqs', 'has_shm_indels', 'qr_gap_seqs', 'gl_gap_seqs', 'naive_seq', 'duplicates', 'leader_seqs', 'c_gene_seqs', 'leaders', 'c_genes'] \
                     + [r + '_per_gene_support' for r in regions] \
                     + [e + '_del' for e in all_erosions] + [b + '_insertion' for b in all_boundaries] \
                     + functional_columns + list(input_metafile_keys.values()) \
                     + ['codon_positions', 'tree-info', 'alternative-annotations']
simulation_headers = linekeys['simu'] + [h for h in annotation_headers if h not in linekeys['hmm']]
sw_cache_headers = [h for h in annotation_headers if h not in linekeys['hmm']] + linekeys['sw']
partition_cachefile_headers = ('unique_ids', 'logprob', 'naive_seq', 'naive_hfrac', 'errors')  # these have to match whatever bcrham is expecting (in packages/ham/src/glomerator.cc, ReadCacheFile() and WriteCacheFile())
bcrham_dbgstrs = {
    'partition' : {  # corresponds to stdout from glomerator.cc
        'read-cache' : ['logprobs', 'naive-seqs'],
        'calcd' : ['vtb', 'fwd', 'hfrac'],
        'merged' : ['hfrac', 'lratio'],
        'time' : ['bcrham', ]
    },
    'annotate' : {  # corresponds to stdout from, uh, trellis.cc or something
        'calcd' : ['vtb', 'fwd'],
        'time' : ['bcrham', ]
    },
}
bcrham_dbgstr_types = {
    'partition' : {
        'sum' : ['calcd', 'merged'],  # for these ones, sum over all procs
        'same' : ['read-cache', ],  # check that these are the same for all procs
        'min-max' : ['time', ]
    },
    'annotate' : {  # strict subset of the 'partition' ones
        'sum' : ['calcd', ],
        'same' : [],
        'min-max' : ['time', ]
    },
}

# ----------------------------------------------------------------------------------------
io_column_configs = {
    'ints' : ['n_mutations', 'cdr3_length', 'padlefts', 'padrights'] + [e + '_del' for e in all_erosions],
    'floats' : ['logprob', 'mut_freqs'],
    'bools' : functional_columns + ['has_shm_indels', 'invalid'],
    'literals' : ['indelfo', 'indelfos', 'k_v', 'k_d', 'all_matches', 'alternative-annotations'],  # simulation has indelfo[s] singular, annotation output has it plural... and I think it actually makes sense to have it that way
    'lists-of-lists' : ['duplicates'] + [r + '_per_gene_support' for r in regions],  # NOTE that some keys are in both 'lists' and 'lists-of-lists', e.g. 'duplicates'. This is now okish, since the ordering in the if/else statements is now the same in both get_line_for_output() and process_input_line(), but it didn't used to be so there's probably files lying around that have it screwed up). It would be nice to not have it this way, but (especially for backwards compatibility) it's probably best not to mess with it.
    'lists' : [k for k in linekeys['per_seq'] if k not in ['indelfos', 'all_matches']],  # indelfos and all_matches are lists, but we can't just split them by colons since they have colons within the dict string
}

# NOTE these are all *input* conversion functions (for ouput we mostly just call str())
# also NOTE these only get used on deprecated csv output files
conversion_fcns = {}
for key in io_column_configs['ints']:
    conversion_fcns[key] = int
for key in io_column_configs['floats']:
    conversion_fcns[key] = float
for key in io_column_configs['bools']:
    conversion_fcns[key] = useful_bool
for key in io_column_configs['literals']:
    conversion_fcns[key] = ast.literal_eval
for region in regions:
    conversion_fcns[region + '_per_gene_support'] = get_str_float_pair_dict
conversion_fcns['duplicates'] = get_list_of_str_list

# ----------------------------------------------------------------------------------------
def nullval(key):  # adding this very late, so could probably stand to be used in lots of other places
    from . import indelutils
    if key == 'indelfos':
        return indelutils.get_empty_indel()
    else:
        return None

# ----------------------------------------------------------------------------------------
def cluster_size_str(partition, split_strs=False, only_passing_lengths=False, clusters_to_emph=None, short=False):  # NOTE doesn't need to be a partition, i.e. can have duplicate clusters
    def cl_str(c):
        rstr = str(len(c))
        if clusters_to_emph is not None and c in clusters_to_emph:
            rstr = color('blue', rstr)
        return rstr
    if split_strs:  # set this if you're passing in a list of ':'.joined clusters, rather than actual clusters
        partition = [cstr.split(':') for cstr in partition]
        if clusters_to_emph is not None:
            clusters_to_emph = [cstr.split(':') for cstr in clusters_to_emph]
    if only_passing_lengths:  # <partition> just contains cluster lengths, not actual clusters
        partition = [['x' for _ in range(l)] for l in partition]
    ptn_strs = [cl_str(c) for c in sorted([c for c in partition if len(c)>1], key=len, reverse=True)]
    fstr = ' '.join(ptn_strs)
    singletons = [c for c in partition if len(c)==1]
    if len(singletons) > 0:
        emph_str = ''
        if clusters_to_emph is not None and any(len(c)==1 for c in clusters_to_emph) and any(c in clusters_to_emph for c in singletons):  # second term is just a speed optimization
            emph_str = ', %s' % color('blue', '%d emph'%len([c for c in singletons if c in clusters_to_emph]))
        fstr = '%s (+%d%s%s)' % (fstr, len(singletons), '' if short else ' singletons', emph_str)
    return fstr

# ----------------------------------------------------------------------------------------
# avoid numpy's "ragged nested sequences" warning (since it converts everything to ndarrays)
def np_rand_choice(tarray, size=None, replace=True, probs=None):
    ichosen = numpy.random.choice(range(len(tarray)), size=size, replace=replace, p=probs)
    return [tarray[i] for i in ichosen]

# ----------------------------------------------------------------------------------------
# this apportionment is basically copied from the subcluster annotation code in partitiondriver.py, although there we start with group size rather than <n_subsets>
def evenly_split_list(inlist, n_subsets, dbgstr='items'):
    n_per_group = math.floor(len(inlist) / n_subsets)
    n_seq_list = [n_per_group for _ in range(n_subsets)]
    for iextra in range(len(inlist) % n_subsets):  # spread out any extra ones
        n_seq_list[iextra] += 1
    list_of_lists, iglobal = [], 0
    for isub in range(n_subsets):
        list_of_lists.append(inlist[iglobal : iglobal + n_seq_list[isub]])
        iglobal += n_seq_list[isub]
    print('  splitting %d %s into %d groups with sizes: %s' % (len(inlist), dbgstr, n_subsets, ' '.join('%d'%len(g) for g in list_of_lists)))
    assert sum(len(g) for g in list_of_lists) == len(inlist)
    return list_of_lists

# ----------------------------------------------------------------------------------------
# chooses [sequences from] N droplets according to <n_max_queries> or <n_random_queries> (has to first group them by droplet id so that we choose both h and l seq together)
# n_max_queries: take only this many (after random shuffling, see notes below), n_random_queries: take this many at random, n_subsets: split into this many groups
# NOTE returns list of lists of seqfos if n_subsets is set (rather than a single list of seqfos)
def subset_paired_queries(seqfos, droplet_id_separators, droplet_id_indices, n_max_queries=-1, n_random_queries=None, n_subsets=None, input_info=None,
                          seed_unique_ids=None, queries_to_include=None, debug=False):  # yes i hate that they have different defaults, but it has to match the original partis arg, which i don't want to change
    # ----------------------------------------------------------------------------------------
    def get_final_seqfos(final_qlists, dbgarg, dbgstr):
        final_sfos = [sfdict[u] for l in final_qlists for u in l]  # this loses the original order from <seqfos>, but the original way I preserved that order was super slow, and whatever who cares
        print('  %s: chose %d / %d droplets %s which had %d / %d seqs' % (dbgarg, len(final_qlists), len(drop_query_lists), dbgstr, len(final_sfos), sum(len(ql) for ql in drop_query_lists)))
        if seed_unique_ids is not None:  # kind of similar to seqfileopener.post_process()
            missing_ids = sorted(set(seed_unique_ids) - set(s['name'] for s in final_sfos))
            final_sfos += [sfdict[u] for u in missing_ids]
            print('      --seed-unique-id: added %d/%d seed ids (others already there)' % (len(missing_ids), len(seed_unique_ids)))
        if queries_to_include is not None:
            missing_ids = sorted(set(queries_to_include) - set(s['name'] for s in final_sfos))
            final_sfos += [sfdict[u] for u in missing_ids]
            print('      --queries-to-include: added %d/%d queries (others already there)' % (len(missing_ids), len(queries_to_include)))
        return final_sfos
    # ----------------------------------------------------------------------------------------
    sfdict = {s['name'] : s for s in seqfos}
    pvals = [n_max_queries, n_random_queries, n_subsets]
    if pvals == [-1, None, None] or pvals.count(None) + pvals.count(-1) != 2:  # not really correct (-1 isn't default for the second two) but why would you set them to that, anyway?
        raise Exception('have to set exactly 1 of [n_max_queries, n_random_queries, n_subsets], but got %s %s %s' % (n_max_queries, n_random_queries, n_subsets))
    if input_info is None:  # default: group seqfos by splitting each uid into droplet id etc
        _, drop_query_lists = get_droplet_groups([s['name'] for s in seqfos], droplet_id_separators, droplet_id_indices, return_lists=True, debug=debug)
    else:  # but if we already have pair info in <input_info>
        _, drop_query_lists = get_droplet_groups_from_pair_info(input_info, droplet_id_separators, droplet_id_indices, return_lists=True)
    if any(u not in sfdict for ql in drop_query_lists for u in ql):  # shouldn't happen any more (split loci wasn't properly removing uids), but probably still a good check to have
        missing_ids = [u for ql in drop_query_lists for u in ql if u not in sfdict]
        raise Exception('%d drop query list uid[s] missing from sfdict (size %d): %s' % (len(missing_ids), len(sfdict), ' '.join(missing_ids)))
    random.shuffle(drop_query_lists)  # we always want to shuffle, since sometimes e.g. all unpaired igh are first when alphabetized by droplet id
    if n_max_queries != -1:  # NOTE not same order as input file, since that wouldn't/doesn't make any sense
        final_qlists = drop_query_lists[:n_max_queries]
        returnfo = get_final_seqfos(final_qlists, '--n-max-queries', '(first %d, after random shuffling)' % n_max_queries)
    elif n_random_queries is not None:
        final_qlists = np_rand_choice(drop_query_lists, size=n_random_queries, replace=False)
        returnfo = get_final_seqfos(final_qlists, '--n-random-queries', '(uniform randomly)')
    elif n_subsets is not None:
        list_of_fqlists = evenly_split_list(drop_query_lists, n_subsets, dbgstr='droplets')
        returnfo = []
        for isub, final_qlists in enumerate(list_of_fqlists):
            returnfo.append(get_final_seqfos(final_qlists, '  isub %d'%isub, 'in order (after random shuffling)'))
    else:
        assert False
    return returnfo

# ----------------------------------------------------------------------------------------
# NOTE duplicates code in paircluster.clean_pair_info()
# TODO better name
def get_droplet_groups_from_pair_info(input_info, droplet_id_separators, droplet_id_indices, return_lists=False):
    drop_ids = []   # all droplet ids
    pid_groups = []  # list of pid groups, i.e. each element is the uids from a single droplet (for 10x)
    pid_ids = {}  # map from each uid to the index of its pid group
    n_no_info = 0
    for uid, tline in input_info.items():
        assert len(tline['unique_ids']) == 1
        pids = tline.get('paired-uids', [[]])[0]
        if any(p not in input_info for p in pids):  # make sure there aren't any ids in paired-uids that don't have their own entry in input_info (probably they were incompletely removed when e.g. they failed an alignment/quality control step)
            miss_pids = [p for p in pids if p not in input_info]
            pids = [p for p in pids if p in input_info]
            # print('    ignoring %d pids not in input info: %s' % (len(miss_pids), ' '.join(miss_pids)))
        if len(pids) == 0:
            n_no_info += 1
        pset = set([uid] + pids)
        found = False
        for pd in pids:  # even if we assume the pair info from all seqs is consistent (which this more or less does), we need to do this to get the <ipg>
            if pd in pid_ids:
                ipg = pid_ids[pd]
                found = True
                pid_groups[ipg] |= pset  # don't really need this, but i guess it catches one possible way in which pair info could be inconsistent? Although we should really probably crash if so
                break
        if not found:
            pid_groups.append(pset)
            ipg = len(pid_groups) - 1
            drop_ids.append(get_droplet_id(uid, droplet_id_separators, droplet_id_indices)) #ipg if droplet_id_separators is None else get_droplet_id(uid, droplet_id_separators, droplet_id_indices)
        assert ipg is not None
        for pid in pset:  # this does a lot of overwriting, but it shouldn't matter
            pid_ids[pid] = ipg
    print('  grouped %d input seqs into %d droplet groups using existing pair info (first few droplet ids: %s)' % (len(input_info), len(pid_groups), ' '.join(drop_ids[:3])))
    if n_no_info > 0:
        print('     %d/%d (%.2f) had no pair info' % (n_no_info, len(input_info), n_no_info / float(len(input_info))))
    for ipg, pg in enumerate(pid_groups):
        #     print '  %3d %s' % (ipg, ' '.join(pg))
        pid_groups[ipg] = sorted(pg)  # need to sort (and convert from set to list) for replicability
    if return_lists:
        return drop_ids, pid_groups
    else:
        return [(did, dql) for did, dql in zip(drop_ids, pid_groups)]

# ----------------------------------------------------------------------------------------
def get_droplet_groups(all_uids, droplet_id_separators, droplet_id_indices, return_lists=False, debug=False):  # group all queries into droplets
    def kfcn(u): return get_droplet_id(u, droplet_id_separators, droplet_id_indices)
    idg_it = itertools.groupby(sorted(all_uids, key=kfcn), key=kfcn)  # iterate over result with: 'for dropid, drop_queries in '
    drop_ids, drop_query_lists = [list(x) for x in zip(*[(did, list(qiter)) for did, qiter in idg_it])]
    if '' in drop_ids:
        n_before = len(drop_ids)
        i_did = drop_ids.index('')
        for qry in drop_query_lists[i_did]:
            drop_ids.append(qry)
            drop_query_lists.append([qry])
        drop_ids.pop(i_did)
        drop_query_lists.pop(i_did)
        print('    %s found empty droplet id when grouping sequences, so using these sequences\' uid as their droplet id (added %d droplet ids)' % (wrnstr(), len(drop_ids) - n_before))
    if debug:
        print('  grouped %d uids into %d droplet groups by splitting uids, then sorting droplet ids (first few droplet ids: %s)' % (len(all_uids), len(drop_ids), ' '.join(drop_ids[:3])))
    if return_lists:  # return separate lists
        return drop_ids, drop_query_lists
    else:  # return list of (droplet id, query list) pairs
        return [(did, dql) for did, dql in zip(drop_ids, drop_query_lists)]

# ----------------------------------------------------------------------------------------
def set_did_vals(uid):  # NOTE this is a pretty shitty way to do this, but i haven't thought of anything better
    if '_contig_' in uid:  # 10x data
        return '_', [0]
    elif uid.count('-') > 0 and uid.split('-')[-1] in loci:  # simulation, e.g. 7842229920653554932-igl
        if '_contig_' in uid:
            print('  %s found string \'_contig_\' in uid \'%s\' that also has locus \'%s\' at the end (which likely means pair info guessing is not working' % (wrnstr(), uid, uid.split('-')[-1]))
        return '-', list(range(uid.count('-')))  # make the last bit the contig id, and keep everything else as the droplet id (kind of arbitrary)
    else:  # default (same as 10x atm)
        return '_', [0]

# ----------------------------------------------------------------------------------------
did_help = {
    'guess' : 'Try to guess pairing info based on paired sequences having similar names, then append the locus to each uid (so e.g. \'2834_contig_1\', \'2834_contig_3\' --> \'2834-igh\', \'2834-igk\'). If set, this *must* be set both when split-loci.py is run (e.g. caching parameters) *and* during later inference (since e.g. --treefname or --input-partition-fname must be modified).',
    'seps' : 'String consisting of characters that should be used to find the droplet id in each uid string (with --droplet-id-index), i.e. calls <uid>.split(<--droplet-id-separators>)[--droplet-id-index] (recursively, if there\'s more than one character). For example \'-._\' would split on each instance of \'-\', \'.\', or \'_\'. Default is set in utils.set_did_vals(), and should work for 10x data and partis simulation (e.g. the 10x uid AAACGGGCAAGCGAGT-1_contig_2 has a droplet id of AAACGGGCAAGCGAGT-1, and the simulation uid of 7842229920653554932-igl has a droplet id of 7842229920653554932).',
    'indices' : 'Colon-separated list of zero-based indices used with --droplet-id-separators (see that arg\'s help). Default is set in utils.set_did_vals(), and should work for 10x data and partis simulation.',
}

# ----------------------------------------------------------------------------------------
# did_seps: string consisting of characters with which to split (i.e. if '_.-' we'll split on all three of them)
# did_indices: list of zero-based indices (using <did_seps>) of droplet id
def get_droplet_id(uid, did_seps, did_indices, return_contigs=False, return_colored=False, debug=False):
    auto_set = False
    if did_seps is None or did_indices is None:
        if did_seps is not None or did_indices is not None:
            print('  %s only one of did_seps/did_indices was set, so still need to guess both of them: %s %s' % (wrnstr(), did_seps, did_indices))
        did_seps, did_indices = set_did_vals(uid)
        auto_set = True
    if not any(s in uid for s in did_seps):  # probably not paired data
        if debug:
            print('  %s no separators %s in uid %s' % (wrnstr(), did_seps, uid))
        return (uid, uid) if return_contigs else uid
    ulist = [uid]
    for sep in did_seps:  # recursively split by each sep character
        ulist = [s for u in ulist for s in u.split(sep)]
    # update: i'm no longer requiring that the drople indices don't include the last chunk, now i'm letting the contig id be empty
    # if any(i > len(ulist) - 2 for i in did_indices):
    #     raise Exception('droplet id indices %s (out of %s) greater than len-1 (%d) of list %s after splitting by separators \'%s\' for uid \'%s\'' % ([i for i in did_indices if i > len(ulist) - 2], did_indices, len(ulist) - 1, ulist, did_seps, uid))
    did = did_seps[0].join(ulist[i] for i in did_indices)  # rejoin with just the first sep (if there was more than one), since doing otherwise would be complicated and i don't think it matters
    if did == '':
        raise Exception('empty droplet id for uid \'%s\'' % uid)
    cid = ulist[-1] if max(did_indices) < len(ulist) - 1 else ''  # normally just set the contig id to the last element, which is correct for [current] 10x data, except when we need the last element as part of the droplet id
    # , and we don't care about it otherwise (NOTE if you change this, you'll have to update the return_colored stuff below)
    if debug:
        print('    droplet id separators%s: %s  indices: %s' % (' (set automatically)' if auto_set else '', did_seps, did_indices))
        print('       e.g. uid \'%s\' --> droplet id \'%s\' contig id \'%s\'' % (uid, did, cid))
    if return_colored:
        return did_seps[0].join(color('red' if i in did_indices else ('blue' if i==len(ulist)-1 else None), ulist[i]) for i, s in enumerate(ulist))
    elif return_contigs:
        return did, cid  # NOTE returning cid as string (rather than int)
    else:
        return did

# ----------------------------------------------------------------------------------------
def get_contig_id(uid, did_seps, did_indices, debug=False):
    return get_droplet_id(uid, did_seps, did_indices, return_contigs=True, debug=debug)[1]

# ----------------------------------------------------------------------------------------
def is_correctly_paired(uid, pid):  # return true if <pid> is actually the correct partner for <uid> (i.e. this is simulation)
    # could use the droplet id fcn, but if you accidentally ran on data it'd cause havoc so i like this better
    def rmlocus(utmp):
        ltmp = get_single_entry([l for l in loci if '-'+l in utmp])
        return utmp.replace('-'+ltmp, '')
    return rmlocus(uid) == rmlocus(pid)  # true if everything except '-'+(any locus string) is the same

# ----------------------------------------------------------------------------------------
def extract_pairing_info(seqfos, droplet_id_separators=None, droplet_id_indices=None, input_metafname=None, debug=True):
    # ----------------------------------------------------------------------------------------
    def did_fcn(uid, return_colored=False, debug=False):  # shorthand that only requires uid
        return get_droplet_id(uid, did_seps=droplet_id_separators, did_indices=droplet_id_indices, return_colored=return_colored, debug=debug)
    # ----------------------------------------------------------------------------------------
    droplet_ids = {}
    for sfo in seqfos:
        did = did_fcn(sfo['name'])
        if did not in droplet_ids:
            droplet_ids[did] = []
        droplet_ids[did].append(sfo['name'])

    print('  extract_pairing_info(): read %d sequences with %d droplet ids' % (len(seqfos), len(droplet_ids)))
    if debug:
        did_fcn(seqfos[0]['name'], debug=True)
    count_info = {}
    for did, dlist in droplet_ids.items():
        if len(dlist) not in count_info:
            count_info[len(dlist)] = 0
        count_info[len(dlist)] += 1
    print('    contigs per')
    print('      droplet     count   fraction')
    total = sum(count_info.values())
    for size, count in sorted(list(count_info.items()), key=operator.itemgetter(0)):
        frac = count / float(total)
        print('       %2d        %5d     %s' % (size, count, ('%.3f'%frac) if frac > 0.001 else ('%.1e'%frac)))

    input_metafos = {}
    if input_metafname is not None:
        with open(input_metafname) as mfile:
            input_metafos = json.load(mfile)

    metafos = {}
    for sfo in seqfos:
        if sfo['name'] in metafos:
            raise Exception('duplicate uid \'%s\' when extracting pairing info' % sfo['name'])
        pids = [u for u in droplet_ids[did_fcn(sfo['name'])] if u != sfo['name']]
        metafos[sfo['name']] = {'paired-uids' : pids}
        if sfo['name'] in input_metafos:
            assert 'paired-uids' not in input_metafos[sfo['name']]  # don't want to adjudicate between alternative versions here
            metafos[sfo['name']].update(input_metafos[sfo['name']])

    if debug > 1:
        print('           did       N     uids')
        max_len = max(len(u) for u in metafos)
        dgpairs = get_droplet_groups(list(metafos.keys()), droplet_id_separators, droplet_id_indices)
        for did, dgroup in sorted(dgpairs, key=lambda x: len(x[1]), reverse=True):
            print('      %10s   %3d   %s' % (did, len(dgroup), '   '.join([did_fcn(u, return_colored=True)+(max_len-len(u))*' ' for u in dgroup])))

    return metafos

# ----------------------------------------------------------------------------------------
def check_concordance_of_cpath_and_annotations(cpath, annotation_list, annotation_dict, use_last=False, debug=False):
    if len(cpath.partitions) == 0:
        return
    if len(annotation_list) == 0:
        print('    %s zero length annotation list in check_concordance_of_cpath_and_annotations()' % color('yellow', 'warning'))
        return
    pfcn = 'last' if use_last else 'best'
    partition = getattr(cpath, pfcn)()
    miss_clusts = [c for c in partition if ':'.join(c) not in annotation_dict]
    extra_antns = [l for l in annotation_list if l['unique_ids'] not in partition]
    if len(miss_clusts) == 0 and debug:
        estr = '' if len(extra_antns)==0 else ' (%d extra%s with size%s: %s)' % (len(extra_antns), plural(len(extra_antns)), plural(len(extra_antns)), ' '.join(str(len(l['unique_ids'])) for l in extra_antns))
        print('    check_concordance_of_cpath_and_annotations(): annotations for all %d clusters in %s partition%s' % (len(partition), pfcn, estr))
    if len(miss_clusts) > 0:
        def pstr(ptn): return ', '.join(str(len(c)) for c in sorted(ptn, key=len, reverse=True))
        present_str = '' if len(partition)>30 else ' (present: %s)' % pstr([c for c in partition if c not in miss_clusts])
        print('    %s check_concordance_of_cpath_and_annotations(): missing annotations for %d/%d clusters in %s partition with size%s: %s%s' % (color('yellow', 'warning'), len(miss_clusts), len(partition), pfcn, plural(len(miss_clusts)), pstr(miss_clusts), present_str))
        if len(extra_antns) > 0:
            print('        %d extra annotations (i.e. they\'re for clusters not in %s partition) with sizes: %s' % (len(extra_antns), pfcn, ' '.join(str(len(l['unique_ids'])) for l in extra_antns)))
    for mclust in miss_clusts:
        overlap_antns = [l for l in annotation_list if len(set(mclust) & set(l['unique_ids'])) > 0]
        if len(overlap_antns) > 0:
            print('          missing \'%s\' cluster with size %d %s overlaps with the following (overlap/size): %s   (%s)' % (pfcn, len(mclust), ':'.join(mclust), '  '.join('%d/%d'%(len(set(mclust) & set(ol['unique_ids'])), len(ol['unique_ids'])) for ol in overlap_antns), ', '.join(':'.join(l['unique_ids']) for l in overlap_antns)))
    if len(miss_clusts) > 0:  # don't really care about extra ones unless there's missing ones that we might want to account for
        for eline in extra_antns:
            overlap_clusts = [c for c in partition if len(set(eline['unique_ids']) & set(c)) > 0]
            if len(overlap_clusts) > 0:
                print('          extra annotation with size %d %s overlaps with the following %s clusters (overlap/size): %s   (%s)' % (len(eline['unique_ids']), ':'.join(eline['unique_ids']), pfcn, '  '.join('%d/%d'%(len(set(eline['unique_ids']) & set(c)), len(c)) for c in overlap_clusts), ' '.join(':'.join(c) for c in overlap_clusts)))

# ----------------------------------------------------------------------------------------
def get_annotation_dict(annotation_list, duplicate_resolution_key=None, cpath=None, use_last=False, ignore_duplicates=False, quiet=False):
    annotation_dict, dup_keys = OrderedDict(), []
    for line in annotation_list:
        uidstr = ':'.join(line['unique_ids'])
        if uidstr in annotation_dict or duplicate_resolution_key is not None:  # if <duplicate_resolution_key> was specified, but it wasn't specified properly (e.g. if it's not sufficient to resolve duplicates), then <uidstr> won't be in <annotation_dict>
            if ignore_duplicates:
                dup_keys.append(uidstr)
            elif duplicate_resolution_key is None:
                print('  %s multiple annotations for the same cluster, but no duplication resolution key was specified, so returning None for annotation dict (which is fine, as long as you\'re not then trying to use it)' % color('yellow', 'warning'))
                return None
            else:
                uidstr = '%s-%s' % (uidstr, line[duplicate_resolution_key])
                if uidstr in annotation_dict:
                    raise Exception('duplicate key even after adding duplicate resolution key: \'%s\'' % uidstr)
        assert ignore_duplicates or uidstr not in annotation_dict
        annotation_dict[uidstr] = line
    if len(dup_keys) > 0 and not quiet:
        print('    %d duplicate annotations when getting annotation dict (ignore_duplicates was set, so this is probably fine)' % len(dup_keys))

    if cpath is not None:  # check that all the clusters in <cpath>.best() are in the annotation dict
        if len(cpath.best()) > 10000:
            print('     note: to save time, not checking concordance of cpath and annotations since partition has more than 10,000 clusters')
        else:
            check_concordance_of_cpath_and_annotations(cpath, annotation_list, annotation_dict, use_last=use_last)
    return annotation_dict

# ----------------------------------------------------------------------------------------
def get_non_vj_len(line):
    return line['regional_bounds']['j'][0] - line['regional_bounds']['v'][1]

# ----------------------------------------------------------------------------------------
def per_seq_val(line, key, uid, use_default=False, default_val=None):  # get value for per-sequence key <key> corresponding to <uid> NOTE now I've written this, I should really go through and use it in all the places where I do it by hand (for search: iseq)
    if key not in linekeys['per_seq']:
        raise Exception('key \'%s\' not in per-sequence keys' % key)
    if use_default and key not in line or uid not in line['unique_ids']:
        return default_val
    return line[key][line['unique_ids'].index(uid)]  # NOTE just returns the first one, idgaf if there's more than one (and maybe I won't regret that...)

# ----------------------------------------------------------------------------------------
def antnval(antn, key, iseq, use_default=False, default_val=None, add_xtr_col=False):  # generalizes per_seq_val(), and maybe they should be integrated? but adding this long afterwards so don't want to mess with that fcn
    from . import treeutils
    # ----------------------------------------------------------------------------------------
    def rtnval():
        # assert key in linekeys['per_seq']  # input metafile keys (e.g. chosens) are no longer always added to per_seq keys
        if key in linekeys['per_family']:
            return antn[key]
        else:
            if iseq > len(antn[key]) - 1:
                raise Exception('per-seq list too long for %s (uid len %d): %d > %d-1: %s' % (key, len(antn['unique_ids']), iseq, len(antn[key]), antn[key]))
            return antn[key][iseq]
    # ----------------------------------------------------------------------------------------
    # NOTE this is starting to duplicate code in add_extra_column()
    if key == 'aa-cdist':
        key = 'cons-dist-aa'  # arggggggh (when we write output info when choosing abs, we write it as 'aa-cdist', and can't change it now)
    if key == 'aa-cfrac':
        key = 'cons-frac-aa'  # arggggggh (when we write output info when choosing abs, we write it as 'aa-cfrac', and can't change it now)
    if key in antn:
        return rtnval()
    elif 'tree-info' in antn and key in antn['tree-info']['lb']:
        return treeutils.smvals(antn, key, iseq=iseq)
    elif key in ['cons-frac-aa', 'cons-dist-aa']:  # NOTE this doesn't check to see if it's there (i think because we don't store it), it just calculates it
        return treeutils.lb_cons_dist(antn, iseq, aa=True, frac='frac' in key)  # we look to see if it's in the tree info in the previous lines, but if it isn't there, whatever, just recalculate
    elif key in ['shm-aa', 'shm_aa']:  # NOTE this doesn't either, which [also] may be worth fixing eventually
        return shm_aa(antn, iseq=iseq)
    elif key == 'shm':  # NOTE this probably shouldn't work like this, but i just want the h+l selection metrics to work atm
        return antn['n_mutations'][iseq]
    # elif key == 'aa-cdist':
    #     cd = treeutils.smvals(antn, 'cons-dist-aa', iseq=iseq)
    #     return cd if cd is None else -cd
    elif key == 'multipy':  # multiplicity
        return get_multiplicity(antn, iseq=iseq)
    elif key == 'cdr3_seq':
        return get_cdr3_seq(antn, iseq)
    elif key == 'cdr3_seq_aa':
        return ltranslate(get_cdr3_seq(antn, iseq))
    elif key == 'imgt_cdr3_length_aa':  # ick
        return antn['cdr3_length'] // 3 - 2
    elif key == 'frame':
        return get_frame(antn)
    else:
        if add_xtr_col:
            rval = add_extra_column(key, antn, None)  # try to add it NOTE this fcn even existing is pretty hackey, it was originally for transferring from info to outfo 
            if rval is not None:
                antn[key] = rval
        if key in antn:
            return rtnval()
        elif use_default:
            return default_val
        else:
            raise Exception('key \'%s\' not found in line' % key)

# ----------------------------------------------------------------------------------------
def n_dups(line):  # number of duplicate sequences summed over all iseqs (note: duplicates do not include the actual uids/sequences in the annotation) UPDATE: they also don't include 'multiplicities'
    return len([u for dlist in line['duplicates'] for u in dlist])

# ----------------------------------------------------------------------------------------
def uids_and_dups(line):  # NOTE it's kind of weird to have this ignore 'multiplicites' (if it's set), but I'm pretty sure the places this fcn gets used in treeutils only want 'duplicates' included
    return line['unique_ids'] + [u for dlist in line['duplicates'] for u in dlist]

# ----------------------------------------------------------------------------------------
def get_multiplicity(line, uid=None, iseq=None):  # if 'multiplicities' is set (i.e. from input meta info), we return that, since it *should* have been updated in waterer to include duplicates. Otherwise we use the 'duplicates' key NOTE we can't combine them here, since we have no way of being sure whether they have previously been combined, i.e. if 'multiplicities' is there, it *has* to be correct
    if uid is None:
        def ifcn(k): return line[k][iseq]
    elif iseq is None:
        def ifcn(k): return per_seq_val(line, k, uid)
    else:
        assert False
    if 'multiplicities' in line:  # if there was input meta info passed in
        return None if ifcn('multiplicities') == 'None' else ifcn('multiplicities')  # stupid old style csv files not converting correctly (i know i could fix it in the conversion function but i don't want to touch that stupid old code)
    elif 'duplicates' in line:
        return len(ifcn('duplicates')) + 1
    else: # this can probably only happen on very old files (e.g. it didn't used to get added to simulation)
        return 1

# ----------------------------------------------------------------------------------------
def get_multiplicities(line):  # combines duplicates with any input meta info multiplicities (well, the 'multiplicities' key in <line> should already have them combined [see waterer])
    return [get_multiplicity(line, iseq=i) for i in range(len(line['unique_ids']))]

# ----------------------------------------------------------------------------------------
def get_abundance_info(antn_list):
    abdnvals, max_abvals = [], []
    for atn in antn_list:
        tabvals = []  # abundance values for this annotation
        for tseq, sgroup in itertools.groupby(sorted(atn['seqs'])):
            tabvals.append(len(list(sgroup)))
        abdnvals += tabvals
        max_abvals.append(max(tabvals))
    return abdnvals, max_abvals

# ----------------------------------------------------------------------------------------
# NOTE see get_non_implicit_copy() below (if you're thinking of again trying to add it here)

# ----------------------------------------------------------------------------------------
# translate the uids in each line in <antn_list> using translation dict <trns>
# specify *either* <trns> (a dict from old to new uid) or <trfcn> (a fcn from old to new uid)
def translate_uids(antn_list, trns=None, trfcn=None, cpath=None, failstr='translation', no_fail=False, translate_pids=False, expect_missing=False, debug=False):
    from . import treeutils
    # ----------------------------------------------------------------------------------------
    def tr_tree(line, treestr, dbgstr):
        dtree = treeutils.get_dendro_tree(treestr=treestr)
        treeutils.translate_labels(dtree, [(u, trfn(u)) for u in line['unique_ids']], dont_fail=True, expect_missing=expect_missing, dbgstr=dbgstr, debug=debug)
        return treeutils.as_str(dtree)
    # ----------------------------------------------------------------------------------------
    def trfn(uid, pid=False):
        if trfcn is not None:
            tid = trfcn(uid)
            if tid in revrns:
                print('    %s translated uid %s (for %s) already in translations dict for %s' % (wrnstr(), tid, uid, revrns[tid]))
            return tid
        elif uid in trns:
            return trns[uid]
        else:
            if not no_fail and uid in line['unique_ids']:  # we expect to be missing translations from e.g. inferred nodes in the tree, but not from line['unique_ids']
                raise Exception('no %s for %s' % (failstr, uid))  # everybody has to have exactly one paired id at this point
            missing_translations.add(uid)  # this is kind of pointless since it'll crash in other places if things are missing, but oh well
            return uid
    # ----------------------------------------------------------------------------------------
    def tr_lb_info(line):
        lbfo = line['tree-info']['lb']
        if debug:
            print('  translating lb info:')
        for mtmp, mtfo in lbfo.items():
            if mtmp in ['tree', 'aa-tree']:
                lbfo[mtmp] = tr_tree(line, lbfo[mtmp], 'aa inf' if mtmp=='aa-tree' else 'nuc inf')
            else:
                lbfo[mtmp] = {trfn(u) : mtfo[u] for u in mtfo}
    # ----------------------------------------------------------------------------------------
    if 'paired-uids' in antn_list[0] and not translate_pids:
        print('  %s paired uids in line, but <translate_pids> wasn\'t set, so they won\'t be translated' % wrnstr())
    missing_translations = set()
    revrns = {}  # reverse translations
    for line in antn_list:
        if debug:
            print('translating cluster len %d' % len(line['unique_ids']))
        if 'tree' in line:  # simulation
            line['tree'] = tr_tree(line, line['tree'], 'true')
        if 'tree-info' in line and 'lb' in line['tree-info']:
            tr_lb_info(line)
        for iseq, old_id in enumerate(line['unique_ids']):
            new_id = trfn(old_id)
            line['unique_ids'][iseq] = new_id
            revrns[new_id] = old_id
            if 'paired-uids' in line and translate_pids:
                line['paired-uids'][iseq] = [trfn(p, pid=True) for p in line['paired-uids'][iseq]]
    if cpath is not None:
        if cpath.seed_unique_id is not None:
            cpath.seed_unique_id = trfn(cpath.seed_unique_id)
        for iptn, partition in enumerate(cpath.partitions):
            cpath.partitions[iptn] = [[trfn(u) for u in c] for c in partition]
    if len(missing_translations) > 0 and not expect_missing:
        print('    %s missing translations for %d: %s' % (color('yellow', 'warning'), len(missing_translations), ' '.join(missing_translations)))

    return revrns

# ----------------------------------------------------------------------------------------
def choose_non_dup_id(uid, n_duplicate_uids, input_info):  # NOTE calling fcn has to add the new id to <input_info> (which must be some sort of iterable)
    new_uid = uid
    iid = 2
    while new_uid in input_info:
        if any('-'+l in uid for l in loci):  # it's nice to insert <iid> before the locus in simulation uids
            ltmp = [l for l in loci if '-'+l in uid][0]  # would be nice to use get_single_entry() but i really don't want this to ever crash
            new_uid = uid.replace('-'+ltmp, '-%d-%s'%(iid, ltmp))
        else:
            new_uid = uid + '-' + str(iid)
        iid += 1
    if n_duplicate_uids == 0:
        print('  %s duplicate uid(s) in input file, so renaming by appending integer string, e.g. \'%s\' --> \'%s\'' % (color('yellow', 'warning'), uid, new_uid))
    n_duplicate_uids += 1
    return new_uid, n_duplicate_uids  # if you decide you want to change it also in <reco_info>, don't forget to also modify the tree (and maybe other stuff, hence why I don't want to do it)

# ----------------------------------------------------------------------------------------
# NOTE the consensus seqs will (obviously) be *different* afterwards
def synthesize_single_seq_line(line, iseq, dont_deep_copy=False):  # setting dont_deep_copy is obviously *really* *dangerous*
    """ without modifying <line>, make a copy of it corresponding to a single-sequence event with the <iseq>th sequence """
    singlefo = {}
    for key in line:
        if key in linekeys['per_seq']:
            singlefo[key] = [line[key][iseq], ]  # used to also deepcopy the per-seq value, but it's really slow and i really think there's no reason to
        else:
            singlefo[key] = line[key] if dont_deep_copy else copy.deepcopy(line[key])
    return singlefo

# ----------------------------------------------------------------------------------------
# <reco_info> just needs to be an annotation dict with single-sequence annotations
def synthesize_multi_seq_line_from_reco_info(uids, reco_info, warn=False): #, dont_deep_copy=False):  # assumes you already added all the implicit info
    assert len(uids) > 0
    # invalid = False
    # if any(reco_info[u]['invalid'] for u in uids):
    #     print('    %s invalid events in synthesize_multi_seq_line_from_reco_info()' % wrnstr())
    #     invalid = True
    multikeys = set(k for u in uids for k in reco_info[u])  # need to have any key that's in any annotation
    multifo = {}
    for col in sorted(multikeys):
        if col in linekeys['per_seq']:
            cval = [per_seq_val(reco_info[u], col, u, use_default=True, default_val=nullval(col)) for u in uids]
        else:  # would be nice to use linekeys['per_family'], but not all keys are in one or the other (e.g. linekeys['sw'] aren't in either)
            cvals = [reco_info[u][col] for u in uids if col in reco_info[u]]
            if warn and any(v != cvals[0] for v in cvals):  # eh, maybe i don't need to raise an exception here? Although generally you sure shouldn't be merging them if they differ
                print('    %s multiple values for key \'%s\' when synthesizing multi-seq line (just using the first one): %s' % (wrnstr(), col, '')) #cvals))
            cval = cvals[0]
        multifo[col] = copy.deepcopy(cval)
    return multifo

# ----------------------------------------------------------------------------------------
# add seqs in <seqfos_to_add> to the annotation in <line>, aligning new seqs against <line>'s naive seq with mafft if necessary (see bin/add-seqs-to-outputs.py)
# NOTE see also replace_seqs_in_line()
# NOTE also that there's no way to add shm indels (or other additional keys) for seqs in <seqfos_to_add>, you have to replace the info afterwards (see e.g. combine_events())
def add_seqs_to_line(line, seqfos_to_add, glfo, try_to_fix_padding=False, refuse_to_align=False, print_added_str='', extra_str='      ', dont_print=False, debug=False):
    from . import indelutils
    # ----------------------------------------------------------------------------------------
    def align_sfo_seqs(sfos_to_align):
        sfos_to_align['naive_seq'] = line['naive_seq']
        msa_info = align_many_seqs([{'name' : n, 'seq' : s} for n, s in sfos_to_align.items()])
        aligned_naive_seq = get_single_entry([sfo['seq'] for sfo in msa_info if sfo['name'] == 'naive_seq'])
        msa_info = [sfo for sfo in msa_info if sfo['name'] != 'naive_seq']
        if debug:
            print('  aligned %d seq%s with length different to naive sequence:' % (len(sfos_to_align) - 1, plural(len(sfos_to_align) - 1)))
            for iseq, sfo in enumerate(msa_info):
                color_mutants(aligned_naive_seq, sfo['seq'], print_result=True, ref_label='naive seq ', seq_label=sfo['name']+' ', extra_str='        ', only_print_seq=iseq>0)
        for sfo in msa_info:
            trimmed_seq = []  # it could be padded too, but probably it'll mostly be trimmed, so we just call it that
            for naive_nuc, new_nuc in zip(aligned_naive_seq, sfo['seq']):
                if naive_nuc in gap_chars:  # probably extra crap on the ends of the new sequence
                    continue
                elif new_nuc in gap_chars:  # naive seq probably has some N padding
                    trimmed_seq.append(ambig_base)
                else:
                    trimmed_seq.append(new_nuc)
            sfos_to_align[sfo['name']] = ''.join(trimmed_seq)
            assert len(sfos_to_align[sfo['name']]) == len(line['naive_seq'])

    # ----------------------------------------------------------------------------------------
    def getseq(sfo):  # could use recursive .get(, .get()), but i just feel like making it more explicit
        # trimmed_seqfos.get(sfo['name'], sfos_to_align.get(sfo['name'], sfo['seq'])
        if sfo['name'] in trimmed_seqfos:
            return trimmed_seqfos[sfo['name']]
        elif sfo['name'] in sfos_to_align:
            return sfos_to_align[sfo['name']]
        else:
            return sfo['seq']

    # ----------------------------------------------------------------------------------------
    original_len, original_uids = len(line['unique_ids']), copy.copy(line['unique_ids'])
    trimmed_seqfos, sfos_to_align = {}, {}
    if try_to_fix_padding:
        for sfo in [s for s in seqfos_to_add if len(s['seq']) > len(line['naive_seq'])]:
            trimmed_seq = remove_ambiguous_ends(sfo['seq'])
            if len(trimmed_seq) == len(line['naive_seq']):  # presumably the naive seq matches any seqs that are already in <line> (and inserts and deletions and whatnot), so we can probably only really fix it if the new seqs are padded but the naive seq isn't
                trimmed_seqfos[sfo['name']] = trimmed_seq
        if debug:
            print('    trimmed %d seq%s to same length as naive seq' % (len(trimmed_seqfos), plural(len(trimmed_seqfos))))
    if not refuse_to_align:
        sfos_to_align = {sfo['name'] : sfo['seq'] for sfo in seqfos_to_add if sfo['name'] not in trimmed_seqfos and len(sfo['seq']) != len(line['naive_seq'])}  # implicit info adding enforces that the naive seq is the same length as all the seqs
        if len(sfos_to_align) > 0:
            align_sfo_seqs(sfos_to_align)
    aligned_seqfos = [{'name' : sfo['name'], 'seq' : getseq(sfo)} for sfo in seqfos_to_add]  # NOTE needs to be in same order as <seqfos_to_add>

    if not line.get('is_fake_paired', False):
        remove_all_implicit_info(line)

    add_per_seq_keys(line)
    for key in sorted(set(line) & set(linekeys['per_seq'])):
        if key == 'unique_ids':
            line[key] += [s['name'] for s in aligned_seqfos]
        elif key == 'input_seqs' or key == 'seqs':  # i think elsewhere these end up pointing to the same list of string objects, but i think that doesn't matter? UPDATE of fuck yes it does
            assert len(line[key]) in [original_len, original_len + len(aligned_seqfos)]
            if len(line[key]) == original_len:  # if they're pointing at the same list, make sure to only add em once (yes this sucks, and yes they shouldn't be pointing at the same list, but if it does happen this keeps it from being a disaster)
                dup_ids = [s['name'] for s in aligned_seqfos if s['name'] in original_uids]
                if len(dup_ids) > 0:  # not really sure i need this dups stuff (didn't need it for the case for which I implemented it), but leaving in, which may be dangerous, hence the warning
                    print('    %s duplicate ids in aligned seqfos (should be handled, but also shouldn\'t happen, i.e. should fix eventually)' % wrnstr())
                line[key] += [s['seq'] for s in aligned_seqfos if s['name'] not in original_uids]
        elif key == 'duplicates':
            line[key] += [[] for _ in aligned_seqfos]
        elif key == 'multiplicities':
            line[key] += [1 for _ in aligned_seqfos]
        elif key == 'has_shm_indels':
            line[key] += [False for _ in aligned_seqfos]  # uh, I think?
        elif key == 'indelfos':
            line[key] += [indelutils.get_empty_indel() for _ in aligned_seqfos]
        else:  # I think this should only be for input meta keys like multiplicities, affinities, and timepoints, and hopefully they can all handle None?
            line[key] += [None for _ in aligned_seqfos]
        if len(line[key]) != original_len + len(aligned_seqfos):
            print(key, original_len, len(line[key]), original_len + len(aligned_seqfos))
            raise Exception('value for key \'%s\' not correct length when adding seqs to line (see previous line)' % key)

    if line.get('is_fake_paired', False):
        if 'seqs_aa' in line:  # fill in the Nones that will have been added above
            add_seqs_aa(line)
        if 'n_mutations' in line:
            line['n_mutations'] = [hamming_distance(line['naive_seq'], line['seqs'][i]) if n is None else n for i, n in enumerate(line['n_mutations'])]
        if 'mut_freqs' in line:
            line['mut_freqs'] = [hamming_fraction(line['naive_seq'], line['seqs'][i]) if f is None else f for i, f in enumerate(line['mut_freqs'])]
        if 'shm_aa' in line:
            line['shm_aa'] = [shm_aa(line, iseq=i) if n is None else n for i, n in enumerate(line['shm_aa'])]
    else:
        add_implicit_info(glfo, line)

    if print_added_str:
        print('%sadded %d %s seqs to line (originally with %d)%s' % (extra_str, len(seqfos_to_add), print_added_str, original_len, '' if len(seqfos_to_add)>16 else ': %s' % ' '.join(s['name'] for s in seqfos_to_add)))
    if debug and not line.get('is_fake_paired', False) and not dont_print:  # dont_print sucks, but combine_events() needs to be able to not print here (since we haven't yet fixed some keys, e.g. indel info)
        print_reco_event(line, label='after adding %d seq%s:'%(len(aligned_seqfos), plural(len(aligned_seqfos))), extra_str=extra_str, queries_to_emphasize=[s['name'] for s in aligned_seqfos])

# ----------------------------------------------------------------------------------------
# same as add_seqs_to_line(), except this removes all existing seqs first, so the final <line> only contains the seqs in <seqfos_to_add>
# NOTE this removes any existing per-seq info in <line>, e.g. shm indels (of course, since they pertain only to seqs that we're removing)
def replace_seqs_in_line(line, seqfos_to_add, glfo, try_to_fix_padding=False, refuse_to_align=False, debug=False):
    n_seqs_to_remove = len(line['unique_ids'])
    add_seqs_to_line(line, seqfos_to_add, glfo, try_to_fix_padding=try_to_fix_padding, refuse_to_align=refuse_to_align, debug=debug)
    iseqs_to_keep = list(range(n_seqs_to_remove, len(line['unique_ids'])))
    restrict_to_iseqs(line, iseqs_to_keep, glfo)

# ----------------------------------------------------------------------------------------
def combine_events(glfo, evt_list, meta_keys=None, extra_str='      ', debug=False):  # combine events in <evt_list> into a single annotation (could also [used to] pass in meta info values, but don't need it atm)
    from . import indelutils
    ldists, rdists = zip(*[get_pad_parameters(l, glfo) for l in evt_list])  # similar to code in re_pad_hmm_seqs()
    if len(set(ldists)) > 1 or len(set(rdists)) > 1:  # if different annotations have different padding, fix it
        if debug:
            print('            different padding parameters (left dists %s  right dists %s) when combining events, so padding out to max on each side' % (ldists, rdists))
        for ia, tatn in enumerate(evt_list):  # pad all the annotations to the max padding (on each side) from any annotation
            leftpad, rightpad = max(ldists) - ldists[ia], max(rdists) - rdists[ia]
            re_pad_atn(leftpad, rightpad, tatn, glfo, extra_str='              ') #, debug=debug)
    combo_evt = get_full_copy(evt_list[0], glfo)
    if debug > 1:
        print('%scombining %d annotations:' % (extra_str, len(evt_list)))
        for tevt in evt_list:
            print_reco_event(tevt, label=color('yellow', 'before:'), extra_str=extra_str)
    seqfos_to_add, added_atn_indices, added_uids = [], [], []  # handle uids that appear in more than one event
    for ia, atn in enumerate(evt_list[1:]):
        for uid, seq in zip(atn['unique_ids'], atn['seqs']):
            if uid in combo_evt['unique_ids'] or uid in added_uids:
                continue
            seqfos_to_add.append({'name' : uid, 'seq' : seq})
            added_atn_indices.append(ia + 1)  # NOTE <ia> is first *added* event, i.e. not combo_evt (hence +1)
            added_uids.append(uid)
    # NOTE add_seqs_to_line() just adds *seqs* while assuming they have no associated extra info (e.g. indel info), which is wrong here, so we have to fix (replace the default info that is added) afterward (in the next lines)
    add_seqs_to_line(combo_evt, seqfos_to_add, glfo, extra_str=extra_str, dont_print=True, debug=debug)
    for tkey in set(combo_evt) & set(linekeys['per_seq']) - set(implicit_linekeys) - set(['unique_ids', 'seqs']):  # TODO this doesn't account for if we align stuff in add_seqs_to_line()
        combo_evt[tkey] = evt_list[0][tkey] + [per_seq_val(evt_list[i], tkey, u) for i, u in zip(added_atn_indices, added_uids)]  # NOTE i here is index in <evt_list> (see note above)
    if meta_keys is not None:
        for mkey in meta_keys:
            combo_evt[mkey] = [v for l in evt_list for v in l[mkey]]
    if 'tree' in combo_evt:  # would not be easy to merge the trees from different events, so give up for now
        combo_evt['tree'] = None
    if debug > 1:
        print_reco_event(combo_evt, label='%s adding %d seq%s:'%(color('green', 'after'), len(seqfos_to_add), plural(len(seqfos_to_add))), extra_str=extra_str, queries_to_emphasize=[s['name'] for s in seqfos_to_add])
    return combo_evt

# ----------------------------------------------------------------------------------------
def get_repfracstr(csize, repertoire_size):  # return a concise string representing <csize> / <repertoire_size>
    if csize > repertoire_size:
        return '1.'  # this seems to mean the repertoire is just one cluster (maybe with one sequence?) and i can't be bothered to fix it
    repfrac = float(csize) / repertoire_size
    denom = int(1. / repfrac)
    estimate = 1. / denom
    frac_error = (estimate - repfrac) / repfrac
    if frac_error > 0.10:  # if it's more than 10% off just use the float str
        # print 'pretty far off: (1/denom - repfrac) / repfrac = (1./%d - %f) / %f = %f' % (denom, repfrac, repfrac, frac_error)
        repfracstr = '%.2f' % repfrac
    elif denom > 1000:
        repfracstr = '%.0e' % repfrac
    else:
        repfracstr = '1/%d' % denom
    return repfracstr

# ----------------------------------------------------------------------------------------
def generate_dummy_v(d_gene):
    pv, sv, al = split_gene(d_gene)
    return get_locus(d_gene).upper() + 'VxDx' + pv + '-' + sv + '*' + al

# ----------------------------------------------------------------------------------------
def get_hash(instr):
    return hashlib.md5(instr.encode()).hexdigest()  # , usedforsecurity=False (can't use this arg til python 3.9)

# ----------------------------------------------------------------------------------------
def uidhashstr(instr, max_len=10):
    return get_hash(instr)[:max_len]

# ----------------------------------------------------------------------------------------
def hashint(instr, max_int=2**31):  # max random seed val is 32, so reduce by 1 to get some leeway
    return int(get_hash(instr), base=16) % max_int

# ----------------------------------------------------------------------------------------
# NOTE see seqfileopener.py or treeutils.py for example usage (both args should be set to None the first time through)
def choose_new_uid(potential_names, used_names, initial_length=1, n_initial_names=None, available_chars=string.ascii_lowercase, repeat_chars=False, shuffle=False, dont_extend=False):
    # NOTE only need to set <initial_length> for the first call -- after that if you're reusing the same <potential_names> and <used_names> there's no need (but it's ok to set it every time, as long as it has the same value)
    # NOTE setting <shuffle> will shuffle every time, i.e. it's designed such that you call with shuffle *once* before starting
    def get_potential_names(length):
        itobj = itertools.product(available_chars, repeat=length) if repeat_chars else itertools.combinations(available_chars, length)
        if n_initial_names is None:
            return [''.join(ab) for i, ab in enumerate(itobj)]
        else:
            rstrs = []
            for icombo, chars in enumerate(itobj):
                rstrs.append(''.join(chars))
                if len(rstrs) >= n_initial_names:
                    break
            return rstrs
    if potential_names is None:  # first time through
        potential_names = get_potential_names(initial_length)
        if n_initial_names is not None and len(potential_names) != n_initial_names:
            raise Exception('couldn\'t make enough potential names (asked for %d, made %d) with length %s from chars %s' % (n_initial_names, len(potential_names), initial_length, available_chars))
        used_names = []
    if len(potential_names) == 0:  # ran out of names
        if dont_extend:
            raise Exception('ran out of names, but <dont_extend> was set')
        potential_names = get_potential_names(len(used_names[-1]) + 1)
    if len(potential_names[0]) < initial_length:
        raise Exception('choose_new_uid(): next potential name \'%s\' is shorter than the specified <initial_length> %d (this is probably only possible if you called this several times with different <initial_length> values [which you shouldn\'t do])' % (potential_names[0], initial_length))
    if shuffle:
        random.shuffle(potential_names)
    new_id = potential_names.pop(0)
    used_names.append(new_id)
    return new_id, potential_names, used_names

# ----------------------------------------------------------------------------------------
def convert_from_adaptive_headers(glfo, line, uid=None, only_dj_rearrangements=False):
    from . import indelutils
    from . import glutils
    newline = {}
    print_it = False

    for head, ahead in adaptive_headers.items():
        newline[head] = line[ahead]
        if head in io_column_configs['lists']:
            newline[head] = [newline[head], ]

    if uid is not None:
        newline['unique_ids'] = [uid, ]

    for erosion in real_erosions:
        newline[erosion + '_del'] = int(newline[erosion + '_del'])
    newline['v_5p_del'] = 0
    newline['j_3p_del'] = 0

    for region in regions:
        if newline[region + '_gene'] == 'unresolved':
            newline[region + '_gene'] = None
            continue
        if region == 'j' and 'P' in newline[region + '_gene']:
            newline[region + '_gene'] = newline[region + '_gene'].replace('P', '')

        if '*' not in newline[region + '_gene']:
            # print uid
            # tmpheads = ['dMaxResolved', 'dFamilyName', 'dGeneName', 'dGeneAllele', 'dFamilyTies', 'dGeneNameTies', 'dGeneAlleleTies']
            # for h in tmpheads:
            #     print '   %s: %s' % (h, line[h]),
            # print ''
            if line['dGeneAlleleTies'] == '':
                newline['failed'] = True
                return newline
            d_alleles = line['dGeneAlleleTies'].split(',')
            newline[region + '_gene'] += '*' + d_alleles[0]
        primary_version, sub_version, allele = split_gene(newline[region + '_gene'])
        primary_version, sub_version = primary_version.lstrip('0'), sub_version.lstrip('0')  # alleles get to keep their leading zero (thank you imgt for being consistent)
        if region == 'j':  # adaptive calls every j sub_version 1
            sub_version = None
        gene = rejoin_gene(glfo['locus'], region, primary_version, sub_version, allele)
        if gene not in glfo['seqs'][region]:
            gene = glutils.convert_to_duplicate_name(glfo, gene)
        if gene not in glfo['seqs'][region]:
            raise Exception('couldn\'t rebuild gene name from adaptive data: %s' % gene)
        newline[region + '_gene'] = gene

    seq = newline['seqs'][0]
    boundlist = ['vIndex', 'n1Index', 'dIndex', 'n2Index', 'jIndex']
    qrbounds, glbounds = {}, {}
    for region in regions:
        if only_dj_rearrangements and region == 'v':
            if newline['d_gene'] is None:
                newline['failed'] = True
                return newline
            newline['v_gene'] = generate_dummy_v(newline['d_gene'])
            line['vIndex'] = 0
            line['n1Index'] = int(line['dIndex']) - 1  # or is it without the -1?
            glfo['seqs']['v'][newline['v_gene']] = seq[line['vIndex'] : line['n1Index']]
            glfo['cyst-positions'][newline['v_gene']] = len(glfo['seqs']['v'][newline['v_gene']]) - 3
        if newline[region + '_gene'] is None:
            newline['failed'] = True
            return newline
        qrb = [int(line[region + 'Index']),
                    int(line[boundlist[boundlist.index(region + 'Index') + 1]]) if region != 'j' else len(seq)]
        glseq = glfo['seqs'][region][newline[region + '_gene']]
        glb = [newline[region + '_5p_del'],
                    len(glseq) - newline[region + '_3p_del']]
        if region == 'j' and glb[1] - glb[0] > qrb[1] - qrb[0]:  # extra adaptive stuff on right side of j
            old = glb[1]
            glb[1] = glb[0] + qrb[1] - qrb[0]
            newline['j_3p_del'] = old - glb[1]

        if qrb[0] == -1 or qrb[1] == -1 or qrb[1] < qrb[0]:  # should this also be equals?
            newline['failed'] = True
            return newline
        if qrb[1] - qrb[0] != glb[1] - glb[0]:
            newline['failed'] = True
            return newline
        qrbounds[region] = qrb
        glbounds[region] = glb

    for bound in boundaries:
        newline[bound + '_insertion'] = seq[qrbounds[bound[0]][1] : qrbounds[bound[1]][0]]  # end of lefthand region to start of righthand region

    newline['fv_insertion'] = ''
    newline['jf_insertion'] = seq[qrbounds['j'][1]:]

    # print seq
    # print seq[:qrbounds['d'][0]],
    # print seq[qrbounds['d'][0] : qrbounds['d'][1]],
    # print seq[qrbounds['d'][1] : qrbounds['j'][0]],
    # print seq[qrbounds['j'][0] : qrbounds['j'][1]],
    # print seq[qrbounds['j'][1] :]

    newline['indelfos'] = [indelutils.get_empty_indel(), ]

    if print_it:
        add_implicit_info(glfo, newline)
        print_reco_event(newline, label=uid)

    # still need to convert to integers/lists/whatnot (?)

    newline['failed'] = False

    return newline

# ----------------------------------------------------------------------------------------
# definitions here: http://clip.med.yale.edu/changeo/manuals/Change-O_Data_Format.pdf
presto_headers = OrderedDict([  # enforce this ordering so the output files are easier to read
    ('SEQUENCE_ID', 'unique_ids'),
    ('V_CALL', 'v_gene'),
    ('D_CALL', 'd_gene'),
    ('J_CALL', 'j_gene'),
    ('JUNCTION_LENGTH', None),
    ('SEQUENCE_INPUT', 'input_seqs'),
    ('SEQUENCE_IMGT', 'aligned_v_plus_unaligned_dj'),
])

# reference: https://docs.airr-community.org/en/stable/datarep/rearrangements.html#fields
airr_headers = OrderedDict([  # enforce this ordering so the output files are easier to read
    # required:
    ('sequence_id', 'unique_ids'),
    ('sequence', 'input_seqs'),
    ('rev_comp', None),
    ('productive', None),
    ('v_call', 'v_gene'),
    ('d_call', 'd_gene'),
    ('j_call', 'j_gene'),
    ('sequence_alignment', None),
    ('germline_alignment', None),
    ('junction', None),  # NOTE this is *actually* the junction, whereas what partis calls the cdr3 is also actually the junction (which is terrible, but i swear it's not entirely my fault, but either way it's just too hard to change now)
    ('junction_aa', None),
    ('v_cigar', None),
    ('d_cigar', None),
    ('j_cigar', None),
    ('clone_id', None),
    # optional:
    ('vj_in_frame', 'in_frames'),
    ('stop_codon', 'stops'),
    ('locus', None),
    ('np1', 'vd_insertion'),
    ('np2', 'dj_insertion'),
    ('duplicate_count', None),
    ('cdr3_start', None),
    ('cdr3_end', None),
    ('cell_id', None),
    ('v_sequence_alignment', None),
    ('d_sequence_alignment', None),
    ('j_sequence_alignment', None),
    ('v_germline_alignment', None),
    ('d_germline_alignment', None),
    ('j_germline_alignment', None),
])
for rtmp in regions:
    airr_headers[rtmp+'_support'] = None      # NOTE not really anywhere to put the alternative annotation, which is independent of this and maybe more accurate
    airr_headers[rtmp+'_identity'] = None
    airr_headers[rtmp+'_sequence_start'] = None
    airr_headers[rtmp+'_sequence_end'] = None
    airr_headers[rtmp+'_germline_start'] = None
    airr_headers[rtmp+'_germline_end'] = None

invalid_airr_heads = OrderedDict([
    ('sequence_id', 'unique_ids'),
    ('sequence', 'input_seqs'),
])

linearham_headers = OrderedDict((
    ('Iteration', None),
    ('RBLogLikelihood', None),
    ('Prior', None),
    ('alpha', None),
    ('er[1]', None), ('er[2]', None), ('er[3]', None), ('er[4]', None), ('er[5]', None), ('er[6]', None),
    ('pi[1]', None), ('pi[2]', None), ('pi[3]', None), ('pi[4]', None),
    ('tree', None),
    ('sr[1]', None), ('sr[2]', None), ('sr[3]', None), ('sr[4]', None),
    ('LHLogLikelihood', None),
    ('LogWeight', None),
    ('NaiveSequence', 'naive_seq'),
    ('VGene', 'v_gene'),
    ('V5pDel', 'v_5p_del'),
    ('V3pDel', 'v_3p_del'),
    ('VFwkInsertion', 'fv_insertion'),
    ('VDInsertion', 'vd_insertion'),
    ('DGene', 'd_gene'),
    ('D5pDel', 'd_5p_del'),
    ('D3pDel', 'd_3p_del'),
    ('DJInsertion', 'dj_insertion'),
    ('JGene', 'j_gene'),
    ('J5pDel', 'j_5p_del'),
    ('J3pDel', 'j_3p_del'),
    ('JFwkInsertion', 'jf_insertion'),
    ('VJInsertion', 'dj_insertion'),  # NOTE this has to get handled specially
))


# ----------------------------------------------------------------------------------------
def get_line_with_presto_headers(line):  # NOTE doesn't deep copy
    """ convert <line> to presto csv format """
    if len(line['unique_ids']) > 1:  # has to happen *before* utils.get_line_for_output()  UPDATE wtf does this mean?
        raise Exception('multiple seqs not handled for presto output')

    presto_line = {}
    for phead, head in presto_headers.items():
        if head == 'aligned_v_plus_unaligned_dj':
            presto_line[phead] = line['aligned_v_seqs'][0] + line['vd_insertion'] + line['d_qr_seqs'][0] + line['dj_insertion'] + line['j_qr_seqs'][0]
        elif phead == 'JUNCTION_LENGTH':
            presto_line[phead] = line['cdr3_length']  # + 6  oops, no +6... what I call cdr3 length is properly junction length (but would be a colossal clusterfuck to change)
        elif head == 'unique_ids' or head == 'input_seqs':
            presto_line[phead] = line[head][0]
        else:
            presto_line[phead] = line[head]

    return presto_line

# ----------------------------------------------------------------------------------------
def write_presto_annotations(outfname, annotation_list, failed_queries):
    print('   writing presto annotations to %s' % outfname)
    assert getsuffix(outfname) == '.tsv'  # already checked in processargs.py
    with open(outfname, csv_wmode()) as outfile:
        writer = csv.DictWriter(outfile, list(presto_headers.keys()), delimiter=str('\t'))
        writer.writeheader()

        for line in annotation_list:
            if len(line['unique_ids']) == 1:
                writer.writerow(get_line_with_presto_headers(line))
            else:
                for iseq in range(len(line['unique_ids'])):
                    writer.writerow(get_line_with_presto_headers(synthesize_single_seq_line(line, iseq)))

        # and write empty lines for seqs that failed either in sw or the hmm
        if failed_queries is not None:
            for failfo in failed_queries:
                assert len(failfo['unique_ids']) == 1
                writer.writerow({'SEQUENCE_ID' : failfo['unique_ids'][0], 'SEQUENCE_INPUT' : failfo['input_seqs'][0]})

# ----------------------------------------------------------------------------------------
def get_airr_cigar_str(line, iseq, region, qr_gap_seq, gl_gap_seq, debug=False):
    from . import indelutils
    if debug:
        if region == 'v':
            print(line['unique_ids'][iseq])
        print('  ', region)
    assert len(qr_gap_seq) == len(gl_gap_seq)  # this should be checked in a bunch of other places, but it's nice to see it here
    istart, istop = line['regional_bounds'][region]  # start/stop indices in indel-reversed seq
    if indelutils.has_indels(line['indelfos'][iseq]):  # convert to indices in gap seq
        istart += count_gap_chars(gl_gap_seq, unaligned_pos=istart)  # note: v_5p and j_3p dels appear as dots/gaps in the debug printing below, but in the gap seqs here they're still Ns
        istop += count_gap_chars(gl_gap_seq, unaligned_pos=istop)
    assert istop <= len(qr_gap_seq)  # in theory, i fixed this and it can't happen any more...
    regional_qr_gap_seq = istart * gap_chars[0] + qr_gap_seq[istart : istop] + (len(qr_gap_seq) - istop) * gap_chars[0]
    regional_gl_gap_seq = istart * gap_chars[0] + gl_gap_seq[istart : istop] + (len(qr_gap_seq) - istop) * gap_chars[0]
    if debug:
        print('      ', regional_qr_gap_seq)
        print('      ', regional_gl_gap_seq)
    cigarstr = indelutils.get_cigarstr_from_gap_seqs(regional_qr_gap_seq, regional_gl_gap_seq, debug=debug)
    return cigarstr

# ----------------------------------------------------------------------------------------
def get_airr_line(pline, iseq, extra_columns=None, skip_columns=None, args=None, glfo=None, debug=False):
    from . import indelutils
    # ----------------------------------------------------------------------------------------
    def getrgn(tk):  # get region from key name
        rgn = tk.split('_')[0]
        assert rgn in regions
        return rgn
    # ----------------------------------------------------------------------------------------
    if pline['invalid']:
        return {akey : pline[pkey][iseq] if pkey in linekeys['per_seq'] else pline[pkey] for akey, pkey in invalid_airr_heads.items()}

    qr_gap_seq = pline['seqs'][iseq]
    gl_gap_seq = pline['naive_seq']
    if 'indelfos' not in pline:  # I'm not sure exactly which cases result in <pline> missing some keys (maybe read directly from file?) but it fixes it to add implicit info, so whatevs
        add_implicit_info(glfo, pline)
    if indelutils.has_indels(pline['indelfos'][iseq]):
        qr_gap_seq = pline['indelfos'][iseq]['qr_gap_seq']
        gl_gap_seq = pline['indelfos'][iseq]['gl_gap_seq']

    aline = {}
    for akey, pkey in airr_headers.items():
        if skip_columns is not None and akey in skip_columns:
            continue
        if pkey is not None:  # if there's a direct correspondence to a partis key
            aline[akey] = pline[pkey][iseq] if pkey in linekeys['per_seq'] else pline[pkey]
        elif akey == 'rev_comp':
            aline[akey] = False
        elif '_cigar' in akey and akey[0] in regions:
            aline[akey] = get_airr_cigar_str(pline, iseq, akey[0], qr_gap_seq, gl_gap_seq, debug=debug)
        elif akey == 'productive':
            aline[akey] = is_functional(pline, iseq)
        elif akey == 'sequence_alignment':
            aline[akey] = qr_gap_seq
        elif akey == 'germline_alignment':
            aline[akey] = gl_gap_seq
        elif any(akey == r+'_sequence_alignment' for r in regions):
            gap_bounds = indelutils.get_regional_bounds_with_indels_reinstated(pline, iseq)[getrgn(akey)]  # pline['regional_bounds'][getrgn(akey)]
            aline[akey] = qr_gap_seq[gap_bounds[0] : gap_bounds[1]]
        elif any(akey == r+'_germline_alignment' for r in regions):
            gap_bounds = indelutils.get_regional_bounds_with_indels_reinstated(pline, iseq)[getrgn(akey)]  # pline['regional_bounds'][getrgn(akey)]
            aline[akey] = gl_gap_seq[gap_bounds[0] : gap_bounds[1]]
        elif akey == 'junction':
            aline[akey] = get_cdr3_seq(pline, iseq)
        elif akey == 'junction_aa':
            aline[akey] = ltranslate(aline.get('junction', get_cdr3_seq(pline, iseq)))  # should already be in there, since we're using an ordered dict and the previous elif block should've added it
        elif akey == 'clone_id':
            aline[akey] = get_clone_id(pline['unique_ids'])
        elif akey == 'locus':
            aline[akey] = get_locus(pline['v_gene']).upper()  # scoper at least requires upper case
        elif '_support' in akey and akey[0] in regions:  # NOTE not really anywhere to put the alternative annotation, which is independent of this and maybe more accurate
            pkey = akey[0] + '_per_gene_support'
            gcall = pline[akey[0] + '_gene']
            if pkey not in pline or gcall not in pline[pkey]:
                continue
            aline[akey] = pline[pkey][gcall]
        elif akey == 'duplicate_count':
            aline[akey] = get_multiplicity(pline, iseq=iseq)
        elif '_identity' in akey:
            aline[akey] = 1. - get_mutation_rate(pline, iseq, restrict_to_region=getrgn(akey))
        elif any(akey == r+'_sequence_start' for r in regions):
            aline[akey] = pline['regional_bounds'][getrgn(akey)][0] + 1  # +1 to switch to 1-based indexing
        elif any(akey == r+'_sequence_end' for r in regions):
            aline[akey] = pline['regional_bounds'][getrgn(akey)][1]  # +1 to switch to 1-based indexing, -1 to switch to closed intervals, so net zero
        elif any(akey == r+'_germline_start' for r in regions):
            aline[akey] = pline[getrgn(akey)+'_5p_del'] + 1  # +1 to switch to 1-based indexing
        elif any(akey == r+'_germline_end' for r in regions):
            aline[akey] = len(pline[getrgn(akey)+'_gl_seq']) + pline[getrgn(akey)+'_5p_del']  # +1 to switch to 1-based indexing, -1 to switch to closed intervals, so net zero
        elif akey == 'cdr3_start':  # airr uses the imgt (correct) cdr3 definition, which excludes both conserved codons, so we add 3 (then add 1 to switch to 1-based indexing)
            aline[akey] = pline['codon_positions']['v'] + 3 + 1
        elif akey == 'cdr3_end':
            aline[akey] = pline['codon_positions']['j']
        elif akey == 'cell_id':
            did_seps, did_indices = (None, None) if args is None else (args.droplet_id_separators, args.droplet_id_indices)
            aline[akey] = get_droplet_id(pline['unique_ids'][iseq], did_seps, did_indices)
        else:
            raise Exception('unhandled airr key / partis key \'%s\' / \'%s\'' % (akey, pkey))

    if extra_columns is not None:
        for key in extra_columns:
            if key in pline:
                aline[key] = pline[key][iseq] if key in linekeys['per_seq'] else pline[key]
            else:
                add_extra_column(key, pline, aline)
                if key in linekeys['per_seq']:  # have to restrict to iseq
                    aline[key] = aline[key][iseq]

    return aline

# ----------------------------------------------------------------------------------------
def convert_airr_line(aline, glfo):
    from . import indelutils
    from . import glutils
    pline = {}
    for aky, pky in airr_headers.items():
        if pky is not None:  # if there's a direct correspondence to a partis key
            pline[pky] = [aline[aky]] if pky in linekeys['per_seq'] else aline[aky]  # NOTE/TODO all end up as single-sequence annotations
    pline['duplicates'] = [[]]

    # print aline['v_call'], aline['d_call'], aline['j_call']
    if aline['d_germline_start'] == '':  # igblast leaves this blank for light chain (and sometimes for heavy)
        d_start, d_end = [int(aline['v_germline_end']) + 1 for _ in range(2)]
        if has_d_gene(glfo['locus']):
            pline['d_gene'] = list(glfo['seqs']['d'])[0]
            d_end += 1
        else:
            pline['d_gene'] = glutils.dummy_d_genes[glfo['locus']]
        aline['d_germline_start'] = d_start
        aline['d_germline_end'] = d_end

    for rgn in regions:
        if ',' in pline[rgn+'_gene']:  # wtf they just put multiple genes, separated by commas
            pline[rgn+'_gene'] = pline[rgn+'_gene'].split(',')[0]

        pline[rgn+'_5p_del'] = int(aline[rgn+'_germline_start']) - 1
        pline[rgn+'_3p_del'] = len(gseq(glfo, pline[rgn+'_gene'])) - int(aline[rgn+'_germline_end'])
        cigars = indelutils.split_cigarstrs(aline[rgn+'_cigar'])
        if rgn == 'v':
            ctype, clen = cigars[0]
            pline['fv_insertion'] = ambig_base * (clen if ctype=='I' else 0)
        elif rgn == 'j':
            ctype, clen = cigars[-1]
            pline['jf_insertion'] = ambig_base * (clen if ctype=='I' else 0)
    assert len(aline['sequence_alignment']) == len(aline['germline_alignment'])
    ir_seq = []  # need to calculate indel reversed seq (NOTE this duplicates code, i think somewhere in waterer)
    for qch, gch in zip(aline['sequence_alignment'], aline['germline_alignment']):
        if gch in gap_chars:  # shm insertion
            continue
        elif qch in gap_chars:  # shm deletion
            ir_seq.append(gch)
        else:
            ir_seq.append(qch)
    pline['seqs'] = [''.join(ir_seq)]
    pline['qr_gap_seqs'] = [aline['sequence_alignment']]
    pline['gl_gap_seqs'] = [aline['germline_alignment']]
    try:
        add_implicit_info(glfo, pline)  # NOTE can't set check_line_keys=True since it crashes bc we add 'indelfos'
    except:
        elines = traceback.format_exception(*sys.exc_info())
        print(pad_lines(''.join(elines)))
        print('      convert_airr_line(): implicit info adding failed for %s (see above)' % pline['unique_ids'])
        pline['invalid'] = True

    return pline

# ----------------------------------------------------------------------------------------
def write_airr_output(outfname, annotation_list, cpath=None, failed_queries=None, extra_columns=None, skip_columns=None, args=None, glfo=None, debug=False):  # NOTE similarity to add_regional_alignments() (but I think i don't want to combine them, since add_regional_alignments() is for imgt-gapped aligments, whereas airr format doesn't require imgt gaps, and we really don't want to deal with imgt gaps if we don't need to)
    if extra_columns is None:
        extra_columns = []
    print('   writing airr annotations to %s' % outfname)
    assert getsuffix(outfname) == '.tsv'  # already checked in processargs.py
    with open(outfname, csv_wmode()) as outfile:
        oheads = list(airr_headers.keys()) + extra_columns
        if skip_columns is not None:
            oheads = [h for h in oheads if h not in skip_columns]
        writer = csv.DictWriter(outfile, oheads, delimiter=str('\t'))
        writer.writeheader()
        if len(annotation_list) == 0 and cpath is not None:
            print('    writing partition with no annotations')
            for cluster in cpath.best():
                clone_id = get_clone_id(cluster)
                for uid in cluster:
                    writer.writerow({'sequence_id' : uid, 'clone_id' : clone_id})
        for line in annotation_list:
            for iseq in range(len(line['unique_ids'])):
                aline = get_airr_line(line, iseq, extra_columns=extra_columns, skip_columns=skip_columns, args=args, glfo=glfo, debug=debug)
                writer.writerow(aline)

        # and write empty lines for seqs that failed either in sw or the hmm
        if failed_queries is not None:
            for failfo in failed_queries:
                assert len(failfo['unique_ids']) == 1
                writer.writerow({'sequence_id' : failfo['unique_ids'][0], 'sequence' : failfo['input_seqs'][0]})

# ----------------------------------------------------------------------------------------
def read_airr_output(fname, glfo=None, locus=None, glfo_dir=None, skip_other_locus=False, clone_id_field='clone_id', sequence_id_field='sequence_id', delimiter='\t', skip_annotations=False):
    from . import clusterpath
    from . import glutils
    if glfo is None and glfo_dir is not None:
        glfo = glutils.read_glfo(glfo_dir, locus)  # TODO this isn't right
    failed_queries, clone_ids, plines, other_locus_ids = [], {}, [], []
    with open(fname) as afile:
        reader = csv.DictReader(afile, delimiter=str(delimiter))
        for aline in reader:
            if clone_id_field in reader.fieldnames:
                # print '  note: no clone ids in airr file %s' % fname
                clone_ids[aline[sequence_id_field]] = aline[clone_id_field]
            if skip_annotations:
                continue
            if 'sequence' not in aline:
                continue
            if aline['v_call'] == '' or aline['j_call'] == '':
                failed_queries.append({'unique_ids' : [aline[sequence_id_field]], 'input_seqs' : [aline['sequence']], 'invalid' : True})
                continue
            if skip_other_locus and get_locus(aline['v_call']) != glfo['locus']:
                other_locus_ids.append(aline[sequence_id_field])
                continue
            plines.append(convert_airr_line(aline, glfo))
    if len(clone_ids) > 0:
        partition = group_seqs_by_value(list(clone_ids.keys()), lambda q: clone_ids[q])
    else:
        partition = [[l['unique_ids'][0]] for l in plines]
    if skip_other_locus:
        partition = [[u for u in c if u not in other_locus_ids] for c in partition]
        partition = [c for c in partition if len(c) > 0]
    if len(plines) > 0:
        sorted_ids = [l['unique_ids'][0] for l in plines]
        partition = sorted(partition, key=lambda c: min(sorted_ids.index(u) for u in c))  # sort by min index in <sorted_ids> of any uid in each cluster
    antn_list = []
    if len(plines) > 0:
        antn_dict = get_annotation_dict(plines)
        for iclust, cluster in enumerate(partition):  # may as well sort by length, otherwise order is just random
            cluster = [u for u in sorted_ids if u in cluster]  # it's nice to try to keep them in the same order, and if partis wrote the single-seq lines this'll put them back in the same order
            partition[iclust] = cluster
            multi_line = synthesize_multi_seq_line_from_reco_info(cluster, antn_dict)
            # print_reco_event(multi_line, extra_str='  ')
            antn_list.append(multi_line)
    for failfo in failed_queries:
        antn_list.append(failfo)
    return glfo, antn_list, clusterpath.ClusterPath(partition=partition)

# ----------------------------------------------------------------------------------------
def fix_linearham_insertion(lh_line, partis_line):
    v_end = partis_line['regional_bounds']['v'][1]
    j_start = partis_line['regional_bounds']['j'][0]
    for dln in ['v_3p_del', 'j_5p_del']:
        xtra_del = int(lh_line[dln]) - partis_line[dln]
        if xtra_del > 0:  # if linearham made the deletion larger
            if 'v_3p' in dln:
                v_end -= xtra_del
            elif 'j_5p' in dln:
                j_start += xtra_del
            else:
                assert False
    if j_start > v_end:
        seq = partis_line['seqs'][0]
        print('    %s read boolean for VJInsertion in linearham output file, so resetting light chain insertion by hand: \'%s\' --> \'%s\'' % (color('yellow', 'warning'), lh_line['dj_insertion'], seq[v_end : j_start]))
        print('      %s' % (seq[ : v_end] + color('red', seq[v_end : j_start]) + seq[j_start : ]))
        # print_reco_event(partis_line)
        lh_line['dj_insertion'] = seq[v_end : j_start]

# ----------------------------------------------------------------------------------------
def process_input_linearham_line(lh_line, partis_line, glfo):
    """ convert <lh_line> (see linearham_headers). Modifies <lh_line>. """
    for lhk in set(lh_line): # set() used because keys are removed from the dict while iterating
        if lhk not in linearham_headers or linearham_headers[lhk] is None: # limit to the ones with a direct partis correspondence
            del lh_line[lhk] #remove keys not in linearham_headers
            continue
        lh_line[linearham_headers[lhk]] = lh_line[lhk]
        del lh_line[lhk] #remove lh_line keys once corresponding linearham_headers key added
    if not has_d_gene(glfo['locus']) and lh_line['dj_insertion'] == 'TRUE':  # not sure why it writes a boolean for light chain insertion rather than the actual inserted str, but we have to go and work it out since otherwise it's an invalid event
        fix_linearham_insertion(lh_line, partis_line)
    process_input_line(lh_line)

# ----------------------------------------------------------------------------------------
def get_parameter_fname(column=None, deps=None, column_and_deps=None):
    """ return the file name in which we store the information for <column>. Either pass in <column> and <deps> *or* <column_and_deps> """
    if column == 'all':
        return 'all-probs.csv'
    if column_and_deps is None:
        if column is None or deps is None:
            raise Exception('you have to either pass <column_and_deps>, or else pass both <column> and <deps>')
        column_and_deps = [column]
        column_and_deps.extend(deps)
    outfname = 'probs.csv'
    for ic in column_and_deps:
        outfname = ic + '-' + outfname
    return outfname

# ----------------------------------------------------------------------------------------
def from_same_event(reco_info, query_names):  # putting are_clonal in a comment, since that's what I always seem to search for when I'm trying to remember where this fcn is
    if len(query_names) > 1:
        # reco_id = reco_info[query_names[0]]['reco_id']  # the first one's reco id
        # for iq in range(1, len(query_names)):  # then loop through the rest of 'em to see if they're all the same
        #     if reco_id != reco_info[query_names[iq]]['reco_id']:
        #         return False
        # return True

        # darn it, this isn't any faster
        reco_ids = [reco_info[query]['reco_id'] for query in query_names]
        return reco_ids.count(reco_ids[0]) == len(reco_ids)  # they're clonal if every reco id is the same as the first one

    else:  # one or zero sequences
        return True

# ----------------------------------------------------------------------------------------
# bash color codes
ansi_color_table = collections.OrderedDict((
    # 'head' : '95'  # not sure wtf this was?
    ('end', 0),
    ('bold', 1),
    ('reverse_video', 7),
    ('grey', 90),
    ('red', 91),
    ('green', 92),
    ('yellow', 93),
    ('blue', 94),
    ('purple', 95),
    ('grey_bkg', 100),
    ('red_bkg', 41),
    ('green_bkg', 42),
    ('yellow_bkg', 43),
    ('blue_bkg', 44),
    ('light_blue_bkg', 104),
    ('purple_bkg', 45),
))
Colors = {c : '\033[%sm'%i for c, i in ansi_color_table.items()}

def color(col, seq, width=None, padside='left'):
    return_str = [seq]
    if col is not None:
        return_str = [Colors[col]] + return_str + [Colors['end']]
    if width is not None:  # make sure final string prints to correct width
        n_spaces = max(0, width - len(seq))  # if specified <width> is greater than uncolored length of <seq>, pad with spaces so that when the colors show up properly the colored sequences prints with width <width>
        if padside == 'left':
            return_str.insert(0, n_spaces * ' ')
        elif padside == 'right':
            return_str.insert(len(return_str), n_spaces * ' ')
        else:
            assert False
    return ''.join(return_str)

def cyclecolor(cid, clist=None):  # cycle through the (actual) colors in <Colors> (note: no control/warning over when it's wrapping around)
    if clist is None:
        clist = [c for c in Colors if c not in ['head', 'end']]
    return clist[cid % len(clist)]

def len_excluding_colors(seq):  # NOTE this won't work if you inserted a color code into the middle of another color code
    for color_code in Colors.values():
        seq = seq.replace(color_code, '')
    return len(seq)

def len_only_letters(seq):  # usually the same as len_excluding_colors(), except it doesn't count gap chars or spaces
    return len(list(filter((alphabet).__contains__, seq)))

def wrnstr():  # adding this very late, so could use it in a *lot* of places
    return color('yellow', 'warning')

def errstr():  # adding this very late, so could use it in a *lot* of places
    return color('red', 'error')

# ----------------------------------------------------------------------------------------
def color_chars(chars, col, seq):
    if sum([seq.count(c) for c in chars]) == 0:  # if <chars> aren't present, immediately return
        return seq
    return_str = [color(col, c) if c in chars else c for c in seq]
    return ''.join(return_str)

# ----------------------------------------------------------------------------------------
# returns a multiple sequence alignemnt from mafft (see run_blastn() below to align against a db of targets)
def align_many_seqs(seqfos, outfname=None, existing_aligned_seqfos=None, ignore_extra_ids=False, aa=False, no_gaps=False, extra_str='', debug=False):  # if <outfname> is specified, we just tell mafft to write to <outfname> and then return None
    def outfile_fcn():
        if outfname is None:
            return tempfile.NamedTemporaryFile(mode='w')
        else:
            return open(outfname, 'w')
    if existing_aligned_seqfos is not None and len(existing_aligned_seqfos) == 0:
        existing_aligned_seqfos = None

    with tempfile.NamedTemporaryFile(mode='w') as fin, outfile_fcn() as fout:
        for seqfo in seqfos:
            fin.write('>%s\n%s\n' % (seqfo['name'], seqfo['seq']))
        fin.flush()
        cmd = 'mafft'
        if no_gaps:
            cmd += ' --op 10'
        if existing_aligned_seqfos is None:  # default: align all the sequences in <seqfos>
            # subprocess.check_call('mafft --quiet %s >%s' % (fin.name, fout.name), shell=True)
            outstr, errstr = simplerun('%s --quiet %s >%s' % (cmd, fin.name, fout.name), shell=True, return_out_err=True, debug=False)
        else:  # if <existing_aligned_seqfos> is set, we instead add the sequences in <seqfos> to the alignment in <existing_aligned_seqfos>
            with tempfile.NamedTemporaryFile(mode='w') as existing_alignment_file:  # NOTE duplicates code in glutils.get_new_alignments()
                biggest_length = max(len(sfo['seq']) for sfo in existing_aligned_seqfos)
                for sfo in existing_aligned_seqfos:
                    dashstr = '-' * (biggest_length - len(sfo['seq']))
                    existing_alignment_file.write('>%s\n%s\n' % (sfo['name'], sfo['seq'].replace('.', '-') + dashstr))
                existing_alignment_file.flush()
                # subprocess.check_call('mafft --keeplength --add %s %s >%s' % (fin.name, existing_alignment_file.name, fout.name), shell=True)  #  --reorder
                outstr, errstr = simplerun('%s --keeplength --add %s %s >%s' % (cmd, fin.name, existing_alignment_file.name, fout.name), shell=True, return_out_err=True, debug=False)  #  --reorder

        if outfname is not None:
            if debug:
                print('  align_many_seqs(): wrote aligned seqs to %s' % outfname)
            return None

        msa_info = read_fastx(fout.name, ftype='fa')
        if existing_aligned_seqfos is not None:  # this may not be necessary, but may as well stay as consistent as possible
            for sfo in msa_info:
                sfo.update({'seq' : sfo['seq'].replace('-', '.')})

        input_ids = set([sfo['name'] for sfo in seqfos])
        output_ids = set([sfo['name'] for sfo in msa_info])
        missing_ids = input_ids - output_ids
        extra_ids = output_ids - input_ids
        if len(missing_ids) > 0 or (not ignore_extra_ids and len(extra_ids) > 0):
            print('  %d input ids not in output: %s' % (len(missing_ids), ' '.join(missing_ids)))
            print('  %d extra ids in output: %s' % (len(extra_ids), ' '.join(extra_ids)))
            print('  mafft out/err:')
            print(pad_lines(outstr))
            print(pad_lines(errstr))
            raise Exception('error reading mafft output from %s (see previous lines)' % fin.name)

    if debug:
        w = max(len(s['name']) for s in msa_info)
        for sfo in msa_info:
            print(color_mutants(msa_info[0]['seq'], sfo['seq'], ref_label=wfmt(msa_info[0]['name'], w)+' ', seq_label=wfmt(sfo['name'], w)+' ', amino_acid=aa, extra_str=extra_str))

    return msa_info

# ----------------------------------------------------------------------------------------
def align_seqs(ref_seq, seq):  # should eventually change name to align_two_seqs() or something
    with tempfile.NamedTemporaryFile(mode='w') as fin, tempfile.NamedTemporaryFile(mode='w') as fout:
        fin.write('>%s\n%s\n' % ('ref', ref_seq))
        fin.write('>%s\n%s\n' % ('new', seq))
        fin.flush()
        subprocess.check_call('mafft --quiet %s >%s' % (fin.name, fout.name), shell=True)
        msa_info = {sfo['name'] : sfo['seq'] for sfo in read_fastx(fout.name, ftype='fa')}
        if 'ref' not in msa_info or 'new' not in msa_info:
            subprocess.check_call(['cat', fin.name])
            raise Exception('incoherent mafft output from %s (cat\'d on previous line)' % fin.name)
    return msa_info['ref'], msa_info['new']

# ----------------------------------------------------------------------------------------
# darn it, maybe there was no reason to add this? I forgot that run_vsearch() seems to do basically the same thing? (although it would have needed some coding to use an arbitrary database)
def run_blastn(queryfos, targetfos, baseworkdir, diamond=False, short=False, aa=False, print_all_matches=False, debug=True):  # if short isn't set, it seems to ignore matches less than like 10 bases
    wkdir = '%s/blastn' % baseworkdir
    tgn = 'targets'
    tgfn, qrfn, ofn = [('%s/%s'%(wkdir, fstr)) for fstr in ['%s.fa'%tgn, 'queries.fa', 'results.out']]
    prep_dir(wkdir)
    write_fasta(tgfn, targetfos)
    write_fasta(qrfn, queryfos)
    qdict, tdict = [{s['name'] : s['seq'] for s in sfos} for sfos in [queryfos, targetfos]]  # should really check for duplicates
    if debug:
        print('    running blast on %s sequences with %d targets' % (len(queryfos), len(targetfos)))
    if diamond:
        assert False  # didn't end up implementing this since it turned out blastn does the same thing, and seemed fast enough for now (binary is still in bin/ though, and it should be easy) (note: could also use https://github.com/soedinglab/mmseqs2)
        dbcmd = './bin/diamond makedb --in reference.fasta -d reference'
    else:
        dbcmd = 'makeblastdb -in %s -out %s/%s -dbtype %s -parse_seqids' % (tgfn, wkdir, tgn, 'prot' if aa else 'nucl')
    _ = simplerun(dbcmd, return_out_err=True, debug=debug)
    shstr = (' -task blast%s-short'%('p' if aa else 'n')) if short else ''
    _ = simplerun('blast%s%s%s -db %s/%s -query %s -out %s -outfmt \"7 std qseq sseq btop\"' % ('p' if aa else 'n', shstr, '' if aa else ' -strand plus', wkdir, tgn, qrfn, ofn), shell=True, return_out_err=True, debug=debug)  # i'm just copying the format from somewhere else, there's probably a better one
    matchfos = collections.OrderedDict()
    if debug:
        max_qlen, max_tlen = [max(len(s['name']) for s in sfos) for sfos in [queryfos, targetfos]]
        max_tlen = max([max_tlen, len('(for first target)')])
        print('             %% id  mism. len    n gaps    %s            %s   target match seq' % (wfmt('target', max_tlen, jfmt='-'), wfmt('query', max_qlen, jfmt='-')))
    fieldnames = ['query', 'subject', '% identity', 'alignment length', 'mismatches', 'gap opens', 'q. start', 'q. end', 's. start', 's. end']
    with open(ofn) as ofile:
        reader = csv.DictReader(filter(lambda row: row[0]!='#', ofile), delimiter=str('\t'), fieldnames=fieldnames)
        for line in reader:
            qry, tgt, pct_id, alen, n_gaps = line['query'], line['subject'], float(line['% identity']), int(line['alignment length']), int(line['gap opens'])
            qbounds, tbounds = [[int(line['%s. start'%s]) - 1, int(line['%s. end'%s])] for s in ['q', 's']]
            if any(b[0] > b[1] for b in [qbounds, tbounds]):
                raise Exception('start bound larger than end bound: %s' % (list(b for b in [qbounds, tbounds] if b[0] > b[1])))
            def gseq(sq, bnd, obd, osq):
                return (obd[0] - bnd[0]) * gap_chars[0] + bnd[0] * gap_chars[0] + sq[bnd[0] : bnd[1]] + (len(sq) - bnd[1]) * gap_chars[0] + ((len(osq) - obd[1]) - (len(sq) - bnd[1])) * gap_chars[0]
            qseq = gseq(qdict[qry], qbounds, tbounds, tdict[tgt])
            tseq = gseq(tdict[tgt], tbounds, qbounds, qdict[qry])
            mfo = {'query' : qry, 'targets' : tgt, 'pct_ids' : pct_id, 'mismatches' : line['mismatches'], 'alens' : alen, 'n_gaps' : n_gaps, 'qbounds' : qbounds, 'tbounds' : tbounds, 'qseqs' : qseq, 'tseqs' : tseq}
            is_best = False
            if qry not in matchfos:  # it always puts the best match first, and i guess that could change but seems really unlikely
                matchfos[qry] = {k : v if k=='query' else [v] for k, v in mfo.items()}
                is_best = True
                if debug and print_all_matches:
                    print('')
            else:
                for key, val in mfo.items():
                    if key == 'query':
                        continue
                    matchfos[qry][key].append(val)
            if debug and (is_best or print_all_matches):
                if n_gaps == 0:
                    tgstr = color_mutants(qseq, tseq, amino_acid=aa)
                else:
                    tgstr = '  %s  ' % color('blue', 'GAPS')
                if is_best:
                    print('          %3s %3s  %3s   %3s    %s      %s      %s   %s' % (color('blue', '-->'), '', '', '', 3 * ' ', wfmt('(for first target)', max_tlen), wfmt(line['query'], max_qlen), qseq))
                print('          %3s %3.0f  %3d   %3d    %s      %s      %-s   %s' % ('', pct_id, int(line['mismatches']), alen, color(None if n_gaps==0 else 'red', '%d'%n_gaps, width=3), wfmt(tgt, max_tlen), max_qlen * ' ', tgstr))

    if len(matchfos) == 0:
        print('      no blast matches')
    mstats = {}  # ends up as a sorted list of pairs <target name, list of matched seqfos>
    def kfn(q): return q['targets'][0]
    for tgt, tgroup in itertools.groupby(sorted(list(matchfos.values()), key=kfn), key=kfn):
        mstats[tgt] = list(tgroup)
    mstats = sorted(list(mstats.items()), key=lambda x: len(x[1]), reverse=True)

    dbfns = ['%s/%s.%s%s'%(wkdir, tgn, 'p' if aa else 'n', s) for s in ('hr', 'in', 'og', 'sd', 'si', 'sq')]
    for fn in [tgfn, qrfn, ofn] + dbfns:
        os.remove(fn)
    os.rmdir(wkdir)

    return matchfos, mstats

# ----------------------------------------------------------------------------------------
def print_cons_seq_dbg(seqfos, cons_seq, aa=False, align=False, tie_resolver_seq=None, tie_resolver_label=None, extra_str='  ', dont_print_cons_seq=False):
    ckey = 'cons_dists_' + ('aa' if aa else 'nuc')
    hfkey = ckey.replace('cons_dists_', 'cons_fracs_')
    if len(seqfos) > 0 and ckey in seqfos[0]:
        seqfos = sorted(seqfos, key=lambda s: s[ckey])
        if hfkey in seqfos[0]:
            seqfos = sorted(seqfos, key=lambda s: s[hfkey])
    post_ref_str = '  hdist hfrac%s' % (('  %s'%color('blue', ' N')) if 'multiplicity' in seqfos[0] else '')
    for iseq, sfo in enumerate(seqfos):
        mstr = ''
        if 'multiplicity' in sfo:
            mstr = '%-3d' % sfo['multiplicity']
            if sfo['multiplicity'] > 1:
                mstr = color('blue', mstr)
        post_str = '   %2d  %5.2f    %s' % (sfo.get(ckey, -1), sfo.get(hfkey, -1), mstr)
        if 'cell-types' in sfo:  # if it's in one, it should be in all of 'em
            cl = str(max(len(str(s['cell-types'])) for s in seqfos))  # wasteful (but cleaner) to do this for eery sfo (also, str call is in case it's None)
            post_str += (' %'+cl+'s') % sfo['cell-types']
        color_mutants(cons_seq, sfo['seq'], align=align, amino_acid=aa, print_result=True, only_print_seq=iseq>0 or dont_print_cons_seq, ref_label=' consensus ', extra_str=extra_str, post_str='%s  %s'%(post_str, sfo['name']), post_ref_str=post_ref_str) #, print_n_snps=True, print_hfrac=True)
    if tie_resolver_seq is not None:
        color_mutants(cons_seq, tie_resolver_seq, align=align, amino_acid=aa, print_result=True, only_print_seq=True, seq_label=' '*len(' consensus '),
                      post_str='    tie resolver%s'%('' if tie_resolver_label is None else (' (%s)'%tie_resolver_label)), extra_str=extra_str, print_n_snps=True)

# ----------------------------------------------------------------------------------------
# return consensus of either aligned or unaligned sequences, in chunks of length <codon_len>, with tied positions *not* ambiguous but instead chosen as the alphabetically first character[s]
# if doing a nuc cons seq with codon_len=3, and <aa_ref_seq> is set, we look for the most common nuc codon *only* among those that code for the aa that appears at that position in <aa_ref_seq>
def cons_seq(aligned_seqfos=None, unaligned_seqfos=None, aa=False, codon_len=1, aa_ref_seq=None, extra_str='', debug=False):  # should maybe call it "chunk" rather than "codon", but if len is 3 it's a codon
    if aligned_seqfos is not None:
        assert unaligned_seqfos is None
        seqfos = aligned_seqfos
    elif unaligned_seqfos is not None:
        assert aligned_seqfos is None
        seqfos = align_many_seqs(unaligned_seqfos, aa=aa)
    else:
        assert False

    # for sfo in seqfos:  # this is pretty darn fast, but probably doesn't really need to be turned on
    #     check_alphabet(sfo['seq'], aa=aa, only_warn=True)

    seq_len = get_single_entry(list(set([len(s['seq']) for s in seqfos])))  # if this fails then the seqs that were passed in aren't all the same length
    if aa_ref_seq is not None:
        assert codon_len == 3 and not aa
        if 3 * len(aa_ref_seq) != len(pad_nuc_seq(seqfos[0]['seq'])):
            # print '\n'.join(s['seq'] for s in seqfos)  # TODO this probably means we're using input seqs for a family with lots of *different* indels, which needs to be fixed/avoided
            # raise Exception('aa ref seq length doesn\'t correspond to padded nuc seq:\n    %s\n    %s' % ('  '.join(aa_ref_seq), pad_nuc_seq(seqfos[0]['seq'])))
            print('%s aa ref seq length doesn\'t correspond to padded nuc seq, so turning off aa_ref_seq (this may mean that this family has lots of different indels, which isn\'t really properly handled):\n    %s\n    %s' % (color('red', 'error'), '  '.join(aa_ref_seq), pad_nuc_seq(seqfos[0]['seq'])))
            aa_ref_seq = None

    if debug:
        print('%staking consensus of %d seqs with len %d in chunks of len %d' % (extra_str, len(seqfos), seq_len, codon_len))
        dbgfo = []
        all_counts = {}  # keep track of usage of each codon/base/aa over full sequence for sorting of dbg info at end
    cseq = []
    for ipos in range(0, seq_len, codon_len):
        pos_counts = {}
        for sfo in seqfos:
            mtpy = sfo['multiplicity'] if 'multiplicity' in sfo else 1
            chnk = sfo['seq'][ipos : ipos + codon_len]  # if there's a partial codon at the end, it'll just stay as partial here, which seems fine (we could pad with pad_nuc_seq() if we wanted)
            if any(g in chnk for g in gap_chars):  # ok i probably shouldn't just skip them, but whatever
                continue
            if chnk not in pos_counts:
                pos_counts[chnk] = 0
            pos_counts[chnk] += mtpy
            if debug:
                if chnk not in all_counts:
                    all_counts[chnk] = 0
                all_counts[chnk] += mtpy
        if len(pos_counts) == 0:  # every sequence has a gap char here
            cseq.append(gap_chars[0])
            continue
        srt_chunks = sorted(list(pos_counts.items()), key=operator.itemgetter(1), reverse=True)
        if aa_ref_seq is not None:  # (try to) remove any that don't code for the residue in aa_ref_seq
            aa_match_chunks = [c for c, _ in srt_chunks if ltranslate(c) == aa_ref_seq[ipos // 3]]
            if debug > 1: removed_chunks = []  # arg
            if len(aa_match_chunks) > 0:  # if none match, we just have to keep all of them (this should be rare)
                if debug > 1:
                    removed_chunks = [c for c, _ in srt_chunks if c not in aa_match_chunks]
                srt_chunks = [(c, n) for c, n in srt_chunks if c in aa_match_chunks]
        n_max = srt_chunks[0][1]
        best_chunks = [c for c, n in srt_chunks if n == n_max]
        if debug:
            # print '      %3d  %s  %s' % (ipos, '  '.join('%s %d'%(color('red' if c in best_chunks and len(best_chunks) > 1 else None, c), n) for c, n in srt_chunks), '(removed %s)'%' '.join(removed_chunks) if aa_ref_seq is not None and debug>1 and len(removed_chunks)>0 else '')
            dbgfo.append({'best' : best_chunks, 'cnts' : pos_counts, 'rmd' : removed_chunks if aa_ref_seq is not None and debug>1 else []})
        cseq.append(sorted(best_chunks)[0])  # if there's more than one tied for most, take the first sorted alphabetically (arbitrary, but at least repeatable)

    if debug:
        prlen = max(4, codon_len) if debug > 1 else codon_len
        print('%s       ties: %s  (%s)' % (extra_str, ''.join(wfmt(len(dfo['best']), prlen, fmt='d') if len(dfo['best'])>1 else ' '*prlen for dfo in dbgfo), ',  '.join('%d: %s'%(codon_len*i, ' '.join(sorted(dfo['best']))) for i, dfo in enumerate(dbgfo) if len(dfo['best'])>1)))
        if codon_len == 3:
            print('%s             %s' % (extra_str, ''.join(wfmt(ltranslate(c), prlen) for c in cseq)))
        print('%s   cons seq: %s' % (extra_str, ''.join(color('blue_bkg' if len(dfo['best'])>1 else None, c, width=prlen if debug>1 else 1) for c, dfo in zip(cseq, dbgfo))))
        if debug > 1:
            def tcol(cstr, dfo): return 'blue' if cstr in dfo['best'] else ('red' if cstr in dfo['rmd'] and dfo['cnts'][cstr]>=dfo['cnts'][dfo['best'][0]] else None)
            print('        ipos %s' % ''.join(wfmt(codon_len*i, prlen, fmt='d') for i in range(len(dbgfo))))
            for cstr, _ in sorted(list(all_counts.items()), key=operator.itemgetter(1), reverse=True):
                print('        %s %s' % (wfmt(cstr, prlen), ''.join(color(tcol(cstr, dfo), str(dfo['cnts'][cstr]), width=prlen) if cstr in dfo['cnts'] else ' '*prlen for dfo in dbgfo)))
            print('      colors: %s %s %s' % (color('blue_bkg', 'ties'), color('blue', 'best'), color('red', 'removed (not aa_ref_seq)') if aa_ref_seq is not None else ''))

    cseq = ''.join(cseq)

    if debug:
        print_cons_seq_dbg(seqfos, cseq, aa=aa, align=aligned_seqfos is None, extra_str='  '+extra_str, dont_print_cons_seq=True)

    return cseq

# ----------------------------------------------------------------------------------------
def old_bio_cons_seq(threshold, aligned_seqfos=None, unaligned_seqfos=None, aa=False, tie_resolver_seq=None, tie_resolver_label=None, debug=False):
    """ return consensus sequence from either aligned or unaligned seqfos """
    # <threshold>: If the percentage*0.01 of the most common residue type is greater then the passed threshold, then we will add that residue type, otherwise an ambiguous character will be added. e.g. 0.1 means if fewer than 10% of sequences have the most common base, it gets an N.
    # <tie_resolver_seq>: in case of ties, use the corresponding base from this sequence (we used to use the naive sequence for this, but now I don't think that makes sense) NOTE if you *don't* set this argument, all tied bases will be Ns
    from io import StringIO
    from Bio.Align import AlignInfo
    import Bio.AlignIO

    if aligned_seqfos is not None:
        assert unaligned_seqfos is None
        seqfos = aligned_seqfos
    elif unaligned_seqfos is not None:
        assert aligned_seqfos is None
        seqfos = align_many_seqs(unaligned_seqfos, aa=aa)
    else:
        assert False

    def fstr(sfo): return '>%s\n%s' % (sfo['name'], sfo['seq'])
    if 'multiplicity' in seqfos[0]:
        fastalist = [fstr(sfo) for sfo in seqfos for _ in range(sfo['multiplicity'])]
    else:
        fastalist = [fstr(sfo) for sfo in seqfos]
    alignment = Bio.AlignIO.read(StringIO('\n'.join(fastalist) + '\n'), 'fasta')
    cons_seq = str(AlignInfo.SummaryInfo(alignment).gap_consensus(threshold, ambiguous=ambiguous_amino_acids[0] if aa else ambig_base))  # NOTE this is reliably *ten* goddamn times slower than my silly fcn ^ (*just* the gap_consensus() call, i.e. excluding everything else)

    if tie_resolver_seq is not None:  # huh, maybe it'd make more sense to just pass in the tie-breaker sequence to the consensus fcn?
        assert len(tie_resolver_seq) == len(cons_seq)
        cons_seq = list(cons_seq)
        for ipos in range(len(cons_seq)):
            if cons_seq[ipos] in all_ambiguous_bases:
                cons_seq[ipos] = tie_resolver_seq[ipos]
        cons_seq = ''.join(cons_seq)

    if debug:
        print_cons_seq_dbg(seqfos, cons_seq, aa=aa, align=aligned_seqfos is None, tie_resolver_seq=tie_resolver_seq, tie_resolver_label=tie_resolver_label)

    return cons_seq

# ----------------------------------------------------------------------------------------
def seqfos_from_line(line, aa=False, extra_keys=None, use_input_seqs=False, add_sfos_for_multiplicity=False, prepend_naive=False, naive_name='naive'):  # for search: get_seqfos()
    tstr = '_aa' if aa else ''
    skey = 'input_seqs' if use_input_seqs else 'seqs'
    seqfos = []
    if prepend_naive:
        seqfos.append({'name' : naive_name, 'seq' : line['naive_seq']})
    for uid, seq, mtp in zip(line['unique_ids'], line[skey+tstr], get_multiplicities(line)):
        if add_sfos_for_multiplicity:  # add the sfo <mtp> times
            for _ in range(mtp if mtp is not None else 1):
                seqfos.append({'name' : uid, 'seq' : seq})
        else:  # add it once, but with the 'multiplicity' key in the sfo
            seqfos.append({'name' : uid, 'seq' : seq, 'multiplicity' : mtp if mtp is not None else 1})
    if extra_keys is not None:
        for sfo in seqfos:
            iseq = None if sfo['name']==naive_name else line['unique_ids'].index(sfo['name'])
            for ekey in extra_keys:
                sfo[ekey] = line[ekey][iseq] if ekey in line and iseq is not None else None
    return seqfos

# ----------------------------------------------------------------------------------------
# NOTE does *not* add either 'consensus_seq' or 'consensus_seq_aa' to <line> (we want that to happen in the calling fcns)
def cons_seq_of_line(line, aa=False, use_input_seqs=False, codon_len=1, aa_ref_seq=None, debug=False):
# NOTE at least for now I'm not by default adding the codon_len=3 argument for the new cons seq fcn (on the theory that here we want the "regular" nuc cons seq, whereas the weird/fancy codon-based one is atm just when choosing abs)
    aligned_seqfos = seqfos_from_line(line, aa=aa, use_input_seqs=use_input_seqs)
    unaligned_seqfos = None
    # if any(indelutils.has_indels_line(line, i) for i in range(len(line['unique_ids']))):  NOTE *don't* use this, only align if you need to
    if len(set(len(s['seq']) for s in aligned_seqfos)) > 1:  # this will probably only happen if use_input_seqs is set and there's shm indels
        unaligned_seqfos = aligned_seqfos
        aligned_seqfos = None
    return cons_seq(aligned_seqfos=aligned_seqfos, unaligned_seqfos=unaligned_seqfos, aa=aa, codon_len=codon_len, aa_ref_seq=aa_ref_seq, debug=debug) # NOTE if you turn the naive tie resolver back on, you also probably need to uncomment in treeutils.add_cons_dists(), tie_resolver_seq=line['naive_seq'], tie_resolver_label='naive seq')
# ----------------------------------------------------------------------------------------
# and after that (below) we used to do this:
    # return old_bio_cons_seq(0.01, aligned_seqfos=aligned_seqfos, unaligned_seqfos=unaligned_seqfos, aa=aa, codon_len=codon_len) # NOTE if you turn the naive tie resolver back on, you also probably need to uncomment in treeutils.add_cons_dists(), tie_resolver_seq=line['naive_seq'], tie_resolver_label='naive seq')
# ----------------------------------------------------------------------------------------
    # Leaving the old version below for the moment just for reference.
    # It got the aa cons seq just by translating the nuc one, which is *not* what we want, since it can give you spurious ambiguous bases in the aa cons seq, e.g. if A and C tie at a position (so nuc cons seq has N there), but with either base it still codes for the same aa.
    # if aa:
    #     cseq = line['consensus_seq'] if 'consensus_seq' in line else cons_seq_of_line(line, aa=False)  # get the nucleotide cons seq, calculating it if it isn't already there NOTE do *not* use .get(), since in python all function arguments are evaluated *before* the call is excecuted, i.e. it'll call the consensus fcn even if the key is already there
    #     return ltranslate(cseq)
    # else:  # this is fairly slow
    #     aligned_seqfos = [{'name' : u, 'seq' : s, 'multiplicity' : m} for u, s, m in zip(line['unique_ids'], line['seqs'], get_multiplicities(line))]
    #     return cons_seq(aligned_seqfos=aligned_seqfos) # NOTE if you turn the naive tie resolver back on, you also probably need to uncomment in treeutils.add_cons_dists(), tie_resolver_seq=line['naive_seq'], tie_resolver_label='naive seq')
# ----------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------
def get_cons_seq_accuracy_vs_n_sampled_seqs(line, n_sample_min=7, n_sample_step=None, debug=False):  # yeah yeah the name is too long, but it's clear, isn't it?
    # NOTE i started to write this so that for small n_sampled it took several different subsamples, but i think that isn't a good idea since we want to make sure they're independent
    def get_step(n_sampled):
        if n_sample_step is not None:  # calling fcn can set it explicitly
            return n_sample_step
        if n_sampled < 10:
            return 1
        elif n_sampled < 20:
            return 3
        elif n_sampled < 30:
            return 4
        elif n_sampled < 40:
            return 5
        elif n_sampled < 85:
            return 10
        else:
            return 30
    def get_n_sample_list(n_total):
        nslist = [n_sample_min]
        while True:
            n_sampled = nslist[-1] + get_step(nslist[-1])
            if n_sampled > n_total:
                break
            nslist.append(n_sampled)
        return nslist

    n_total = len(line['unique_ids'])
    if n_total < n_sample_min:
        if debug:
            print('  small cluster %d' % n_total)
        return

    if debug:
        print('  getting cons seq accuracy for cluster with size %d' % n_total)
    ctypes = ['nuc', 'aa']
    info = {ct : {'n_sampled' : [], 'cseqs' : [], 'hdists' : None} for ct in ctypes}
    for n_sampled in get_n_sample_list(n_total):
        seqfos = [{'name' : line['unique_ids'][i], 'seq' : line['seqs'][i]} for i in numpy.random.choice(range(n_total), size=n_sampled, replace=False)]
        cseq = cons_seq(aligned_seqfos=seqfos)  # if we're sampling few enough that ties are frequent, then this doesn't really do what we want (e.g. there's a big oscillation for odd vs even n_sampled, I think because even ones have more ambiguous bases which don't count against them). But it doesn't matter, since ties are only at all frequent for families smaller than we care about (like, less than 10).
        for ctype in ctypes:
            info[ctype]['n_sampled'].append(n_sampled)
        info['nuc']['cseqs'].append(cseq)
        info['aa']['cseqs'].append(ltranslate(cseq))

    for ctype in ctypes:
        best_cseq = info[ctype]['cseqs'][-1]
        info[ctype]['hdists'] = [hamming_distance(cs, best_cseq, amino_acid=ctype=='aa') for cs in info[ctype]['cseqs']]  # it might make more sense for this to use hamming fraction, since small families get perhaps too much credit for ties that end up as ambiguous characters, but whatever
        if debug:
            print('  %s    N sampled   hdist   ' % color('green', ctype, width=6))
            for n_sampled, cseq, hdist in zip(info[ctype]['n_sampled'], info[ctype]['cseqs'], info[ctype]['hdists']):
                print('             %4d      %5.2f     %s' % (n_sampled, hdist, color_mutants(best_cseq, cseq, amino_acid=ctype=='aa')))

    return info

# ----------------------------------------------------------------------------------------
def color_mutants(ref_seq, seq, print_result=False, extra_str='', ref_label='', seq_label='', post_str='', post_ref_str='',
                  print_hfrac=False, print_isnps=False, return_isnps=False, print_n_snps=False, emphasis_positions=None, use_min_len=False,
                  only_print_seq=False, align=False, align_if_necessary=False, return_ref=False, amino_acid=False):  # NOTE if <return_ref> is set, the order isn't the same as the input sequence order
    """ default: return <seq> string with colored mutations with respect to <ref_seq> """

    # NOTE now I've got <return_ref>, I can probably remove a bunch of the label/whatever arguments and do all the damn formatting in the caller

    if use_min_len:
        min_len = min(len(ref_seq), len(seq))
        ref_seq = ref_seq[:min_len]
        seq = seq[:min_len]

    if align or (align_if_necessary and len(ref_seq) != len(seq)):  # it would be nice to avoid aligning when we don't need to... but i'm not sure how to identify cases where multiple indels result in the same length
        ref_seq, seq = align_seqs(ref_seq, seq)

    if len(ref_seq) != len(seq):
        raise Exception('unequal lengths in color_mutants()\n    %s\n    %s' % (ref_seq, seq))

    if amino_acid:
        tmp_ambigs = ambiguous_amino_acids
        tmp_gaps = gap_chars #[]
    else:
        tmp_ambigs = all_ambiguous_bases
        tmp_gaps = gap_chars

    return_str, isnps = [], []
    for inuke in range(len(seq)):  # would be nice to integrate this with hamming_distance() (especially the isnp stuff)
        rchar = ref_seq[inuke]
        char = seq[inuke]
        if char in tmp_ambigs or rchar in tmp_ambigs:
            char = color('blue', char)
        elif char in tmp_gaps or rchar in tmp_gaps:
            char = color('blue', char)
        elif char != rchar:
            char = color('red', char)
            isnps.append(inuke)
        if emphasis_positions is not None and inuke in emphasis_positions:
            char = color('reverse_video', char)
        return_str.append(char)

    isnp_str, n_snp_str = '', ''
    if len(isnps) > 0:
        if print_isnps:
            isnp_str = '   %d snp%s' % (len(isnps), plural(len(isnps)))
            if len(isnps) < 10:
                isnp_str +=  ' at: %s' % ' '.join([str(i) for i in isnps])
    if print_n_snps:
        n_snp_str = ' %3d' % len(isnps)
        if len(isnps) == 0:
            n_snp_str = color('blue', n_snp_str)
    hfrac_str = ''
    if print_hfrac:
        hfrac_str = '   hfrac %.3f' % hamming_fraction(ref_seq, seq, amino_acid=amino_acid)
    if print_result:
        lwidth = max(len_excluding_colors(ref_label), len_excluding_colors(seq_label))
        if not only_print_seq:
            ref_print_str = ''.join([color('blue' if c in tmp_ambigs + tmp_gaps else None, c) for c in ref_seq])
            print('%s%s%s%s%s' % (extra_str, ('%' + str(lwidth) + 's') % ref_label, ref_print_str, '  hdist' if print_n_snps else '', post_ref_str))
        print('%s%s%s' % (extra_str, ('%' + str(lwidth) + 's') % seq_label, ''.join(return_str) + n_snp_str + post_str + isnp_str + hfrac_str))

    return_list = [extra_str + seq_label + ''.join(return_str) + post_str + isnp_str + hfrac_str]
    if return_ref:
        return_list.append(''.join([ch if ch not in tmp_gaps else color('blue', ch) for ch in ref_seq]))
    if return_isnps:
        return_list.append(isnps)
    return return_list[0] if len(return_list) == 1 else return_list

# ----------------------------------------------------------------------------------------
def plural_str(pstr, count):
    if count == 1:
        return pstr
    else:
        return pstr + 's'

# ----------------------------------------------------------------------------------------
def plural(count, prefix=''):  # I should really combine these
    if count == 1:
        return ''
    else:
        return prefix + 's'

# ----------------------------------------------------------------------------------------
def summarize_gene_name(gene):
    region = get_region(gene)
    primary_version, sub_version, allele = split_gene(gene)
    return ' '.join([region, primary_version, sub_version, allele])

# ----------------------------------------------------------------------------------------
def color_gene(gene, width=None, leftpad=False, allow_constant=False):
    """ color gene name (and remove extra characters), eg IGHV3-h*01 --> hv3-h1 """
    default_widths = {'v' : 15, 'd' : 9, 'j' : 6}  # wide enough for most genes

    locus = get_locus(gene)
    locus = locus[2]  # hmm... maybe?
    region = get_region(gene, allow_constant=allow_constant)
    primary_version, sub_version, allele = split_gene(gene, allow_constant=allow_constant)

    n_chars = len(locus + region + primary_version)  # number of non-special characters
    return_str = color('purple', locus) + color('red', region) + color('purple', primary_version)
    if sub_version is not None:
        n_chars += 1 + len(sub_version)
        return_str += color('purple', '-' + sub_version)
    n_chars += len(allele)
    return_str += color('yellow', allele)
    if width is not None:
        if width == 'default':
            width = default_widths[region]
        if leftpad:
            return_str = (width - n_chars) * ' ' + return_str
        else:
            return_str += (width - n_chars) * ' '
    return return_str

# ----------------------------------------------------------------------------------------
def color_genes(genelist, allow_constant=False):  # now that I've added this fcn, I should really go and use this in all the places where the list comprehension is written out
    return ' '.join([color_gene(g, allow_constant=allow_constant) for g in genelist])

#----------------------------------------------------------------------------------------
def int_to_nucleotide(number):
    """ Convert between (0,1,2,3) and (A,C,G,T) """
    if number == 0:
        return 'A'
    elif number == 1:
        return 'C'
    elif number == 2:
        return 'G'
    elif number == 3:
        return 'T'
    else:
        print('ERROR nucleotide number not in [0,3]')
        sys.exit()

# ----------------------------------------------------------------------------------------
def count_n_separate_gaps(seq, exclusion_5p=None, exclusion_3p=None):  # NOTE compare to count_gap_chars() below (and gap_len() at top)
    if exclusion_5p is not None:
        seq = seq[exclusion_5p : ]
    if exclusion_3p is not None:
        seq = seq[ : len(seq) - exclusion_3p]

    n_gaps = 0
    within_a_gap = False
    for ch in seq:
        if ch not in gap_chars:
            within_a_gap = False
            continue
        elif within_a_gap:
            continue
        within_a_gap = True
        n_gaps += 1

    return n_gaps

# ----------------------------------------------------------------------------------------
def count_gap_chars(aligned_seq, aligned_pos=None, unaligned_pos=None):  # NOTE compare to count_n_separate_gaps() above (and gap_len() at top)
    """ return number of gap characters up to, but not including a position, either in unaligned or aligned sequence """
    if aligned_pos is not None:
        assert unaligned_pos is None
        aligned_seq = aligned_seq[ : aligned_pos]
        return sum([aligned_seq.count(gc) for gc in gap_chars])
    elif unaligned_pos is not None:
        assert aligned_pos is None
        ipos = 0  # position in unaligned sequence
        n_gaps_passed = 0  # number of gapped positions in the aligned sequence that we pass before getting to <unaligned_pos> (i.e. while ipos < unaligned_pos)
        while ipos < unaligned_pos or (ipos + n_gaps_passed < len(aligned_seq) and aligned_seq[ipos + n_gaps_passed] in gap_chars):  # second bit handles alignments with gaps immediately before <unaligned_pos>
            if ipos + n_gaps_passed == len(aligned_seq):  # i.e. <unaligned_pos> is just past the end of the sequence, i.e. slice notation for end of sequence (it would be nice to switch the while to a 'while True' and put all the logic into break statements, but I don't want to change the original stuff a.t.m.
                break
            if aligned_seq[ipos + n_gaps_passed] in gap_chars:
                n_gaps_passed += 1
            else:
                ipos += 1
        return n_gaps_passed
    else:
        assert False

# ----------------------------------------------------------------------------------------
def get_codon_pos_in_alignment(codon, aligned_seq, seq, pos, gene):  # NOTE see gap_len() and accompanying functions above
    """ given <pos> in <seq>, find the codon's position in <aligned_seq> """
    if not codon_unmutated(codon, seq, pos):  # this only gets called on the gene with the *known* position, so it shouldn't fail
        print('  %s mutated %s before alignment in %s' % (color('yellow', 'warning'), codon, gene))
    pos_in_alignment = pos + count_gap_chars(aligned_seq, unaligned_pos=pos)
    if not codon_unmutated(codon, aligned_seq, pos_in_alignment):
        print('  %s mutated %s after alignment in %s' % (color('yellow', 'warning'), codon, gene))
    return pos_in_alignment

# ----------------------------------------------------------------------------------------
def get_pos_in_alignment(aligned_seq, pos):  # kind of annoying to have this as well as get_codon_pos_in_alignment(), but I don't want to change that function's signature NOTE see gap_len() and accompanying functions above
    """ given <pos> in <seq>, find position in <aligned_seq> """
    return pos + count_gap_chars(aligned_seq, unaligned_pos=pos)

# ----------------------------------------------------------------------------------------
def both_codons_unmutated(locus, seq, positions, debug=False, extra_str=''):
    both_ok = True
    for region, codon in conserved_codons[locus].items():
        both_ok &= codon_unmutated(codon, seq, positions[region], debug=debug, extra_str=extra_str)
    return both_ok

# ----------------------------------------------------------------------------------------
def codon_unmutated(codon, seq, position, debug=False, extra_str=''):
    if len(seq) < position + 3:
        if debug:
            print('%ssequence length %d less than %s position %d + 3' % (extra_str, len(seq), codon, position))
        return False
    if seq[position : position + 3] not in codon_table[codon]:  # NOTE this allows it to be mutated to one of the other codons that codes for the same amino acid
        if debug:
            print('%s%s codon %s not among expected codons (%s)' % (extra_str, codon, seq[position : position + 3], ' '.join(codon_table[codon])))
        return False
    return True

#----------------------------------------------------------------------------------------
def in_frame_germline_v(seq, cyst_position):  # NOTE duplication with in_frame() (this is for when all we have is the germline v gene, whereas in_frame() is for when we have the whole rearrangement line)
    return cyst_position <= len(seq) - 3 and (cyst_position - count_gap_chars(seq, aligned_pos=cyst_position)) % 3 == 0

#----------------------------------------------------------------------------------------
def in_frame(seq, codon_positions, fv_insertion, v_5p_del, debug=False):  # NOTE I'm not passing the whole <line> in order to make it more explicit that <seq> and <codon_positions> need to correspond to each other, i.e. are either both for input seqs, or both for indel-reversed seqs
    # NOTE duplication with in_frame_germline_v()
    """ return true if the start and end of the cdr3 are both in frame with respect to the start of the V """
    germline_v_start = len(fv_insertion) - v_5p_del  # position in <seq> (the query sequence) to which the first base of the germline sequence aligns
    v_cpos = codon_positions['v'] - germline_v_start
    j_cpos = codon_positions['j'] - germline_v_start  # NOTE I'm actually not sure how necessary it is that the right side of the J is in frame. I mean, I think it's pretty framework-ey, but I'm not sure.
    if debug:
        print('    in frame:   v codon %d   j codon %d  -->  %s' % (v_cpos % 3 == 0, j_cpos % 3 == 0, v_cpos % 3 == 0 and j_cpos % 3 == 0))
    return v_cpos % 3 == 0 and j_cpos % 3 == 0

#----------------------------------------------------------------------------------------
# return list of the codons in <seq> (starting from first in frame position), along with any bits trimmed off the start and end
# NOTE duplicates some of pad_seq_for_translation()
def get_codon_list(seq, fv_insertion, jf_insertion, v_5p_del, debug=False):  # NOTE getting the indexing correct here is extremely non-trivial
    germline_v_start = len(fv_insertion) - v_5p_del  # position in <seq> (the query sequence) to which the first base of the germline sequence aligns
    istart = germline_v_start  # start with the first complete codon after <germline_v_start>
    while istart < len(fv_insertion):  # while staying in frame with the start of the v, skootch up to the first base in the query sequence that's actually aligned to the germline (i.e. up to 0 if no fv_insertion, and further if there is one)
        istart += 3
    germline_j_end = len(seq) - len(jf_insertion)  # position immediately after the end of the germline j (if there's a j_3p_del it's accounted for with len(seq))
    istop = germline_j_end - ((germline_j_end - istart) % 3)
    codons = [seq[i : i + 3] for i in range(istart, istop, 3)]
    if debug:
        print('   getting codons: istart %d  istop %d' % (istart, istop))
        print('     seq: %25s' % seq)
        print('  codons: %s' % ''.join([color('red' if cdn in codon_table['stop'] else None, cdn) for cdn in codons]))
    return codons, [seq[:istart], seq[istop:]]

#----------------------------------------------------------------------------------------
def is_there_a_stop_codon(seq, fv_insertion, jf_insertion, v_5p_del, return_n_stops=False, debug=False):
    """ true if there's a stop codon in frame with respect to the start of the V """
    codons, _ = get_codon_list(seq, fv_insertion, jf_insertion, v_5p_del, debug=debug)
    if debug:
        print('    %d stop codon%s'  % (len(set(codons) & set(codon_table['stop'])), plural(len(set(codons) & set(codon_table['stop'])))))
    if return_n_stops:
        return len([c for c in codons if c in codon_table['stop']])
    else:
        return len(set(codons) & set(codon_table['stop'])) > 0

# ----------------------------------------------------------------------------------------
def fix_stop(cdn):  # mutate codon <cdn> so it's no longer a stop codon
    while cdn in codon_table['stop']:
        imut = numpy.random.choice(range(len(cdn)))
        new_nuc = numpy.random.choice([n for n in nukes if n != cdn[imut]])
        cdn = cdn[:imut] + new_nuc + cdn[imut + 1:]
    return cdn

# ----------------------------------------------------------------------------------------
# mutate any stop codons until there's no stop codons left
def mutate_stop_codons(seq, fv_insertion, jf_insertion, v_5p_del, debug=False):
    codons, trim_bits = get_codon_list(seq, fv_insertion, jf_insertion, v_5p_del, debug=debug)
    for istp in [i for i, c in enumerate(codons) if c in codon_table['stop']]:
        new_cdn = fix_stop(codons[istp])
        if debug:
            print('  fixed %3d: %s --> %s' % (istp, codons[istp], new_cdn))
        codons[istp] = new_cdn
    seq = trim_bits[0] + ''.join(codons) + trim_bits[1]
    if debug:
        assert not is_there_a_stop_codon(seq, fv_insertion, jf_insertion, v_5p_del)
    return seq

# ----------------------------------------------------------------------------------------
def disambiguate_effective_insertions(bound, line, iseq, debug=False):
    # These are kinda weird names, but the distinction is important
    # If an insert state with "germline" N emits one of [ACGT], then the hmm will report this as an inserted N. Which is what we want -- we view this as a germline N which "mutated" to [ACGT].
    # This concept of insertion germline state is mostly relevant for simultaneous inference on several sequences, i.e. in ham we don't want to just say the inserted base was the base in the query sequence.
    # But here, we're trimming off the effective insertions and we have to treat the inserted germline N to [ACGT] differently than we would an insertion which was "germline" [ACGT] which emitted an N,
    # and also differently to a real germline [VDJ] state that emitted an N.
    naive_insertion = line[bound + '_insertion']  # reminder: ham gets this from the last character in the insertion state name, e.g. 'insert_left_A' or 'insert_right_N'
    insert_len = len(line[bound + '_insertion'])
    if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
        mature_insertion = line['seqs'][iseq][ : insert_len]
    elif bound == 'jf':
        mature_insertion = line['seqs'][iseq][len(line['seqs'][iseq]) - insert_len : ]
    else:
        assert False

    if naive_insertion == mature_insertion:  # all is simple and hunky-dory: no insertion 'mutations'
        final_insertion = ''  # leave this bit as an insertion in the final <line>
        insertion_to_remove = naive_insertion  # this bit we'll remove -- it's just Ns (note that this is only equal to the N padding if we correctly inferred the right edge of the J [for jf bound])
    else:
        if len(naive_insertion) != len(mature_insertion):
            raise Exception('naive and mature insertions not the same length\n   %s\n   %s\n' % (naive_insertion, mature_insertion))
        assert naive_insertion.count('N') == len(naive_insertion)  # should be e.g. naive: NNN   mature: ANN
        if bound == 'fv':  # ...but to accomodate multiple sequences, the insert states can emit non-germline states, so the mature bases might be different.
            i_first_non_N = find_first_non_ambiguous_base(mature_insertion)
            final_insertion = mature_insertion[i_first_non_N : ]
            insertion_to_remove = mature_insertion[ : i_first_non_N]
        elif bound == 'jf':
            i_first_N = find_last_non_ambiguous_base_plus_one(mature_insertion)
            final_insertion = mature_insertion[ : i_first_N]
            insertion_to_remove = mature_insertion[i_first_N : ]
        else:
            assert False

    if debug:
        print('     %s      final: %s' % (color('red', bound), color('purple', final_insertion)))
        print('         to_remove: %s' % color('blue', insertion_to_remove))

    return final_insertion, insertion_to_remove

# ----------------------------------------------------------------------------------------
# modify <line> so it has no 'fwk' insertions to left of v or right of j
def trim_fwk_insertions(glfo, line, modify_alternative_annotations=False, debug=False):  # NOTE this is *different* to reset_effective_erosions_and_effective_insertions() (i think kind of, but not entirely, the opposite?)
    from . import indelutils
    if line.get('is_fake_paired', False):
        raise Exception('doesn\'t work for fake paired annotations')
    # NOTE duplicates code in waterer.remove_framework_insertions(), and should really be combined with that fcn
    fv_len = len(line['fv_insertion'])
    jf_len = len(line['jf_insertion'])
    if debug:
        print('trimming fwk insertions: fv %d  jv %d' % (fv_len, jf_len))
        print_reco_event(line, label='before trimming:', extra_str='    ')

    if fv_len == 0 and jf_len == 0:
        return

    remove_all_implicit_info(line)

    for seqkey in ['seqs', 'input_seqs']:
        line[seqkey] = [seq[fv_len : len(seq) - jf_len] for i, seq in enumerate(line[seqkey])]
    for iseq in range(len(line['unique_ids'])):
        if indelutils.has_indels(line['indelfos'][iseq]):
            indelutils.trim_indel_info(line, iseq, line['fv_insertion'], line['jf_insertion'], 0, 0)
    line['fv_insertion'] = ''
    line['jf_insertion'] = ''

    if modify_alternative_annotations and 'alternative-annotations' in line:  # in principle it'd be nice to also generate these alternative naive seqs when we re-add implicit info, but we don't keep around near enough information to be able to do that
        for iseq, (seq, prob) in enumerate(line['alternative-annotations']['naive-seqs']):
            line['alternative-annotations']['naive-seqs'][iseq] = [seq[fv_len : len(seq) - jf_len], prob]

    add_implicit_info(glfo, line)

    if debug:
        print_reco_event(line, label='after trimming:', extra_str='    ')

# ----------------------------------------------------------------------------------------
def reset_effective_erosions_and_effective_insertions(glfo, padded_line, aligned_gl_seqs=None, debug=False):  # , padfo=None
    from . import indelutils
    """
    Ham does not allow (well, no longer allows) v_5p and j_3p deletions -- we instead pad sequences with Ns.
    This means that the info we get from ham always has these effective erosions set to zero, but for downstream
    things we sometimes want to know where the reads stopped (e.g. if we want to mimic them in simulation).
    Note that these effective erosion values will be present in the parameter dir, but are *not* incorporated into
    the hmm yaml files.
    """

    if debug:
        print('resetting effective erosions/insertions for %s' % ' '.join(padded_line['unique_ids']))

    line = {k : copy.deepcopy(padded_line[k]) for k in padded_line if k not in implicit_linekeys}

    assert line['v_5p_del'] == 0  # just to be safe
    assert line['j_3p_del'] == 0
    nseqs = len(line['unique_ids'])  # convenience

    # first disambiguate/remove effective (fv and jf) insertions
    if debug:
        print('   disambiguating effective insertions')
    trimmed_seqs = [line['seqs'][iseq] for iseq in range(nseqs)]
    trimmed_input_seqs = [line['input_seqs'][iseq] for iseq in range(nseqs)]
    final_insertions = [{} for _ in range(nseqs)]  # the effective insertions that will remain in the final info
    insertions_to_remove = [{} for _ in range(nseqs)]  # the effective insertions that we'll be removing, so which won't be in the final info
    for iseq in range(nseqs):
        fin = final_insertions[iseq]
        rem = insertions_to_remove[iseq]
        for bound in effective_boundaries:
            final_insertion, insertion_to_remove = disambiguate_effective_insertions(bound, line, iseq, debug=debug)
            fin[bound] = final_insertion
            rem[bound] = insertion_to_remove
        trimmed_seqs[iseq] = trimmed_seqs[iseq][len(rem['fv']) : len(trimmed_seqs[iseq]) - len(rem['jf'])]
        trimmed_input_seqs[iseq] = trimmed_input_seqs[iseq][len(rem['fv']) : len(trimmed_input_seqs[iseq]) - len(rem['jf'])]
        if debug:
            print('       %s  %s%s%s%s%s' % (' '.join(line['unique_ids']),
                                             color('blue', rem['fv']), color('purple', fin['fv']),
                                             trimmed_seqs[iseq][len(fin['fv']) : len(trimmed_seqs[iseq]) - len(fin['jf'])],
                                             color('purple', fin['jf']), color('blue', rem['jf'])))

    # arbitrarily use the zeroth sequence (in principle v_5p and j_3p should be per-sequence, not per-rearrangement... but that'd be a mess to implement, since the other deletions are per-rearrangement)
    TMPiseq = 0  # NOTE this is pretty hackey: we just use the values from the first sequence. But it's actually not that bad -- we can either have some extra pad Ns showing, or chop off some bases.
    trimmed_seq = trimmed_seqs[TMPiseq]
    final_fv_insertion = final_insertions[TMPiseq]['fv']
    final_jf_insertion = final_insertions[TMPiseq]['jf']
    fv_insertion_to_remove = insertions_to_remove[TMPiseq]['fv']
    jf_insertion_to_remove = insertions_to_remove[TMPiseq]['jf']

    def max_effective_erosion(erosion):  # don't "erode" more than there is left to erode
        region = erosion[0]
        gl_len = len(glfo['seqs'][region][line[region + '_gene']])
        if '5p' in erosion:
            other_del = line[region + '_3p_del']
        elif '3p' in erosion:
            other_del = line[region + '_5p_del']
        return gl_len - other_del - 1

    line['v_5p_del'] = min(max_effective_erosion('v_5p'), find_first_non_ambiguous_base(trimmed_seq))
    line['j_3p_del'] = min(max_effective_erosion('j_3p'), len(trimmed_seq) - find_last_non_ambiguous_base_plus_one(trimmed_seq))

    if debug:
        v_5p = line['v_5p_del']
        j_3p = line['j_3p_del']
        print('     %s:  %d' % (color('red', 'v_5p'), v_5p))
        print('     %s:  %d' % (color('red', 'j_3p'), j_3p))
        for iseq in range(nseqs):
            print('       %s  %s%s%s' % (' '.join(line['unique_ids']), color('red', v_5p * '.'), trimmed_seqs[iseq][v_5p : len(trimmed_seqs[iseq]) - j_3p], color('red', j_3p * '.')))

    for iseq in range(nseqs):
        line['seqs'][iseq] = trimmed_seqs[iseq][line['v_5p_del'] : len(trimmed_seqs[iseq]) - line['j_3p_del']]
        line['input_seqs'][iseq] = trimmed_input_seqs[iseq][line['v_5p_del'] : len(trimmed_input_seqs[iseq]) - line['j_3p_del']]
        if indelutils.has_indels_line(line, iseq):  # if 'indelfos' isn't in <line>, implicit info needs to be added, but maybe it will crash (arg)
            indelutils.trim_indel_info(line, iseq, fv_insertion_to_remove, jf_insertion_to_remove, line['v_5p_del'], line['j_3p_del'])

    line['fv_insertion'] = final_fv_insertion
    line['jf_insertion'] = final_jf_insertion

    # if padfo is None:
    #     line['padlefts'], line['padrights'] = [0 for _ in range(len(line['seqs']))], [0 for _ in range(len(line['seqs']))]
    # else:
    #     line['padlefts'], line['padrights'] = [padfo[uid]['padded']['padleft'] for uid in line['unique_ids']], [padfo[uid]['padded']['padright'] for uid in line['unique_ids']]

    # NOTE fixed the problem we were actually seeing, so this shouldn't fail any more, but I'll leave it in for a bit just in case UPDATE totally saved my ass from an unrelated problem (well, maybe not "saved" -- definitely don't remove the add_implicit_info() call though)
    try:
        add_implicit_info(glfo, line, aligned_gl_seqs=aligned_gl_seqs)
    except:
        exc_type, exc_value, exc_traceback = sys.exc_info()
        lines = traceback.format_exception(exc_type, exc_value, exc_traceback)
        print(pad_lines(''.join(lines)))
        print('  failed adding implicit info to \'%s\' (see above)' % ':'.join(line['unique_ids']))
        line['invalid'] = True

    return line

# ----------------------------------------------------------------------------------------
def add_qr_seqs(line):
    """ Add [vdj]_qr_seq, i.e. the sections of the query sequence which are assigned to each region. """

    starts = {}
    starts['v'] = len(line['fv_insertion'])
    starts['d'] = starts['v'] + len(line['v_gl_seq']) + len(line['vd_insertion'])
    starts['j'] = starts['d'] + len(line['d_gl_seq']) + len(line['dj_insertion'])

    def get_single_qr_seq(region, seq):
        return seq[starts[region] : starts[region] + len(line[region + '_gl_seq'])]

    for region in regions:
        line[region + '_qr_seqs'] = [get_single_qr_seq(region, seq) for seq in line['seqs']]

# ----------------------------------------------------------------------------------------
def is_functional_dbg_str(line, iseq, sep=', '):  # NOTE code duplication with is_functional(
    dbg_str_list = []
    if line['mutated_invariants'][iseq]:
        dbg_str_list.append('mutated invariant codon')
    if not line['in_frames'][iseq]:
        dbg_str_list.append('out of frame cdr3')
    if line['stops'][iseq]:
        dbg_str_list.append('stop codon')
    return sep.join(dbg_str_list)

# ----------------------------------------------------------------------------------------
# NOTE that this only looks at one chain, so it doesn't know anything about "non-pairing" ig chains (should add a link here to a reference)
def is_functional(line, iseq):  # NOTE code duplication with is_functional_dbg_str(
    if line['mutated_invariants'][iseq]:
        return False
    if not line['in_frames'][iseq]:
        return False
    if line['stops'][iseq]:
        return False
    return True

# ----------------------------------------------------------------------------------------
def add_functional_info(locus, line, input_codon_positions):
    nseqs = len(line['seqs'])  # would normally use 'unique_ids', but this gets called during simulation before the point at which we choose the uids
    line['mutated_invariants'] = [not both_codons_unmutated(locus, line['input_seqs'][iseq], input_codon_positions[iseq])
                                  for iseq in range(nseqs)]
    line['in_frames'] = [in_frame(line['input_seqs'][iseq], input_codon_positions[iseq], line['fv_insertion'], line['v_5p_del'])
                         for iseq in range(nseqs)]
    line['stops'] = [is_there_a_stop_codon(line['input_seqs'][iseq], line['fv_insertion'], line['jf_insertion'], line['v_5p_del'])
                     for iseq in range(nseqs)]

# ----------------------------------------------------------------------------------------
def remove_all_implicit_info(line):
    for col in implicit_linekeys:
        if col in line:
            del line[col]

# ----------------------------------------------------------------------------------------
def get_non_implicit_copy(line):  # return a deep copy of <line> with only non-implicit info
    return {col : copy.deepcopy(line[col]) for col in line if col not in implicit_linekeys}

# ----------------------------------------------------------------------------------------
def get_full_copy(line, glfo):  # NOTE this doesn't really make much sense (see next [commented] fcn), it's a placeholder til i get around to writing a cp fcn that avoids so many deepcopy() calls
    new_atn = get_non_implicit_copy(line)
    add_implicit_info(glfo, new_atn)
    return new_atn

# # ----------------------------------------------------------------------------------------
# deepcopy is just so freaking slow, it'd be nice to write something like this eventually
# def get_full_copy(line):
#     new_atn = {}
#     for key, val in line.items():
#         if isinstance(val, str) or isinstance(val, int) or isinstance(val, float):
#             new_atn[key] = val

# ----------------------------------------------------------------------------------------
def process_per_gene_support(line, debug=False):
    for region in regions:
        if debug:
            print(region)
        support = OrderedDict()
        logtotal = float('-inf')
        for gene, logprob in line[region + '_per_gene_support'].items():
            support[gene] = logprob
            logtotal = add_in_log_space(logtotal, logprob)

        for gene in support:
            if debug:
                print('   %5.2f     %5.2f   %s' % (support[gene], math.exp(support[gene] - logtotal), color_gene(gene)))
            support[gene] = math.exp(support[gene] - logtotal)

        if len(list(support.keys())) > 0 and list(support.keys())[0] != line[region + '_gene']:
            print('   %s best-supported gene %s not same as viterbi gene %s' % (color('yellow', 'warning'), color_gene(list(support.keys())[0]), color_gene(line[region + '_gene'])))

        line[region + '_per_gene_support'] = support

# ----------------------------------------------------------------------------------------
def re_sort_per_gene_support(line):
    for region in [r for r in regions if r + '_per_gene_support' in line]:
        if type(line[region + '_per_gene_support']) == type(collections.OrderedDict()):  # already ordered, don't need to do anything
            continue
        elif isinstance(line[region + '_per_gene_support'], dict):  # plain dict, i.e. we just read it from a json dump()'d file
            line[region + '_per_gene_support'] = collections.OrderedDict(sorted(list(line[region + '_per_gene_support'].items()), key=operator.itemgetter(1), reverse=True))

# ----------------------------------------------------------------------------------------
def get_null_linearham_info():
    return {'flexbounds' : None, 'relpos' : None}

# ----------------------------------------------------------------------------------------
def add_linearham_info(sw_info, annotation_list, glfo, min_cluster_size=None, debug=False):
    print('  adding linearham info')
    n_already_there, n_size_skipped, n_added = 0, 0, 0
    for line in annotation_list:
        trim_fwk_insertions(glfo, line)
        if min_cluster_size is not None and len(line['unique_ids']) < min_cluster_size:
            if debug:
                print('       %s adding null linearham info to line: %s, because cluster is less than the passed <min_cluster_size> value of %d' % (color('yellow', 'warning'), ':'.join(line['unique_ids']), min_cluster_size))
            line['linearham-info'] = get_null_linearham_info()
            n_size_skipped += 1
            continue
        if 'linearham-info' in line:
            if debug:
                print('       %s overwriting linearham info that was already in line: %s' % (color('yellow', 'warning'), ':'.join(line['unique_ids'])))
            n_already_there += 1
        else:
            n_added += 1
        line['linearham-info'] = get_linearham_bounds(sw_info, line, debug=debug)  # note that we don't skip ones that fail, since we don't want to just silently ignore some of the input sequences -- skipping should happen elsewhere where it can be more explicit
    if n_already_there > 0:
        print('    %s overwriting %d / %d that already had linearham info' % (color('yellow', 'warning'), n_already_there, len(annotation_list)))
    if n_size_skipped > 0:
        print('    skipped %d / %d with size less than %d' % (n_size_skipped, len(annotation_list), min_cluster_size))
    if n_added > 0:
        print('    added linearham info for %d clusters' % n_added)
    else:
        print('    %s didn\'t add any linearham info' % wrnstr())

# ----------------------------------------------------------------------------------------
def get_linearham_bounds(sw_info, line, vj_flexbounds_shift=10, debug=False):
    """ compute the flexbounds/relpos values and return in a dict """  # NOTE deep copies per_gene_support, and then modifies this copy
    def get_swfo(uid):
        def getmatches(matchfo):  # get list of gene matches sorted by decreasing score
            genes, gfos = zip(*sorted(list(matchfo.items()), key=lambda x: x[1]['score'], reverse=True))
            return genes
        swfo = {'flexbounds' : {}, 'relpos' : {}}
        for region in getregions(get_locus(line['v_gene'])):
            matchfo = sw_info[uid]['all_matches'][0][region]
            if isinstance(matchfo, list):
                raise Exception('\'all_matches\' key in sw info doesn\'t have qr/gl bound information, so can\'t add linearham info (it\'s probably an old sw cache file from before we started storing this info -- re-cache-parameters to rewrite it)')
            sortmatches = getmatches(matchfo)
            bounds_l, bounds_r = zip(*[matchfo[g]['qrbounds'] for g in sortmatches])  # left- (and right-) bounds for each gene
            swfo['flexbounds'][region + '_l'] = dict(zip(sortmatches, bounds_l))
            swfo['flexbounds'][region + '_r'] = dict(zip(sortmatches, bounds_r))
            for gene, gfo in matchfo.items():
                swfo['relpos'][gene] = gfo['qrbounds'][0] - gfo['glbounds'][0]  # position in the query sequence of the start of each uneroded germline match
        return swfo

    fbounds = {}  # initialize the flexbounds/relpos dicts
    rpos = {}

    for region in getregions(get_locus(line['v_gene'])):
        left_region, right_region = region + '_l', region + '_r'
        fbounds[left_region] = {}
        fbounds[right_region] = {}

    # rank the query sequences according to their consensus sequence distances
    cons_seq = ''.join([Counter(site_bases).most_common()[0][0] for site_bases in zip(*line['seqs'])])
    dists_to_cons = {line['unique_ids'][i] : hamming_distance(cons_seq, line['seqs'][i]) for i in range(len(line['unique_ids']))}

    # loop over the ranked query sequences and update the flexbounds/relpos dicts
    while len(dists_to_cons) > 0:
        query_name = min(dists_to_cons, key=dists_to_cons.get)
        swfo = get_swfo(query_name)

        for region in getregions(get_locus(line['v_gene'])):
            left_region, right_region = region + '_l', region + '_r'
            fbounds[left_region] = dict(list(swfo['flexbounds'][left_region].items()) + list(fbounds[left_region].items()))
            fbounds[right_region] = dict(list(swfo['flexbounds'][right_region].items()) + list(fbounds[right_region].items()))
        rpos = dict(list(swfo['relpos'].items()) + list(rpos.items()))

        del dists_to_cons[query_name]

    # 1) restrict flexbounds/relpos to gene matches with non-zero support
    # 2) if the left/right flexbounds overlap in a particular region, remove the worst gene matches until the overlap disappears
    # 3) if neighboring flexbounds overlap across different regions, apportion the flexbounds until the overlap disappears
    # 4) if possible, widen the gap between neighboring flexbounds
    # 5) align the V-entry/J-exit flexbounds to the sequence bounds

    def are_fbounds_empty(fbounds, region, gene_removed, reason_removed):
        left_region, right_region = region + '_l', region + '_r'
        if len(list(fbounds[left_region].values())) == 0 or len(list(fbounds[right_region].values())) == 0:
            print('{}: removed all genes from flexbounds for region {}: {}. The last gene removed was {}. It was removed because {}. Returning null linearham info'.format(color('yellow', 'warning'), region, fbounds, gene_removed, reason_removed))
            return True
        return False

    def span(bound_list):
        return [min(bound_list), max(bound_list)]

    for region in getregions(get_locus(line['v_gene'])):
        # EH: remember when reading this that left_region and right_region are not one of (v,d,j) as variables named like *region* often are in partis. Here they have _l or _r on the end so they are a particlar end of a region
        left_region, right_region = region + '_l', region + '_r'
        per_gene_support = copy.deepcopy(line[region + '_per_gene_support'])
        # remove the gene matches with zero support
        for k in list(fbounds[left_region].keys()):
            support_check1 = (k not in per_gene_support)
            support_check2 = (math.fabs(per_gene_support[k] - 0.0) < eps) if not support_check1 else False
            if support_check1 or support_check2:
                del fbounds[left_region][k]
                del fbounds[right_region][k]
                del rpos[k]
                if debug:
                    print('removing %s from fbounds (and per_gene_support if it was there to begin with) for region %s because it was not in per_gene_support or had too low support.' % (k, region))
                if support_check2:
                    del per_gene_support[k]
            if are_fbounds_empty(fbounds, region, k, '{} was not in per_gene_support or had too low support'.format(k)):
                return get_null_linearham_info()

        # compute the initial left/right flexbounds
        left_flexbounds = span(list(fbounds[left_region].values()))
        right_flexbounds = span(list(fbounds[right_region].values()))
        germ_len = right_flexbounds[0] - left_flexbounds[1]
        # make sure there is no overlap between the left/right flexbounds
        while germ_len < 1:
            k = min(per_gene_support, key=per_gene_support.get)
            del fbounds[left_region][k]
            del fbounds[right_region][k]
            del rpos[k]
            del per_gene_support[k]
            if debug:
                print('removing %s from fbounds and perg_gene_support to resolve a supposed overlap berween left and right flexbounds for region %s' % (k, region))
            # check removing all items from fbounds 
            if are_fbounds_empty(fbounds, region, k, 'right flexbounds was less than left for {}'.format(region)):
                return get_null_linearham_info()
            left_flexbounds = span(list(fbounds[left_region].values()))
            right_flexbounds = span(list(fbounds[right_region].values()))
            germ_len = right_flexbounds[0] - left_flexbounds[1]
        fbounds[left_region] = left_flexbounds
        fbounds[right_region] = right_flexbounds

    # make sure there is no overlap between neighboring flexbounds
    # maybe widen the gap between neighboring flexbounds
    for rpair in region_pairs(get_locus(line['v_gene'])):
        # EH: remember when reading this that left_region and right_region are not one of (v,d,j) as variables named like *region* often are in partis. Here they have _l or _r on the end so they are a particlar end of a region
        left_region, right_region = rpair['left'] + '_r', rpair['right'] + '_l'
        leftleft_region, rightright_region = rpair['left'] + '_l', rpair['right'] + '_r'

        left_germ_len = fbounds[left_region][0] - fbounds[leftleft_region][1]
        junction_len = fbounds[right_region][1] - fbounds[left_region][0]
        right_germ_len = fbounds[rightright_region][0] - fbounds[right_region][1]

        if junction_len < 1:
            if debug:
                print('''
                          Overlap resolution code running in partis utils.get_linearham_bounds.
                          Post Duncan's fix to ig-sw (see *), we really should not have a true overlap of neighboring matches.
                          So ideally this code would not ever get triggered. However, this code does get triggered if there
                          are adjacent matches which share regional bounds, which is possible because partis uses python slice
                          conventions for regional bounds (so this is not a "real" overlap).
                          E.g. if fbounds[left_region] = [x,x] and fbounds[right_region] = [x,x] as well,
                          this code gets executed despite this not being a true overlap.
                          However in such a case this code does nothing so we are not changing it for fear of messing up the logic here.
                          *: https://github.com/psathyrella/partis/commit/471e5eac6d2b0fbdbb2b6024c81af14cdc3d9399
                      ''')
            fbounds[left_region][0] = fbounds[right_region][0]
            fbounds[right_region][1] = fbounds[left_region][1]

            left_germ_len = fbounds[left_region][0] - fbounds[leftleft_region][1]
            junction_len = fbounds[right_region][1] - fbounds[left_region][0]
            right_germ_len = fbounds[rightright_region][0] - fbounds[right_region][1]

            # are the neighboring flexbounds even fixable?
            if left_germ_len < 1 or right_germ_len < 1:
                print('    failed adding linearham info for line %s due to overlapping neighboring fbounds between %s and %s' % (':'.join(line['unique_ids']), left_region, right_region))
                return get_null_linearham_info()

        # EH: This section corresponds to step #4 in the comment earlier in this fcn. Note that both the lower and upper bounds are shifted away from their neighboring gene in all cases here. This might seem odd, since performing such a shift on just one bound might help account for more uncertainty at the junction of each pair of genes, but shifting both bounds in the same direction wouldn't appear to have that effect. However, because linearham only cares about the bound furthest from the neighboring gene when considering a junction between two genes, the end result of shifting both bounds in this logic here is the same as if we had just shifted one bound in each gene (linearham's bound of interest for the gene) and the desired effect of the shift is achieved, i.e. that we just allow for some extra flexibility/uncertainty in these junction regions.
        if rpair['left'] == 'v' and left_germ_len > vj_flexbounds_shift:
            if debug:
                print('shifting lower and uppper fbounds for %s by %d' % (left_region, vj_flexbounds_shift))
            fbounds[left_region][0] -= vj_flexbounds_shift
            fbounds[left_region][1] -= vj_flexbounds_shift

        # the D gene match region is constrained to have a length of 1
        if rpair['left'] == 'd':
            if debug:
                print('shifting lower and uppper fbounds for %s by %d' % (left_region, left_germ_len - 1))
            fbounds[left_region][0] -= (left_germ_len - 1)
            fbounds[left_region][1] -= (left_germ_len - 1)
        if rpair['right'] == 'd':
            if debug:
                print('shifting lower and uppper fbounds for %s by %d' % (right_region, right_germ_len // 2))
            fbounds[right_region][0] += (right_germ_len // 2)
            fbounds[right_region][1] += (right_germ_len // 2)

        if rpair['right'] == 'j' and right_germ_len > vj_flexbounds_shift:
            if debug:
                print('shifting lower and uppper fbounds for %s by %d' % (right_region, vj_flexbounds_shift))
            fbounds[right_region][0] += vj_flexbounds_shift
            fbounds[right_region][1] += vj_flexbounds_shift

    # align the V-entry/J-exit flexbounds to the possible sequence positions
    fbounds['v_l'][0] = 0
    fbounds['j_r'][1] = len(line['seqs'][0])  # remember indel-reversed sequences are all the same length

    # are the V-entry/J-exit flexbounds valid?
    for region in ['v_l', 'j_r']:
        bounds_len = fbounds[region][1] - fbounds[region][0]
        if bounds_len < 0:
            print('    failed adding linearham info for line: %s . fbounds is negative length for region %s' % (':'.join(line['unique_ids']), region))
            return get_null_linearham_info()

    return {'flexbounds' : fbounds, 'relpos' : rpos}

# ----------------------------------------------------------------------------------------
def check_per_seq_lengths(line, bkey='unique_ids'):  # maybe there's something like this elsewhere? adding it very late. In any case there's lots of places i should add this check, i've screwed this up several times before and it's *really* hard to catch
    blen = len(line[bkey])
    for tkey in [k for k in line if k in linekeys['per_seq']]:
        if len(line[tkey]) != blen:
            raise Exception('per seq line keys \'%s\' %d and \'%s\' %d not the same length' % (bkey, blen, tkey, len(line[tkey])))

# ----------------------------------------------------------------------------------------
def add_implicit_info(glfo, line, aligned_gl_seqs=None, check_line_keys=False, reset_indel_genes=False):  # should turn on <check_line_keys> for a bit if you change anything
    from . import indelutils
    from . import glutils
    """ Add to <line> a bunch of things that are initially only implicit. """
    if line['v_gene'] == '':
        raise Exception('can\'t add implicit info to line with failed annotation:\n%s' % (''.join(['  %+20s  %s\n' % (k, v) for k, v in line.items()])))

    if check_line_keys:
        initial_keys = set(line)
        # first make sure there aren't any unauthorized keys
        if len(initial_keys - all_linekeys) > 0:
            raise Exception('unexpected keys: \'%s\'' % '\' \''.join(initial_keys - all_linekeys))
        # then keep track of the keys we got to start with
        pre_existing_implicit_info = {ek : copy.deepcopy(line[ek]) for ek in implicit_linekeys if ek in line}

    for region in regions:  # backwards compatibility with old simulation files should be removed when you're no longer running on them
        if line[region + '_gene'] not in glfo['seqs'][region]:
            alternate_name = glutils.convert_to_duplicate_name(glfo, line[region + '_gene'])
            # print ' using alternate name %s instead of %s' % (alternate_name, line[region + '_gene'])
            line[region + '_gene'] = alternate_name

    # add the regional germline seqs and their lengths
    line['lengths'] = {}  # length of each match (including erosion)
    for region in regions:
        uneroded_gl_seq = glfo['seqs'][region][line[region + '_gene']]
        del_5p = line[region + '_5p_del']
        del_3p = line[region + '_3p_del']
        length = len(uneroded_gl_seq) - del_5p - del_3p  # eroded length
        if length < 0:
            raise Exception('invalid %s lengths passed to add_implicit_info()\n    gl seq: %d  5p: %d  3p: %d' % (region, len(uneroded_gl_seq), del_5p, del_3p))
        line[region + '_gl_seq'] = uneroded_gl_seq[del_5p : del_5p + length]
        line['lengths'][region] = length

    # add codon-related stuff
    line['codon_positions'] = {}
    for region, codon in conserved_codons[glfo['locus']].items():
        eroded_gl_pos = glfo[codon + '-positions'][line[region + '_gene']] - line[region + '_5p_del']
        if region == 'v':
            line['codon_positions'][region] = eroded_gl_pos + len(line['f' + region + '_insertion'])
        elif region == 'j':
            line['codon_positions'][region] = eroded_gl_pos + len(line['fv_insertion']) + line['lengths']['v'] + len(line['vd_insertion']) + line['lengths']['d'] + len(line['dj_insertion'])
        else:
            assert False
    line['cdr3_length'] = line['codon_positions']['j'] - line['codon_positions']['v'] + 3  # i.e. first base of cysteine to last base of tryptophan inclusive

    # add naive seq stuff
    line['naive_seq'] = len(line['fv_insertion']) * ambig_base + line['v_gl_seq'] + line['vd_insertion'] + line['d_gl_seq'] + line['dj_insertion'] + line['j_gl_seq'] + len(line['jf_insertion']) * ambig_base
    for iseq in range(len(line['seqs'])):
        if len(line['naive_seq']) != len(line['seqs'][iseq]):
            color_mutants(line['naive_seq'], line['seqs'][iseq], ref_label='naive ', seq_label=line['unique_ids'][iseq] + ' ', align_if_necessary=True, print_result=True, extra_str='        ')
            raise Exception('naive and mature sequences different lengths %d %d for %d-seq annotation %s (see also previous stdout lines):\n    %s\n    %s' % (len(line['naive_seq']), len(line['seqs'][iseq]), len(line['unique_ids']), ':'.join(line['unique_ids']), line['naive_seq'], line['seqs'][iseq]))

    start, end = {}, {}  # add naive seq bounds for each region (could stand to make this more concise)
    start['v'] = len(line['fv_insertion'])  # NOTE this duplicates code in add_qr_seqs()
    end['v'] = start['v'] + len(line['v_gl_seq'])  # base just after the end of v
    start['d'] = end['v'] + len(line['vd_insertion'])
    end['d'] = start['d'] + len(line['d_gl_seq'])
    start['j'] = end['d'] + len(line['dj_insertion'])
    end['j'] = start['j'] + len(line['j_gl_seq'])
    line['regional_bounds'] = {r : (start[r], end[r]) for r in regions}

    try:
        indelutils.deal_with_indel_stuff(line, reset_indel_genes=reset_indel_genes)
    except indelutils.IndelfoReconstructionError:  # I don't like this here, but see note in the one place it can be raised
        line['invalid'] = True
        return

    input_codon_positions = [indelutils.get_codon_positions_with_indels_reinstated(line, iseq, line['codon_positions']) for iseq in range(len(line['seqs']))]
    if 'indel_reversed_seqs' not in line:  # everywhere internally, we refer to 'indel_reversed_seqs' as simply 'seqs'. For interaction with outside entities, however (i.e. writing files) we use the more explicit 'indel_reversed_seqs'
        line['indel_reversed_seqs'] = line['seqs']

    # add regional query seqs
    add_qr_seqs(line)

    add_functional_info(glfo['locus'], line, input_codon_positions)

    hfracfo = [hamming_fraction(line['naive_seq'], mature_seq, also_return_distance=True) for mature_seq in line['seqs']]
    line['mut_freqs'] = [hfrac for hfrac, _ in hfracfo]
    line['n_mutations'] = [n_mutations for _, n_mutations in hfracfo]

    # set validity (alignment addition [below] can also set invalid)  # it would be nice to clean up this checking stuff
    line['invalid'] = False
    seq_length = len(line['seqs'][0])  # they shouldn't be able to be different lengths
    for chkreg in regions:
        if start[chkreg] < 0 or end[chkreg] < 0 or end[chkreg] < start[chkreg] or end[chkreg] > seq_length:
            line['invalid'] = True
    if end['j'] + len(line['jf_insertion']) != seq_length:
        line['invalid'] = True
    if line['cdr3_length'] < 6:  # i.e. if cyst and tryp overlap  NOTE six is also hardcoded in waterer
        line['invalid'] = True

    # add alignment info (this is only used if presto output has been specified on the command line, which requires specification of your own alignment file)
    if aligned_gl_seqs is None:  # empty/dummy info
        for region in regions:
            line['aligned_' + region + '_seqs'] = ['' for _ in range(len(line['seqs']))]
    else:
        add_alignments(glfo, aligned_gl_seqs, line)

    re_sort_per_gene_support(line)  # in case it was read from json dump()'d file

    # for pskey in set(linekeys['per_seq']) & set(line):  # make sure all the per-seq keys have the right length (it happened once, admittedly cause i was editing by hand, but the consequences were extremely bad)
    #     if len(line[pskey]) != len(line['unique_ids']):
    #         raise Exception('line with %d uids has %d values for key \'%s\'' % (len(line['unique_ids']), len(line[pskey]), pskey))

    if check_line_keys:
        new_keys = set(line) - initial_keys
        if len(new_keys - implicit_linekeys) > 0:
            raise Exception('added new keys that aren\'t in implicit_linekeys: %s' % ' '.join(new_keys - implicit_linekeys))
        for ikey in implicit_linekeys:  # make sure every key/value we added is either a) new or b) the same as it was before
            if ikey in initial_keys:
                if pre_existing_implicit_info[ikey] != line[ikey]:
                    print('%s pre-existing info for \'%s\' in %s\n    %s\n    doesn\'t match new info\n    %s' % (color('yellow', 'warning'), ikey, line['unique_ids'], pre_existing_implicit_info[ikey], line[ikey]))
            else:
                assert ikey in new_keys  # only really checks the logic of the previous few lines

# ----------------------------------------------------------------------------------------
def restrict_to_iseqs(line, iseqs_to_keep, glfo, sw_info=None, remove_tree=False):  # could have called it subset_seqs_in_line, or at least i always seem to search for that when i'm trying to find this (or subset_iseqs or subset_to_iseqs)
    # NOTE if you want to return a new one rather than modifying <line>, call get_non_implicit_copy() on <line> as you pass it in
    """ remove from <line> any sequences corresponding to indices not in <iseqs_to_keep>. modifies line. """
    if len(iseqs_to_keep) < 1:  # NOTE we could also return if we're keeping all of them, but then we have to mess with maybe adding implicit info
        raise Exception('must be called with at least one sequence to keep (got %s)' % iseqs_to_keep)
    remove_all_implicit_info(line)
    add_per_seq_keys(line)
    for tkey in set(linekeys['per_seq']) & set(line):
        line[tkey] = [line[tkey][iseq] for iseq in iseqs_to_keep]
    add_implicit_info(glfo, line)
    if remove_tree and 'tree' in line:  # would need to collapse a bunch of nodes in the tree, too hard to do right now
        line['tree'] = None
    if line.get('linearham-info') is not None:
        if sw_info is not None:
            add_linearham_info(sw_info, [line], glfo)
        else:
            print('% restrict_to_iseqs(line, iseqs_to_keep, glfo, sw_info=None) needs sw_info to re-add existing \'linearham-info\' key to an annotation' % color('yellow', 'warning'))

# ----------------------------------------------------------------------------------------
def print_true_events(glfo, reco_info, inf_line, print_naive_seqs=False, full_true_partition=None, extra_str='    '):
    """ print the true events which contain the seqs in <inf_line> """
    if print_naive_seqs:
        true_naive_seqs = []
    true_partition_of_line_uids = get_partition_from_reco_info(reco_info, ids=uids_and_dups(inf_line))  # *not* in general the same clusters as in the complete true partition, since <inf_line['unique_ids']> may not contain all uids from all clusters from which it contains representatives
    if len(true_partition_of_line_uids) > 1:
        print('%s true clusters %d for this inferred cluster (size %d including duplicates)' % (color('red', 'multiple'), len(true_partition_of_line_uids), len(uids_and_dups(inf_line))))
    if full_true_partition is None:
        full_true_partition = get_partition_from_reco_info(reco_info)
    for itc, tpl_ids in enumerate(true_partition_of_line_uids):  # <tpl_ids>: ids in inf_line, split according to true partition
        full_true_cluster = get_single_entry([c for c in full_true_partition if len(set(c) & set(tpl_ids)) > 0])
        multiple_str = '' if len(true_partition_of_line_uids)==1 else '  true cluster index %d %s for inf cluster of size %d' % (itc, color('red', '(of %d)' % len(true_partition_of_line_uids)), len(uids_and_dups(inf_line)))
        missing_uids = set(full_true_cluster) - set(tpl_ids)
        missing_str = ''
        if len(missing_uids) > 0:
            n_dups = len(uids_and_dups(inf_line)) - len(tpl_ids)
            missing_str = '   %s %d/%d sequences from actual true cluster%s' % (color('red', 'missing'), len(missing_uids), len(full_true_cluster), ' (but includes %d duplicates not shown below)'%n_dups if n_dups>0 else '')

        multiline = synthesize_multi_seq_line_from_reco_info(tpl_ids, reco_info)
        ixstr = extra_str
        if inf_line['fv_insertion'] != '' and multiline['fv_insertion'] == '':
            ixstr = ' '*len(inf_line['fv_insertion']) + ixstr  # aligns true + inferred vertically
        print_reco_event(multiline, extra_str=ixstr, label=color('green', 'true:') + multiple_str + missing_str)
        if print_naive_seqs:
            true_naive_seqs.append(multiline['naive_seq'])

    if print_naive_seqs:
        print('      naive sequences:')
        for tseq in true_naive_seqs:
            color_mutants(tseq, inf_line['naive_seq'], print_result=True, print_hfrac=True, ref_label='true ', extra_str='          ')

# ----------------------------------------------------------------------------------------
def get_uid_extra_strs(line, extra_print_keys, uid_extra_strs, uid_extra_str_label):
    def vstr(val):
        if isinstance(val, float): return ('%.3f'%val)
        elif isinstance(val, list): return ':'.join(str(w) for w in val)
        else: return str(val)
    if uid_extra_strs is None:
        uid_extra_strs = ['' for _ in line['unique_ids']]
    if uid_extra_str_label is None:
        uid_extra_str_label = ''
    for ekey in extra_print_keys:
        vlist = [antnval(line, ekey, i, use_default=True, default_val=color('blue', '-'), add_xtr_col=True) for i in range(len(line['unique_ids']))]
        # tw = str(max(len(ekey), max(len(vstr(v)) for v in vlist)))  # maybe include len of ekey in width?
        tw = max([len(vstr(v)) for v in vlist] + [len(ekey)])
        uid_extra_str_label += '  ' + wfmt(ekey, tw, jfmt='-')
        assert len(vlist) == len(uid_extra_strs)
        uid_extra_strs = [('%s  %s'%(e, wfmt(vstr(v), tw, jfmt='-'))) for v, e in zip(vlist, uid_extra_strs)]
    return uid_extra_strs, uid_extra_str_label

# ----------------------------------------------------------------------------------------
def print_reco_event(line, extra_str='', label='', post_label='', uid_extra_strs=None, uid_extra_str_label=None, extra_print_keys=None, queries_to_emphasize=None):
    if line['invalid']:
        print('%s%s %s tried to print invalid event' % (extra_str, label, wrnstr()))
        return
    from . import prutils
    if extra_print_keys is not None:
        uid_extra_strs, uid_extra_str_label = get_uid_extra_strs(line, extra_print_keys, uid_extra_strs, uid_extra_str_label)
    if uid_extra_strs is not None and len(uid_extra_strs) != len(line['unique_ids']):
        raise Exception('uid_extra_strs %d different length to unique_ids %d' % (len(uid_extra_strs), len(line['unique_ids'])))
    duplicate_counts = [(u, line['unique_ids'].count(u)) for u in line['unique_ids']]
    duplicated_uids = {u : c for u, c in duplicate_counts if c > 1}
    if len(line['unique_ids']) > 1:
        label += '%s%d sequences with %.1f mean mutations (%.1f%%)' % ('' if label == '' else '    ', len(line['unique_ids']), numpy.mean(line['n_mutations']), 100*numpy.mean(line['mut_freqs']))
    for iseq in range(len(line['unique_ids'])):
        prutils.print_seq_in_reco_event(line, iseq, extra_str=extra_str, label=(label + post_label if iseq==0 else ''), one_line=(iseq>0), queries_to_emphasize=queries_to_emphasize, duplicated_uids=duplicated_uids, uid_extra_str=uid_extra_strs[iseq] if uid_extra_strs is not None else '', uid_extra_str_label=uid_extra_str_label)

#----------------------------------------------------------------------------------------
def sanitize_name(name):
    """ Replace characters in gene names that make crappy filenames. """
    saniname = name.replace('*', '_star_')
    saniname = saniname.replace('/', '_slash_')
    return saniname

#----------------------------------------------------------------------------------------
def unsanitize_name(name):
    """ Re-replace characters in gene names that make crappy filenames. """
    unsaniname = name.replace('_star_', '*')
    unsaniname = unsaniname.replace('_slash_', '/')
    return unsaniname

# ----------------------------------------------------------------------------------------
def get_locus(inputstr):
    """ return locus given gene or gl fname """
    # if len(inputstr) < 4:
    #     raise Exception('input str \'%s\' too short for get_locus() (need 4 chars)' % inputstr)
    locus = inputstr[:3].lower()  # only need the .lower() if it's a gene name
    if locus not in loci:
        raise Exception('couldn\'t get locus from input string \'%s\'' % inputstr)
    return locus

# ----------------------------------------------------------------------------------------
def get_region(inputstr, allow_constant=False):
    """ return v, d, or j of gene or gl fname """
    region = inputstr[3].lower()  # only need the .lower() if it's a gene name
    if not allow_constant:
        allowed_regions = regions
    else:
        allowed_regions = regions + constant_regions
    if region not in allowed_regions:
        raise Exception('unexpected region %s from %s (expected one of %s)' % (region, inputstr, allowed_regions))
    return region

# ----------------------------------------------------------------------------------------
def are_alleles(gene1, gene2):
    return primary_version(gene1) == primary_version(gene2) and sub_version(gene1) == sub_version(gene2)

# ----------------------------------------------------------------------------------------
def construct_valid_gene_name(gene, locus=None, region=None, default_allele_str='x', debug=False):  # kind of duplicates too much of split_gene(), but I don't want to rewrite split_gene() to be robust to all the ways a gene name can be broken
    try:  # if it's ok, don't do anything
        split_gene(gene)
        return gene
    except:
        pass

    if debug:
        initial_name = gene

    if len(gene) < 4 or gene[:3].lower() not in loci or gene[3].lower() not in regions:
        if locus is not None or region is not None:
            gene = locus.upper() + region.upper() + gene  # just tack (e.g.) 'IGHV' on the fron of whatever crap was originall there
        else:
            raise Exception('gene name %s doesn\'t have locus/region info, and it wasn\'t passed to us' % gene)

    if debug:
        middle_name = gene

    if gene.count('*') == 0:
        gene = gene + '*' + default_allele_str
    elif gene.count('*') > 1:
        gene = gene.replace('*', '.s.') + '*' + default_allele_str

    if debug:
        print('  %-25s  -->  %-25s  --> %-25s' % (initial_name, middle_name if middle_name != initial_name else '-', gene if gene != middle_name else '-'))

    return gene

# ----------------------------------------------------------------------------------------
def split_gene(gene, allow_constant=False):
    """ returns (primary version, sub version, allele) """
    # make sure {IG,TR}{[HKL],[abgd]}[VDJ] is at the start, and there's a *
    if '_star_' in gene or '_slash_' in gene:
        raise Exception('gene name \'%s\' isn\'t entirely unsanitized' % gene)
    if gene[:4] != get_locus(gene).upper() + get_region(gene, allow_constant=allow_constant).upper():
        raise Exception('unexpected string in gene name %s' % gene)
    if gene.count('*') != 1:
        raise Exception('expected exactly 1 \'*\' in %s but found %d' % (gene, gene.count('*')))

    if '-' in gene and gene.find('-') < gene.find('*'):  # Js (and a few Vs) don't have sub versions
        primary_version = gene[4 : gene.find('-')]  # the bit between the IG[HKL][VDJ] and the first dash (sometimes there's a second dash as well)
        sub_version = gene[gene.find('-') + 1 : gene.find('*')]  # the bit between the first dash and the star
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != get_locus(gene).upper() + get_region(gene).upper() + primary_version + '-' + sub_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s %s' % (gene, primary_version, sub_version, allele))
    else:
        primary_version = gene[4 : gene.find('*')]  # the bit between the IG[HKL][VDJ] and the star
        sub_version = None
        allele = gene[gene.find('*') + 1 : ]  # the bit after the star
        if gene != get_locus(gene).upper() + get_region(gene, allow_constant=allow_constant).upper() + primary_version + '*' + allele:
            raise Exception('couldn\'t build gene name %s from %s %s' % (gene, primary_version, allele))

    return primary_version, sub_version, allele

# ----------------------------------------------------------------------------------------
def shorten_gene_name(name, use_one_based_indexing=False, n_max_mutstrs=3):
    from . import glutils
    if name[:2] != 'IG':
        raise Exception('bad node name %s' % name)

    pv, sv, al = split_gene(name)
    if glutils.is_novel(name):
        _, template_name, mutstrs = glutils.split_inferred_allele_name(name)
        if use_one_based_indexing:
            mutstrs = [('%s%d%s' % (mstr[0], int(mstr[1:-1]) + 1, mstr[-1])) for mstr in mutstrs]
        if mutstrs is None:
            al = '%s (+...)' % (allele(template_name))
        elif len(mutstrs) < n_max_mutstrs:
            al = ('%s+%s' % (allele(template_name), '.'.join(mutstrs)))
        else:
            al = '%s (+%d snp%s)' % (allele(template_name), len(mutstrs), plural(len(mutstrs)))
    if sv is not None:
        return '%s-%s*%s' % (pv, sv, al)
    else:
        return '%s*%s' % (pv, al)

# ----------------------------------------------------------------------------------------
def rejoin_gene(locus, region, primary_version, sub_version, allele):
    """ reverse the action of split_gene() """
    return_str = locus.upper() + region.upper() + primary_version
    if sub_version is not None:  # e.g. J genes typically don't have sub-versions
        return_str += '-' + sub_version
    return return_str + '*' + allele

# ----------------------------------------------------------------------------------------
def primary_version(gene):
    return split_gene(gene)[0]

# ----------------------------------------------------------------------------------------
def gene_family(gene):  # same as primary_version(), except ignore stuff after the slash, e.g. 1/OR15 --> 1
    return primary_version(gene).split('/')[0].replace('D', '')

# ----------------------------------------------------------------------------------------
def sub_version(gene):
    return split_gene(gene)[1]

# ----------------------------------------------------------------------------------------
def allele(gene):
    return split_gene(gene)[2]

# ----------------------------------------------------------------------------------------
def is_constant_gene(gene):
    region = get_region(gene, allow_constant=True)
    if region not in constant_regions:
        return False
    if region != 'd':  # d is in both constant and regular regions
        return True
    pv, sv, allele = split_gene(gene)
    if pv == '' and sv is None:  # constant region d is like IGHD*01 (no pv)
        return True
    return False

# ----------------------------------------------------------------------------------------
def are_same_primary_version(gene1, gene2):
    """
    Return true if the bit up to the dash is the same.
    """
    if get_region(gene1) != get_region(gene2):
        return False
    if primary_version(gene1) != primary_version(gene2):
        return False
    return True

# ----------------------------------------------------------------------------------------
def separate_into_allelic_groups(glfo, allele_prevalence_freqs=None, debug=False):  # prevalence freqs are just for printing
    allelic_groups = {r : {} for r in regions}
    for region in regions:
        for gene in glfo['seqs'][region]:
            primary_version, sub_version, allele = split_gene(gene)
            if primary_version not in allelic_groups[region]:
                allelic_groups[region][primary_version] = {}
            if sub_version not in allelic_groups[region][primary_version]:
                allelic_groups[region][primary_version][sub_version] = set()
            allelic_groups[region][primary_version][sub_version].add(gene)
    if debug:
        for r in regions:
            print('%s%77s' % (color('reverse_video', color('green', r)), 'percent prevalence' if allele_prevalence_freqs is not None else ''))
            for p in sorted(allelic_groups[r]):
                print('    %15s' % p)
                for s in sorted(allelic_groups[r][p]):
                    print('        %15s      %s' % (s, ' '.join([color_gene(g, width=14) for g in allelic_groups[r][p][s]])), end=' ')
                    if len(allelic_groups[r][p][s]) < 2:  # won't work if anybody has more than two alleles
                        print('%14s' % '', end=' ')
                    if allele_prevalence_freqs is not None:
                        print('  %s' % ' '.join([('%4.1f' % (100 *allele_prevalence_freqs[r][g])) for g in allelic_groups[r][p][s]]), end=' ')
                    print('')
    return allelic_groups  # NOTE doesn't return the same thing as separate_into_snp_groups()

# ----------------------------------------------------------------------------------------
def separate_into_snp_groups(glfo, region, n_max_snps, genelist=None, debug=False):  # NOTE <n_max_snps> corresponds to v, whereas d and j are rescaled according to their lengths
    """ where each class contains all alleles with the same length (up to cyst if v), and within some snp threshold (n_max_v_snps for v)"""
    from . import glutils
    def getseq(gene):
        seq = glfo['seqs'][region][gene]
        if region == 'v':  # only go up through the end of the cysteine
            cpos = cdn_pos(glfo, region, gene)
            seq = seq[:cpos + 3]
        return seq
    def in_this_class(classfo, seq):
        for gfo in classfo:
            if len(gfo['seq']) != len(seq):
                continue
            hdist = hamming_distance(gfo['seq'], seq)
            if hdist < n_max_snps:  # if this gene is close to any gene in the class, add it to this class
                snp_groups[snp_groups.index(classfo)].append({'gene' : gene, 'seq' : seq, 'hdist' : hdist})
                return True
        return False  # if we fall through, nobody in this class was close to <seq>

    if genelist is None:
        genelist = list(glfo['seqs'][region].keys())
    snp_groups = []
    for gene in genelist:
        seq = getseq(gene)
        add_new_class = True  # to begin with, assume we'll add a new class for this gene
        for classfo in snp_groups:  # then check if, instead, this gene belongs in any of the existing classes
            if in_this_class(classfo, seq):
                add_new_class = False
                break
        if add_new_class:
            snp_groups.append([{'gene' : gene, 'seq' : seq, 'hdist' : 0}, ])

    if debug:
        print('separated %s genes into %d groups separated by %d snps:' % (region, len(snp_groups), n_max_snps))
        glutils.print_glfo(glfo, gene_groups={region : [(str(isub), {gfo['gene'] : gfo['seq'] for gfo in snp_groups[isub]}) for isub in range(len(snp_groups))]})

    return snp_groups  # NOTE this is a list of lists of dicts, whereas separate_into_allelic_groups() returns a dict of region-keyed dicts

# ----------------------------------------------------------------------------------------
def read_single_gene_count(indir, gene, expect_zero_counts=False, debug=False):
    region = get_region(gene)
    count = 0
    with open(indir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
        reader = csv.DictReader(infile)
        for line in reader:
            if line[region + '_gene'] == gene:
                count = int(line['count'])
                break

    if count == 0 and not expect_zero_counts:
        print('          %s %s not found in %s_gene-probs.csv, returning zero' % (color('red', 'warning'), gene, region))

    if debug:
        print('    read %d observations of %s from %s' % (count, color_gene(gene), indir))

    return count

# ----------------------------------------------------------------------------------------
def read_overall_gene_probs(indir, only_gene=None, normalize=True, expect_zero_counts=False, debug=False):
    """
    Return the observed counts/probabilities of choosing each gene version.
    If <normalize> then return probabilities
    If <only_gene> is specified, just return the prob/count for that gene  NOTE but don't forget read_single_gene_count() above ^, which I think probably does the same thing
    """
    counts, probs = {r : {} for r in regions}, {r : {} for r in regions}
    for region in regions:
        total = 0
        with open(indir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
            reader = csv.DictReader(infile)
            for line in reader:
                line_count = int(line['count'])
                gene = line[region + '_gene']
                total += line_count
                if gene not in counts[region]:
                    counts[region][gene] = 0
                counts[region][gene] += line_count
        if total < 1:
            raise Exception('less than one count in %s' % indir + '/' + region + '_gene-probs.csv')
        for gene in counts[region]:
            probs[region][gene] = float(counts[region][gene]) / total

    if debug:
        for region in regions:
            print('  %s' % color('green', region))
            for gene, count in sorted(list(counts[region].items()), key=operator.itemgetter(1), reverse=True):
                print('    %5d  %5.4f   %s' % (count, probs[region][gene], color_gene(gene, width='default')))

    if only_gene is not None and only_gene not in counts[get_region(only_gene)]:
        if not expect_zero_counts:
            print('      WARNING %s not found in overall gene probs, returning zero' % only_gene)
        if normalize:
            return 0.0
        else:
            return 0

    if only_gene is None:
        if normalize:
            return probs
        else:
            return counts
    else:
        if normalize:
            return probs[get_region(only_gene)][only_gene]
        else:
            return counts[get_region(only_gene)][only_gene]

# ----------------------------------------------------------------------------------------
def get_genes_with_enough_counts(parameter_dir, min_prevalence_fractions, debug=False):
    if debug:
        print('  applying min gene prevalence fractions: %s' % '  '.join(('%s %.4f' % (r, min_prevalence_fractions[r])) for r in regions))
    gene_freqs = read_overall_gene_probs(parameter_dir, normalize=True, debug=debug)
    genes_with_enough_counts = set([g for r in regions for g, f in gene_freqs[r].items() if f > min_prevalence_fractions[r]])  # this is kind of weird because <gene_freqs> of course normalizes within each region, but then we mash all the regions together in the list, but it's all ok
    if debug:
        print('   removed genes:')
        for region in regions:
            genes_without_enough_counts = set(gene_freqs[region]) - genes_with_enough_counts
            print('      %s   %s' % (color('green', region), color_genes(sorted(genes_without_enough_counts))))
    return genes_with_enough_counts

# ----------------------------------------------------------------------------------------
def find_replacement_genes(param_dir, min_counts, gene_name=None, debug=False, all_from_region=''):  # NOTE if <gene_name> isn't in <param_dir>, it won't be among the returned genes
    if gene_name is not None:  # if you specify <gene_name> you shouldn't specify <all_from_region>
        assert all_from_region == ''
        region = get_region(gene_name)
    else:  # and vice versa
        assert all_from_region in regions
        assert min_counts == -1
        region = all_from_region
    lists = OrderedDict()  # we want to try alleles first, then primary versions, then everything and it's mother
    lists['allele'] = []  # list of genes that are alleles of <gene_name>
    lists['primary_version'] = []  # same primary version as <gene_name>
    lists['all'] = []  # give up and return everything
    with open(param_dir + '/' + region + '_gene-probs.csv', 'r') as infile:  # NOTE note this ignores correlations... which I think is actually ok, but it wouldn't hurt to think through it again at some point
        reader = csv.DictReader(infile)
        for line in reader:
            gene = line[region + '_gene']
            count = int(line['count'])
            vals = {'gene':gene, 'count':count}
            if all_from_region == '':
                if are_alleles(gene, gene_name):
                    lists['allele'].append(vals)
                if are_same_primary_version(gene, gene_name):
                    lists['primary_version'].append(vals)
            lists['all'].append(vals)

    if all_from_region != '':
        return [vals['gene'] for vals in lists['all']]
    for list_type in lists:
        total_counts = sum([vals['count'] for vals in lists[list_type]])
        if total_counts >= min_counts:
            return_list = [vals['gene'] for vals in lists[list_type]]
            if debug:
                print('      returning all %s for %s (%d gene%s, %d total counts)' % (list_type + 's', color_gene(gene_name), len(return_list), plural(len(return_list)), total_counts))
            return return_list
        else:
            if debug:
                print('      not enough counts in %s' % (list_type + 's'))

    raise Exception('couldn\'t find enough counts among genes for %s in %s (found %d, needed %d -- to decrease this minimum set --min-observations-per-gene, although note that you\'re probably getting this exception because you have too few events to have very informative distributions)' % (gene_name, param_dir, total_counts, min_counts))

    # print '    \nWARNING return default gene %s \'cause I couldn\'t find anything remotely resembling %s' % (color_gene(hackey_default_gene_versions[region]), color_gene(gene_name))
    # return hackey_default_gene_versions[region]

# ----------------------------------------------------------------------------------------
def hamming_distance(seq1, seq2, extra_bases=None, return_len_excluding_ambig=False, return_mutated_positions=False, align=False, align_if_necessary=False, amino_acid=False):
    if extra_bases is not None:
        raise Exception('not sure what this was supposed to do (or did in the past), but it doesn\'t do anything now! (a.t.m. it seems to only be set in bin/plot-germlines.py, which I think doesn\'t do anything useful any more)')
    if align or (align_if_necessary and len(seq1) != len(seq2)):  # way the hell slower if you have to align, of course
        if align_if_necessary:
            print('  %s unequal length sequences %d %d (so aligning):\n  %s\n  %s' % (color('yellow', 'warning'), len(seq1), len(seq2), seq1, seq2))
        seq1, seq2 = align_seqs(seq1, seq2)
    if len(seq1) != len(seq2):
        s2str, s1str = color_mutants(seq1, seq2, align_if_necessary=True, return_ref=True)
        raise Exception('unequal length sequences (%d vs %d) in hamming distance:\n  %s\n  %s' % (len(seq1), len(seq2), s1str, s2str))
    if len(seq1) == 0:
        if return_len_excluding_ambig:
            return 0, 0
        else:
            return 0

    if amino_acid:
        skip_chars = set(ambiguous_amino_acids + gap_chars)
    else:
        skip_chars = set(all_ambiguous_bases + gap_chars)

    distance, len_excluding_ambig = 0, 0
    mutated_positions = []
    for ich in range(len(seq1)):  # already made sure they're the same length
        if seq1[ich] in skip_chars or seq2[ich] in skip_chars:
            continue
        len_excluding_ambig += 1
        if seq1[ich] != seq2[ich]:
            distance += 1
            if return_mutated_positions:
                mutated_positions.append(ich)

    if return_len_excluding_ambig and return_mutated_positions:
        return distance, len_excluding_ambig, mutated_positions
    elif return_len_excluding_ambig:
        return distance, len_excluding_ambig
    elif return_mutated_positions:
        return distance, mutated_positions
    else:
        return distance

# ----------------------------------------------------------------------------------------
def hamming_fraction(seq1, seq2, extra_bases=None, also_return_distance=False, amino_acid=False, align_if_necessary=False):  # NOTE use hamming_distance() to get the positions (yeah, I should eventually add it here as well)
    distance, len_excluding_ambig = hamming_distance(seq1, seq2, extra_bases=extra_bases, return_len_excluding_ambig=True, amino_acid=amino_acid, align_if_necessary=align_if_necessary)

    fraction = 0.
    if len_excluding_ambig > 0:
        fraction = distance / float(len_excluding_ambig)

    if also_return_distance:
        return fraction, distance
    else:
        return fraction

# ----------------------------------------------------------------------------------------
def get_mut_positions(line):  # for each sequence, return a list of position indices at which it is mutated
    hdistfo = [hamming_distance(line['naive_seq'], mature_seq, return_mutated_positions=True) for mature_seq in line['seqs']]
    return [mps for _, mps in hdistfo]

# ----------------------------------------------------------------------------------------
def get_mut_codes(naive_seq, obs_seq, amino_acid=False, debug=False):  # return list of mutation "events" e.g. [{'pos' : 3, 'initial' : 'A', 'final' : 'C', 'str' : 'A3C'}, ]
    hdist, mutated_positions = hamming_distance(naive_seq, obs_seq, return_mutated_positions=True, amino_acid=amino_acid)
    mcodes = []
    for ipos in mutated_positions:
        mcd = {'pos' : ipos, 'initial' : naive_seq[ipos], 'final' : obs_seq[ipos]}
        mcd['str'] = '%s%d%s' % (mcd['initial'], mcd['pos'], mcd['final'])
        mcodes.append(mcd)

    if debug:
        color_mutants(naive_seq, obs_seq, print_result=True, amino_acid=amino_acid)
        print(''.join(color('red', 'x') if i in mutated_positions else ' ' for i in range(len(naive_seq))))
        print('  '.join('%s-->%s'%(c['initial'], c['final']) for c in mcodes))
    return mcodes

# ----------------------------------------------------------------------------------------
def mean_pairwise_hfrac(seqlist, amino_acid=False, n_max_seqs=None):
    return mean_pairwise_dist(seqlist, hamming_fraction, amino_acid=amino_acid, n_max_seqs=n_max_seqs)

# ----------------------------------------------------------------------------------------
def mean_pairwise_hdist(seqlist, amino_acid=False, n_max_seqs=None):
    return mean_pairwise_dist(seqlist, hamming_distance, amino_acid=amino_acid, n_max_seqs=n_max_seqs)

# ----------------------------------------------------------------------------------------
def mean_pairwise_dist(seqlist, dfcn, amino_acid=False, n_max_seqs=None):
    if len(seqlist) < 2:
        return 0.
    if n_max_seqs is not None and len(seqlist) > n_max_seqs:
        old_len = len(seqlist)
        seqlist = numpy.random.choice(seqlist, n_max_seqs, replace=False)
    return numpy.mean([dfcn(s1, s2, amino_acid=amino_acid) for s1, s2 in itertools.combinations(seqlist, 2)])

# ----------------------------------------------------------------------------------------
def lev_dist(s1, s2, aa=False):  # NOTE does *not* handle ambiguous characters correctly (also NOTE <aa> has no effect
    if 'Levenshtein' not in sys.modules:  # installing this reliably is being super annoying, and we only using it in the ab choice stuff
        import Levenshtein
    return sys.modules['Levenshtein'].distance(s1, s2)

# ----------------------------------------------------------------------------------------
# return list of families in <antn_pairs> sorted by their nearness to <refpair> by 'lev' (levenshtein), 'ham' (hamming) distance, or blast mismathces between naive sequences (automatically skips <refpair> if it's in <antn_pairs>)
# For single chain just set that chain's entry in <refpair> and <antn_pairs> to None
def non_clonal_clusters(refpair, antn_pairs, dtype='lev', aa=False, workdir=None, max_print_dist=16, max_n_print=5, extra_str='', labelstr='', debug=True):
    # ----------------------------------------------------------------------------------------
    def nseq(atn_pair):
        tkey = 'naive_seq'
        if aa:
            tkey += '_aa'
            for tl in [l for l in atn_pair if l is not None]:
                add_naive_seq_aa(tl)
        return ''.join(l[tkey] for l in atn_pair if l is not None)
    # ----------------------------------------------------------------------------------------
    def gids(atn):
        if atn is None:
            return []
        return atn['unique_ids']
    # ----------------------------------------------------------------------------------------
    start = time.time()
    assert dtype in ['lev', 'ham', 'blast']
    if None in refpair:  # assume original value is for both chains together
        max_n_print /= 2
    distances = []
    for iclust, atn_pair in enumerate(antn_pairs):
        if [l['unique_ids'] for l in atn_pair if l is not None] == [l['unique_ids'] for l in refpair if l is not None]:
            continue
        h_atn, l_atn = atn_pair
        distances.append({'i' : iclust, 'h_ids' : gids(h_atn), 'l_ids' : gids(l_atn)})
        if dtype in ['lev', 'ham']:
            distances[-1]['dist'] = (lev_dist if dtype=='lev' else hamming_distance)(nseq(atn_pair), nseq(refpair))
        else:
            distances[-1]['seq'] = nseq(atn_pair)

    if dtype in ['blast']:
        qfos = [{'name' : 'query', 'seq' : nseq(refpair)}]
        tfos = [{'name' : 'iclust-%d'%dfo['i'], 'seq' : dfo['seq']} for dfo in distances]
        if workdir is None:
            workdir = choose_random_subdir('/tmp/%s/non-clonal-clusters'%os.getenv('USER', default='partis-work'))
        matchfos, mstats = run_blastn(qfos, tfos, workdir, aa=aa) #, debug=True, print_all_matches=True)
        mfo = matchfos['query']
        min_len = int(0.9 * max(mfo['alens']))  # require that the matches are at least 90% as long as the longest match (yes, ick)
        matchdict = {mfo['targets'][i] : int(mfo['mismatches'][i]) for i in range(len(mfo['targets'])) if mfo['alens'][i] > min_len}
        for dfo in distances:
            dfo['dist'] = matchdict.get('iclust-%d'%dfo['i'], 999999)

    if len(distances) == 0:
        return []
    sdists = sorted(distances, key=lambda d: d['dist'])
    if debug:
        nearest = sdists[0]
        near_dfos = [d for d in sdists if d['dist'] <= max_print_dist]
        if labelstr!='':
            labelstr = ' %s ' % labelstr
        print('%s%snearest: %d edit%s (%d cluster%s less than %d)' % (extra_str, labelstr, nearest['dist'], plural(nearest['dist']), len(near_dfos), plural(len(near_dfos)), max_print_dist))
        if len(near_dfos) > 0:
            print('    %s%s-dist  iclust   N uids (%s)' % (extra_str, dtype, ('%s only'%('h' if refpair[1] is None else 'l')) if None in refpair else 'h l'))
            for dfo in near_dfos[:max_n_print]:
                # assert len(dfo['h_ids']) == len(dfo['l_ids'])
                print('   %s   %3d     %3d     %3d %3d    %s' % (extra_str, dfo['dist'], dfo['i'], len(dfo['h_ids']), len(dfo['l_ids']), color_mutants(nseq(refpair), nseq(antn_pairs[dfo['i']]), amino_acid=aa, align_if_necessary=True)))
            if len(near_dfos) > max_n_print:
                print('                  (only printing nearest %d)' % max_n_print)
            print('')

    if debug:
        print('      non-clonal cluster time: %.1f' % (time.time() - start))
    return sdists

# ----------------------------------------------------------------------------------------
def subset_sequences(line, restrict_to_region=None, exclusion_3p=None, iseq=None):
    # NOTE don't call with <iseq> directly, instead use subset_iseq() below

    naive_seq = line['naive_seq']  # NOTE this includes the fv and jf insertions
    if iseq is None:
        muted_seqs = copy.deepcopy(line['seqs'])
    else:
        muted_seqs = [line['seqs'][iseq]]

    if restrict_to_region != '':  # NOTE this is very similar to code in performanceplotter. I should eventually cut it out of there and combine them, but I'm nervous a.t.m. because of all the complications there of having the true *and* inferred sequences so I'm punting
        if restrict_to_region in regions:
            bounds = line['regional_bounds'][restrict_to_region]
        elif restrict_to_region == 'cdr3':
            bounds = (line['codon_positions']['v'], line['codon_positions']['j'] + 3)
        else:
            assert False
        if exclusion_3p is not None:  # see NOTE in performanceplotter.hamming_to_true_naive()
            bounds = (bounds[0], max(bounds[0], bounds[1] - exclusion_3p))
        naive_seq = naive_seq[bounds[0] : bounds[1]]
        muted_seqs = [mseq[bounds[0] : bounds[1]] for mseq in muted_seqs]

    return naive_seq, muted_seqs

# ----------------------------------------------------------------------------------------
def subset_iseq(line, iseq, restrict_to_region=None, exclusion_3p=None):
    naive_seq, muted_seqs = subset_sequences(line, iseq=iseq, restrict_to_region=restrict_to_region, exclusion_3p=exclusion_3p)
    return naive_seq, muted_seqs[0]  # if <iseq> is specified, it's the only one in the list (this is kind of confusing, but it's less wasteful)

# ----------------------------------------------------------------------------------------
def get_n_muted(line, iseq, restrict_to_region='', return_mutated_positions=False):
    naive_seq, muted_seq = subset_iseq(line, iseq, restrict_to_region=restrict_to_region)
    return hamming_distance(naive_seq, muted_seq, return_mutated_positions=return_mutated_positions)

# ----------------------------------------------------------------------------------------
def get_mutation_rate(line, iseq, restrict_to_region=''):
    naive_seq, muted_seq = subset_iseq(line, iseq, restrict_to_region=restrict_to_region)
    return hamming_fraction(naive_seq, muted_seq)

# ----------------------------------------------------------------------------------------
def get_mutation_rate_and_n_muted(line, iseq, restrict_to_region='', exclusion_3p=None):
    naive_seq, muted_seq = subset_iseq(line, iseq, restrict_to_region=restrict_to_region, exclusion_3p=exclusion_3p)
    fraction, distance = hamming_fraction(naive_seq, muted_seq, also_return_distance=True)
    return fraction, distance

# ----------------------------------------------------------------------------------------
def get_sfs_occurence_info(line, restrict_to_region=None, debug=False):
    if restrict_to_region is None:
        naive_seq, muted_seqs = line['naive_seq'], line['seqs']  # I could just call subset_sequences() to get this, but this is a little faster since we know we don't need the copy.deepcopy()
    else:
        naive_seq, muted_seqs = subset_sequences(line, restrict_to_region=restrict_to_region)
    if debug:
        print('  %d %ssequences' % (len(muted_seqs), '%s region ' % restrict_to_region if restrict_to_region is not None else ''))
    mutated_positions = [hamming_distance(naive_seq, mseq, return_mutated_positions=True)[1] for mseq in muted_seqs]
    all_positions = sorted(set([p for mp in mutated_positions for p in mp]))
    if debug:
        print('    %.2f mean mutations  %s' % (numpy.mean([len(mpositions) for mpositions in mutated_positions]), ' '.join([str(len(mpositions)) for mpositions in mutated_positions])))
        print('    %d positions are mutated in at least one sequence' % len(all_positions))
    occurence_indices = [[i for i in range(len(line['unique_ids'])) if p in mutated_positions[i]] for p in all_positions]  # for each position in <all_positions>, a list of the sequence indices that have a mutation at that position
    occurence_fractions = [len(iocc) / float(len(line['unique_ids'])) for iocc in occurence_indices]  # fraction of all sequences that have a mutation at each position in <all_positions>
    return occurence_indices, occurence_fractions

# ----------------------------------------------------------------------------------------
def fay_wu_h(line, restrict_to_region=None, occurence_indices=None, n_seqs=None, debug=False):  # from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1461156/pdf/10880498.pdf and https://www.biorxiv.org/content/biorxiv/early/2017/10/19/145052.full.pdf
    if occurence_indices is None:
        occurence_indices, _ = get_sfs_occurence_info(line, restrict_to_region=restrict_to_region, debug=debug)
        n_seqs = len(line['unique_ids'])
    else:
        assert line is None  # don't pass both of 'em
    if n_seqs == 1:
        return 0.
    mutation_multiplicities = [len(oindices) for oindices in occurence_indices]  # <oindices> is a list of the indices of sequences that had this mutation, so this gives the number of sequences that had a mutation at this position
    theta_h = 0.
    for inm in range(1, n_seqs):
        theta_h += mutation_multiplicities.count(inm) * inm * inm
    theta_h *= 2. / (n_seqs * (n_seqs - 1))

    theta_pi = 0.
    for inm in range(1, n_seqs):
        theta_pi += mutation_multiplicities.count(inm) * inm * (n_seqs - inm)
    theta_pi *= 2. / (n_seqs * (n_seqs - 1))

    if debug:
        print('   h for %d seqs:  %6.2f - %6.2f = %6.2f' % (n_seqs, theta_pi, theta_h, theta_pi - theta_h))

    return theta_pi - theta_h

# ----------------------------------------------------------------------------------------
def dot_product(naive_seq, seq1, seq2):
    _, imutes1 = hamming_distance(naive_seq, seq1, return_mutated_positions=True)
    _, imutes2 = hamming_distance(naive_seq, seq2, return_mutated_positions=True)
    both_muted = set(imutes1) & set(imutes2)
    both_muted_to_same_thing = [imut for imut in both_muted if seq1[imut] == seq2[imut]]
    dot_product = len(both_muted_to_same_thing)
    # print '    naive  %s' % naive_seq
    # print '           %s' % utils.color_mutants(naive_seq, seq1)
    # print '           %s' % utils.color_mutants(naive_seq, seq2)
    # print '    dot %d' % dot_product
    return dot_product

# ----------------------------------------------------------------------------------------
def round_to_n_digits(val, n_digits):  # round <val> to <n_digits> significant figures
    if val == 0:
        return val
    return round(val, n_digits - int(math.floor(math.log10(abs(val)))) - 1)

# ----------------------------------------------------------------------------------------
def get_key(names):
    """
    Return a hashable combination of the two query names that's the same if we reverse their order.
    """
    return '.'.join(sorted([str(name) for name in names]))

# ----------------------------------------------------------------------------------------
def split_key(key):
    """
    Reverse the action of get_key().
    NOTE does not necessarily give a_ and b_ in the same order, though
    NOTE also that b_name may not be the same (if 0), and this just returns strings, even if original names were ints
    """
    # assert len(re.findall('.', key)) == 1  # make sure none of the keys had a dot in it
    return key.split('.')

# ----------------------------------------------------------------------------------------
def mkdir(path, isfile=False):  # adding this very late, so could use it in a lot of places
    if isfile:
        path = os.path.dirname(path)
    if path == '':  # if it's a relative path we give up
        return
    if not os.path.exists(path):
        os.makedirs(path)

# ----------------------------------------------------------------------------------------
def makelink(odir, target, link_name, dryrun=False, extra_str='', debug=False):  # <odir> is generally os.path.dirname(link_name)   for search: def link def ln
    if os.path.exists(link_name):
        if os.path.islink(link_name):
            if os.path.isdir(link_name):  # if the link we're making (well, made on a previous call) is to a dir, calling a second time will instead create a new link inside the linked-to dir (which we don't want), so we have to remove it
                os.remove(link_name)
        else:
            raise Exception('link name %s exists and is not link' % link_name)
    target = target.replace(odir+'/', '')
    link_name = link_name.replace(odir+'/', '')
    if '/' in link_name:  # it's a subdir of odir, we need to add at least one ../
        n_slashes = [x[0] for x in itertools.groupby(link_name)].count('/')  # have to collapse any adjacent /s
        target = n_slashes * '../' + target
    mkdir(odir)
    simplerun('cd %s && ln -sf %s %s' % (odir, target, link_name), shell=True, dryrun=dryrun, extra_str=extra_str, debug=debug)

    if not os.path.exists(target if target==fpath(target) else odir+'/'+target):
        raise Exception('linked to missing file in odir %s (target %s)' % (odir, target if target==fpath(target) else odir+'/'+target))

# ----------------------------------------------------------------------------------------
def fpath(path):  # if <path> is relative, add full path from root (path can also be None)
    if path is None or path[0] == '/':
        return path
    else:
        return '%s/%s' % (os.getcwd(), path)

# ----------------------------------------------------------------------------------------
def prep_dir(dirname, wildlings=None, subdirs=None, rm_subdirs=False, fname=None, allow_other_files=False):
    """
    Make <dirname> if it d.n.e.
    Also, if shell glob <wildling> is specified, remove existing files which are thereby matched.
    """
    if fname is not None:  # passed in a file name, and we want to prep the file's dir
        assert dirname is None
        dirname = os.path.dirname(fname)
        if dirname == '' or dirname[0] != '/':
            dirname = '/'.join([pn for pn in [os.getcwd(), dirname] if pn != ''])

    if wildlings is None:
        wildlings = []
    elif isinstance(wildlings, six.string_types):  # allow to pass in just a string, instead of a list of strings
        wildlings = [wildlings, ]

    if subdirs is not None:  # clean out the subdirs first
        for subdir in subdirs:
            prep_dir(dirname + '/' + subdir, wildlings=wildlings, allow_other_files=allow_other_files)
            if rm_subdirs:
                os.rmdir(dirname + '/' + subdir)

    if os.path.exists(dirname):
        for wild in wildlings:
            for ftmp in glob.glob(dirname + '/' + wild):
                if os.path.exists(ftmp):
                    os.remove(ftmp)
                else:
                    print('%s file %s exists but then it doesn\'t' % (color('red', 'wtf'), ftmp))
        remaining_files = [fn for fn in os.listdir(dirname) if subdirs is None or fn not in subdirs]  # allow subdirs to still be present
        if len(remaining_files) > 0 and not allow_other_files:  # make sure there's no other files in the dir
            raise Exception('files remain in %s despite wildlings %s:\n (%s)' % (dirname, wildlings, ' '.join(remaining_files)))
    else:
        os.makedirs(dirname)

# ----------------------------------------------------------------------------------------
def clean_files(fnames, expect_missing=False):  # <fnames> can include dirs, we first sort it to put them after the files
    if len(fnames) != len(set(fnames)):  # remove any duplicates (don't always do it, since we'd rather not change the order if we don't need to, although we sort for files before removing anything below so it doesn't matter too much)
        fnames = list(set(fnames))
    fnames = sorted(fnames, key=len, reverse=True)  # also sort by len, so subdirs are removed before their parents
    fnames = sorted(fnames, key=lambda x: os.path.isdir(x))  # put files first, dirs second (doing this sort second means they're sorted for length *within* file/dir classes)
    missing_files = []
    for fn in fnames:
        if os.path.isfile(fn) or os.path.islink(fn):
            os.remove(fn)
        elif os.path.isdir(fn):
            if len(os.listdir(fn)) > 0:
                raise Exception('can\'t rmdir, files remain in %s (%s)' % (fn, ' '.join(os.listdir(fn))))
            os.rmdir(fn)
        elif not os.path.exists(fn):
            missing_files.append(fn)
        else:  # not sure if there's another possibility
            assert False
    if len(missing_files) > 0 and not expect_missing:
        basedirs = sorted(set(os.path.dirname(f) for f in fnames), key=lambda x: x.count('/'))  # sort by number of slashes, i.e. most-parent dir is first
        basedirs = [basedirs[0]] + [d for d in basedirs if basedirs[0] not in d]  # hackey way to remove any that are just subdirs (if all the files aren't under one parent dir, these shenanigans won't really work)
        print('      %s expected to remove %d/%d files+dirs that weren\'t in %s/: %s' % (color('yellow', 'warning'), len(missing_files), len(fnames), ' '.join(basedirs), ' '.join(f.replace(basedirs[0] + '/', '') for f in missing_files)))

# ----------------------------------------------------------------------------------------
def rmdir(dname, fnames=None):
    if fnames is not None:
        clean_files(fnames)
    remaining_files = os.listdir(dname)
    if len(remaining_files) > 0:
        raise Exception('files remain in %s so can\'t rm dir: %s' % (dname, ' '.join(remaining_files)))
    os.rmdir(dname)

# ----------------------------------------------------------------------------------------
def process_input_line(info, skip_literal_eval=False):
    from . import indelutils
    """
    Attempt to convert all the keys and values in <info> according to the specifications in <io_column_configs> (e.g. splitting lists, casting to int/float, etc).
    """

    if 'v_gene' in info and info['v_gene'] == '':
        return

    if 'seq' in info:  # old simulation files
        for key in ['unique_id', 'seq', 'indelfo']:
            if key not in info:
                continue
            info[key + 's'] = info[key]
            del info[key]
        if 'indelfos' not in info:  # hm, at least some old sim files don't have 'indelfo'
            info['indelfos'] = str(indelutils.get_empty_indel())
        info['indelfos'] = '[' + info['indelfos'] + ']'
        info['input_seqs'] = info['seqs']

    for key in info:
        if info[key] == '':  # handle these below, once we know how many seqs in the line
            continue
        convert_fcn = conversion_fcns.get(key, pass_fcn)
        if skip_literal_eval and convert_fcn is ast.literal_eval:  # it's really slow (compared to the other conversions at least), and it's only for keys that we hardly ever use
            continue
        if key in io_column_configs['lists-of-lists']:
            info[key] = convert_fcn(info[key].split(';'))
        elif key in io_column_configs['lists']:
            info[key] = [convert_fcn(val) for val in info[key].split(':')]
        else:
            info[key] = convert_fcn(info[key])

    # this column is called 'seqs' internally (for conciseness and to avoid rewriting a ton of stuff) but is called 'indel_reversed_seqs' in the output file to avoid user confusion
    if 'indel_reversed_seqs' in info and 'input_seqs' in info:  # new-style csv output and simulation files, i.e. it stores 'indel_reversed_seqs' instead of 'seqs'
        if info['indel_reversed_seqs'] == '':
            info['indel_reversed_seqs'] = ['' for _ in range(len(info['unique_ids']))]
        transfer_indel_reversed_seqs(info)
    elif 'seqs' in info:  # old-style csv output file: just copy 'em into the explicit name
        info['indel_reversed_seqs'] = info['seqs']

    # process things for which we first want to know the number of seqs in the line
    for key in [k for k in info if info[k] == '']:
        if key in io_column_configs['lists']:
            info[key] = ['' for _ in range(len(info['unique_ids']))]
        elif key in io_column_configs['lists-of-lists']:
            info[key] = [[] for _ in range(len(info['unique_ids']))]

    # NOTE indels get fixed up (espeicially/only for old-style files) in add_implicit_info(), since we want to use the implicit info to do it

    if 'all_matches' in info and isinstance(info['all_matches'], dict):  # it used to be per-family, but then I realized it should be per-sequence, so any old cache files lying around have it as per-family
        info['all_matches'] = [info['all_matches']]
    # make sure everybody's the same lengths
    for key in [k for k in info if k in io_column_configs['lists']]:
        if key == 'duplicates' and len(info[key]) != len(info['unique_ids']) and len(info[key]) == 1:  # fix problem caused by 'duplicates' being in both 'lists' and 'lists-of-lists', combined with get_line_for_output() and this fcn having the if/else blocks in opposite order (fixed now, but some files are probably lying around with the whacked out formatting)
            info[key] = info[key][0]
            info[key] = [ast.literal_eval(v) for v in info[key]]  # specifically, in get_line_for_output() 'lists' was first in the if/else block, so the list of duplicates for each sequence was str() converted rather than being converted to colon/semicolon separated string
        if len(info[key]) != len(info['unique_ids']):
            raise Exception('list length %d for %s not the same as for unique_ids %d\n  uids: %s\n  %s: %s' % (len(info[key]), key, len(info['unique_ids']), info['unique_ids'], key, info[key]))

# ----------------------------------------------------------------------------------------
def revcomp(nuc_seq):
    if 'Bio.Seq' not in sys.modules:  # import is frequently slow af
        from Bio.Seq import Seq
    bseq = sys.modules['Bio.Seq']
    return str(bseq.Seq(nuc_seq).reverse_complement())

# ----------------------------------------------------------------------------------------
# NOTE do *not* call this directly on seqs from annotations, since they're not necessarily padded to the start of a codon (i.e. without first calling pad_seq_for_translation())
def ltranslate(nuc_seq, trim=False, debug=False):  # local file translation function
    if 'Bio.Seq' not in sys.modules:  # import is frequently slow af
        from Bio.Seq import Seq
    bseq = sys.modules['Bio.Seq']
    if trim:  # this should probably be the default, but i don't want to change anything that's using the padding (even though it probably wouldn't matter)
        nuc_seq = trim_nuc_seq(nuc_seq.strip(ambig_base))
    final_nseq = pad_nuc_seq(nuc_seq)
    aa_seq = str(bseq.Seq(final_nseq).translate())  # the padding is annoying, but it's extremely common for bcr sequences to have lengths not a multiple of three (e.g. because out out of frame rearrangements), so easier to just always check for it
    if debug:
        assert len(final_nseq) == len(aa_seq*3)
        print(' '.join([final_nseq[ic : ic + 3] for ic in range(0, len(final_nseq), 3)]))
        print('   '.join(aa_seq))
    return aa_seq

# ----------------------------------------------------------------------------------------
def get_cdr3_seq(info, iseq):  # NOTE includeds both codons, i.e. not the same as imgt definition
    return info['seqs'][iseq][info['codon_positions']['v'] : info['codon_positions']['j'] + 3]

# ----------------------------------------------------------------------------------------
def add_naive_seq_aa(line):
    if 'naive_seq_aa' in line:
        return
    line['naive_seq_aa'] = ltranslate(pad_seq_for_translation(line, line['naive_seq']))

# ----------------------------------------------------------------------------------------
# return naive sequence with uncertain bases (e.g. non-templated insertions) replace with Ns
# n_fuzz: number of bases by which to expand non-templated insertions (set to negative value to mask entire cdr3 [excluding conserved codons])
def get_cdr3_masked_naive(antn, n_fuzz=1, debug=True):
    if antn.get('is_fake_paired', False):
        raise Exception('can\'t mask fake paired annotation')
    masked_seq = antn['naive_seq']
    if n_fuzz < 0:
        cdr3_start, cdr3_end = [antn['codon_positions'][r] for r in 'vj']
        masked_seq = masked_seq[ : cdr3_start + 3] + ambig_base * (cdr3_end - cdr3_start - 3) + masked_seq[cdr3_end : ]
    else:
        for tbound in boundaries:
            nti_start, nti_end = [antn['regional_bounds'][r][i] for r, i in zip(tbound, (1, 0))]  # (1, 0) takes the end of first region but start of second, e.g. for vd boundary we want end of v, start of d
            nti_start -= n_fuzz
            nti_end += n_fuzz
            masked_seq = masked_seq[ : nti_start] + ambig_base * (nti_end - nti_start) + masked_seq[nti_end : ]
    if debug:
        print('    masked %s:' % ('entire cdr3, excluding conserved codons' if n_fuzz < 0 else 'non-templated insertions with fuzz %d'%n_fuzz))
        # print_reco_event(antn)
        print(color_mutants(antn['naive_seq'], masked_seq, extra_str='        '))

    return masked_seq

# ----------------------------------------------------------------------------------------
def pad_seq_for_translation(line, tseq, return_n_padded=False, debug=False):  # this duplicates the arithmetic in waterer that pads things to the same length, but we do that after a bunch of places where we might call this fcn, so we need to check for it here as well
    # NOTE this duplicates [is reverse of?] code in n_variable_ambig() (but not really sure how to clean up/combine them)
    # NOTE also duplicates some of get_codon_list()
    old_tseq = tseq  # just for dbg
    fv_xtra, v_5p_xtra = 0, 0
    if len(line['fv_insertion']) % 3 != 0:  # first base in v (i.e. after any fv insertion) is in frame, so have to align with that
        fv_xtra = 3 - len(line['fv_insertion']) % 3  # e.g. if there's 1 fv insertion base, add 2 Ns to left side
        tseq = fv_xtra * ambig_base + tseq
    if line['v_5p_del'] % 3 != 0:
        v_5p_xtra = line['v_5p_del'] % 3  # e.g. if there's 1 v 5p deleted base, add 1 N to left side
        tseq = v_5p_xtra * ambig_base + tseq
    end_amb = ''
    if len(tseq) % 3 != 0:  # have to pad also the righthand side to allow for smashing h+l seqs together (or at least it's easier to do it here than elsewhere)
        end_amb = (3 - len(tseq) % 3) * ambig_base
        tseq += end_amb
    assert len(tseq) % 3 == 0
    if debug:
        print('  fv: 3 - %d%%3: %d  v_5p: %d%%3: %d  end amb: %d  (jf %d)' % (len(line['fv_insertion']), fv_xtra, line['v_5p_del'], v_5p_xtra, len(end_amb), len(line.get('jf_insertion', []))))  # NOTE the first one is kind of wrong, since it's 0 if the %3 is 0
        print('    %s%s%s%s' % (color('blue', fv_xtra * ambig_base), color('blue', v_5p_xtra * ambig_base), color_mutants(old_tseq, old_tseq), color('blue', end_amb)))
    if return_n_padded:
        return tseq, (fv_xtra + v_5p_xtra, len(end_amb))
    else:
        return tseq

# ----------------------------------------------------------------------------------------
def add_seqs_aa(line, debug=False):  # NOTE similarity to block in add_extra_column()
    from . import indelutils
    if 'seqs_aa' in line and line['seqs_aa'].count(None) == 0:  # if we've just added some seqs to the <line> some will have None aa seqs
        return
    line['seqs_aa'] = [ltranslate(pad_seq_for_translation(line, s)) for s in line['seqs']]
    if 'input_seqs' in line:
        line['input_seqs_aa'] = [ltranslate(pad_seq_for_translation(line, inseq)) if indelutils.has_indels_line(line, iseq) else irseq_aa for iseq, (inseq, irseq_aa) in enumerate(zip(line['input_seqs'], line['seqs_aa']))]
    add_naive_seq_aa(line)
    if debug:
        print(pad_lines('\n'.join(line['seqs_aa'])))
        if any(indelutils.has_indels_line(line, i) for i in range(len(line['unique_ids']))):
            print('%s in add_seqs_aa() one has an indel, input_seqs_aa debug needs to be implemented' % color('red', 'error'))

# ----------------------------------------------------------------------------------------
def shm_aa(line, iseq=None, uid=None):  # it's kind of weird to have this fcn separate, whereas the non-aa one we don't, but it only really exists so it can add the aa seqs UPDATE ok using it elsewhere now
    if iseq is None:  # have to specify either iseq or uid (no i don't care if you specified both and they're inconsistent
        iseq = line['unique_ids'].index(uid)
    add_naive_seq_aa(line)
    add_seqs_aa(line)
    return hamming_distance(line['naive_seq_aa'], line['seqs_aa'][iseq], amino_acid=True)

# ----------------------------------------------------------------------------------------
def pad_nuc_seq(nseq, side='right', return_n_padded=False):  # if length not multiple of three, pad on right (by default) with Ns
    padstr = ''
    if len(nseq) % 3 != 0:
        if side == 'right':
            padstr = 'N' * (3 - (len(nseq) % 3))
            nseq += padstr
        elif side == 'left':
            padstr = 'N' * (3 - (len(nseq) % 3))
            nseq = padstr + nseq
        else:
            assert False
    return (nseq, len(padstr)) if return_n_padded else nseq

# ----------------------------------------------------------------------------------------
# duplicates some code in waterer.get_padding_parameters() (NOTE the version there doesn't include the jf insertion length)
def get_pad_parameters(line, glfo, iseq=0):  # iseq value shouldn't matter, but have it here just to make it explicit
    fvstuff = max(0, len(line['fv_insertion']) - line['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
    gl_cpos = glfo['cyst-positions'][line['v_gene']] + fvstuff
    cpos = line['codon_positions']['v']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
    gl_cpos_to_j_end = len(line['seqs'][iseq]) - cpos + line['j_3p_del'] - len(line['jf_insertion'])
    return gl_cpos, gl_cpos_to_j_end  # when using these return values, call them (ldist, rdist)

# ----------------------------------------------------------------------------------------
# this fcn actually does the re-padding, whereas re_pad_hmm_seqs() figures out what padding is needed in orderto make the hmm antn match sw, then calls this fcn
def re_pad_atn(n_fv_pad, n_jf_pad, atn, glfo, extra_str='      ', debug=False):  # NOTE both N pads can be negative (yeah, probably should rename them now)
    from . import indelutils
    remove_all_implicit_info(atn)
    padstrs, n_trims = {'fv' : '', 'jf' : ''}, {'fv' : 0, 'jf' : 0}
    for istr, n_pad in [['fv', n_fv_pad], ['jf', n_jf_pad]]:
        if n_pad == 0:
            continue
        elif n_pad > 0:  # actually pad with Ns
            padstrs[istr] = n_pad * ambig_base
            atn[istr+'_insertion'] += padstrs[istr]
        else:  # trim it
            if abs(n_pad) > len(atn[istr+'_insertion']):  # don't trim off more than there is
                print('    %s tried to trim %s %d more than there was %d when re-padding annotation for %s' % (wrnstr(), istr, n_pad, len(atn[istr+'_insertion']), ':'.join(atn['unique_ids'])))
                n_pad = -len(atn[istr+'_insertion'])
            n_trims[istr] = abs(n_pad)
            atn[istr+'_insertion'] = atn[istr+'_insertion'][ : len(atn[istr+'_insertion']) - n_trims[istr]]
    for seqkey in ['seqs', 'input_seqs']:
        for iseq in range(len(atn['unique_ids'])):  # only 2 of the 4 numbers should be set, but I don't think i really need to check that
            atn[seqkey][iseq] = atn[seqkey][iseq][n_trims['fv'] : len(atn[seqkey][iseq]) - n_trims['jf']]
            atn[seqkey][iseq] = padstrs['fv'] + atn[seqkey][iseq] + padstrs['jf']
    for iseq in [i for i in range(len(atn['unique_ids'])) if indelutils.has_indels_line(atn, i)]:
        if len(padstrs['fv']) > 0 or len(padstrs['jf']) > 0:
            indelutils.pad_indelfo(atn['indelfos'][iseq], padstrs['fv'], padstrs['jf'])
        if any(n>0 for n in n_trims.values()):
            indelutils.trim_padding_from_indelfo(atn['indelfos'][iseq], n_trims['fv'], n_trims['jf'])
    add_implicit_info(glfo, atn)  # NOTE this is probably slow, i could probably save a lot of time by not recalculating e.g. n_mutations
    if debug:
        if len(padstrs['fv']) > 0 or len(padstrs['jf']) > 0:
            print('%spadded with leftstr: %s rightstr: %s' % (extra_str, padstrs['fv'], padstrs['jf']))
        if any(n>0 for n in n_trims.values()):
            print('%strimmed with fv: %d jf: %d' % (extra_str, n_trims['fv'], n_trims['jf']))

# ----------------------------------------------------------------------------------------
# In subset-partition, we end up with hmm-inferred multi-seq annotations that were padded in their individual subset process, i.e. that in some cases aren't padded enough for the full sample.
# Here we add any extra padding in the sw info (which will have been padded on the full sample when reading the subset-merged sw cache file)
# UPDATE: actually also need to *trim* padding in input annotations, at least for the case we're subset partitioning and ignoring small clusters, in which case sw padding in the merge process will be missing all the sequences from small clusters, so it can end up with less padding than the original sw (and thus input annotation)
# NOTE that in most cases this just modifies elements in <input_antn_list>, but sometimes also modifies the list, so in general you need to e.g. remake any associated annotation dict from the calling fcn
def re_pad_hmm_seqs(input_antn_list, input_glfo, sw_info, debug=False):  # NOTE quite similar to pad_seqs_to_same_length() in waterer.py (although there we change a lot more stuff by hand, i think because i didn't yet have remove_all_implicit_info [or maybe bc that would be slower])
    from . import indelutils
    # ----------------------------------------------------------------------------------------
    def has_bad_lengths(iatn):
        return any(len(sw_info[u]['naive_seq']) != len(iatn['naive_seq']) for u in iatn['unique_ids'])
    # ----------------------------------------------------------------------------------------
    def print_aligned_seqs(iatn, xstr='', extra_str='        '):
        print('   %snaive seq length in input annotation %d different from that in sw info %s' % (xstr, len(iatn['naive_seq']), ' '.join(str(len(sw_info[u]['naive_seq'])) for u in iatn['unique_ids'])))
        align_many_seqs([{'name' : 'input-atn', 'seq' : iatn['naive_seq']}] + [{'name' : u, 'seq' : sw_info[u]['naive_seq']} for u in iatn['unique_ids']], extra_str=extra_str, debug=True)
    # ----------------------------------------------------------------------------------------
    if debug:
        print('  re-padding %d annotations' % len(input_antn_list))
    n_padded = 0
    for ia, iatn in enumerate(input_antn_list):
        if not has_bad_lengths(iatn):
            continue
        if debug:
            print_aligned_seqs(iatn)
        sw_ldists, sw_rdists = zip(*[get_pad_parameters(sw_info[iatn['unique_ids'][i]], input_glfo) for i, u in enumerate(iatn['unique_ids'])])
        sw_ldist, sw_rdist = [get_single_entry(list(set(l))) for l in [sw_ldists, sw_rdists]]
        ia_ldist, ia_rdist = get_pad_parameters(iatn, input_glfo)  # they should all be the same, so can use 0 (?)
        leftpad, rightpad = sw_ldist - ia_ldist, sw_rdist - ia_rdist
        if leftpad == 0 and rpad == 0:
            print('    %s lengths don\'t match when re-padding hmm seqs, but padding parameters are the same (will probably crash just below)' % wrnstr())
        else:
            re_pad_atn(leftpad, rightpad, iatn, input_glfo, debug=debug)
            n_padded += 1
        if has_bad_lengths(iatn):
            print_aligned_seqs(iatn, xstr='    failed to fix: ', extra_str='            ')
            print('        %s didn\'t manage to fix lengths for %s (see previous lines)' % (wrnstr(), ':'.join(iatn['unique_ids'])))
            print('           so giving up and making multi-seq line from sw info')
            input_antn_list[ia] = synthesize_multi_seq_line_from_reco_info(iatn['unique_ids'], sw_info)  # this is the only place we modify the annotation list (rather than just modifying annotations)
    if debug:
        print('    fixed padding for %d / %d annotations' % (n_padded, len(input_antn_list)))

# ----------------------------------------------------------------------------------------
def trim_nuc_seq(nseq):  # if length not multiple of three, trim extras from the right side
    if len(nseq) % 3 != 0:
        nseq = nseq[ : len(nseq) - (len(nseq) % 3)]
    return nseq

# ----------------------------------------------------------------------------------------
def get_frame(line):  # return frame of (start of) V gene, e.g. for input to gctree
    return (line['v_5p_del'] - len(line['fv_insertion'])) % 3 + 1

# ----------------------------------------------------------------------------------------
def add_extra_column(key, info, outfo, glfo=None, definitely_add_all_columns_for_csv=False):
    from . import treeutils
    from . import indelutils
    if outfo is None:  # hacking on this behavior that doesn't require you to pass in <outfo>
        outfo = {}
    # NOTE use <info> to calculate all quantities, *then* put them in <outfo>: <outfo> only has the stuff that'll get written to the file, so can be missing things that are needed for calculations
    if key == 'cdr3_seqs':
        outfo[key] = [get_cdr3_seq(info, iseq) for iseq in range(len(info['unique_ids']))]
    if key == 'cdr3_seqs_aa':
        outfo[key] = [ltranslate(get_cdr3_seq(info, iseq)) for iseq in range(len(info['unique_ids']))]
    elif key == 'full_coding_naive_seq':
        if glfo is None:
            raise Exception('have to pass in glfo for extra annotation column \'%s\'' % key)
        delstr_5p = glfo['seqs']['v'][info['v_gene']][ : info['v_5p_del']]  # bit missing from the input sequence
        outfo[key] = delstr_5p + info['naive_seq']
        if info['j_3p_del'] > 0:
            delstr_3p = glfo['seqs']['j'][info['j_gene']][-info['j_3p_del'] : ]
            outfo[key] += delstr_3p
        # print outfo['unique_ids']
        # color_mutants(info['naive_seq'], outfo[key], print_result=True, align=True, extra_str='  ')
    elif key == 'full_coding_input_seqs':
        def inpseq(i): return info['input_seqs'][iseq][len(info['fv_insertion']) : len(info['input_seqs'][iseq]) - len(info['jf_insertion'])]
        full_coding_input_seqs = [info['v_5p_del'] * ambig_base + inpseq(iseq) + info['j_3p_del'] * ambig_base for iseq in range(len(info['unique_ids']))]
        outfo[key] = full_coding_input_seqs
        # for iseq in range(len(info['unique_ids'])):
        #     print info['unique_ids'][iseq]
        #     color_mutants(info['input_seqs'][iseq], full_coding_input_seqs[iseq], print_result=True, align=True, extra_str='  ')
    elif key == 'consensus_seq':
        outfo[key] = cons_seq_of_line(info)
    elif key == 'consensus_seq_aa':
        add_seqs_aa(info)
        outfo[key] = cons_seq_of_line(info, aa=True)
    elif key == 'naive_seq_aa':  # NOTE similarity to add_naive_seq_aa()
        add_naive_seq_aa(info)
        outfo[key] = info['naive_seq_aa']
    elif key == 'seqs_aa':  # NOTE similarity to add_seqs_aa()
        outfo[key] = [ltranslate(s) for s in info['seqs']]
    elif key in ['cons_dists_nuc', 'cons_dists_aa']:
        treeutils.add_cons_dists(info, aa='_aa' in key)
        outfo[key] = info[key]  # has to be done in two steps (see note at top of fcn)
    elif key in linekeys['hmm'] + linekeys['sw'] + linekeys['simu']:  # these are added elsewhere
        if definitely_add_all_columns_for_csv:
            if key in io_column_configs['lists']:
                outfo[key] = [None for _ in info['unique_ids']]
            elif key in io_column_configs['lists-of-lists']:
                outfo[key] = [[] for _ in info['unique_ids']]
            else:
                outfo[key] = None
        else:
            return  # only here to remind you that nothing happens
    elif key in list(input_metafile_keys.values()):  # uh, not really sure what's the best thing to do, but this only gets called on deprecated csv files, so oh well
        outfo[key] = [None for _ in info['unique_ids']]
    elif key == 'has_shm_indels':
        outfo[key] = [indelutils.has_indels_line(info, i) for i in range(len(info['unique_ids']))]
    else:  # this happens all the time now
        # raise Exception('column \'%s\' missing from annotation' % key)
        # print '    %s column \'%s\' missing from annotation' % (wrnstr(), key)
        pass
    return outfo.get(key)

# ----------------------------------------------------------------------------------------
def transfer_indel_reversed_seqs(line):
    # if there's no indels, we will have just written 'input_seqs' to the file, and left 'indel_reversed_seqs' empty
    line['seqs'] = [line['indel_reversed_seqs'][iseq] if line['indel_reversed_seqs'][iseq] != '' else line['input_seqs'][iseq] for iseq in range(len(line['unique_ids']))]

# ----------------------------------------------------------------------------------------
def transfer_indel_info(info, outfo):  # NOTE reverse of this happens in indelutils.deal_with_indel_stuff()
    from . import indelutils
    """
    for keys in special_indel_columns_for_output: in memory, I need the indel info under the 'indelfos' key (for historical reasons that I don't want to change a.t.m.), but I want to mask that complexity for output
    """
    if special_indel_columns_for_output[0] in info:  # they're only already transferred if we're reading simulation files for merging
        for sicfo in special_indel_columns_for_output:
            outfo[sicfo] = info[sicfo]
    else:
        if info['invalid']:
            return
        outfo['has_shm_indels'] = [indelutils.has_indels(ifo) for ifo in info['indelfos']]
        outfo['qr_gap_seqs'] = [ifo['qr_gap_seq'] for ifo in info['indelfos']]
        outfo['gl_gap_seqs'] = [ifo['gl_gap_seq'] for ifo in info['indelfos']]
        outfo['indel_reversed_seqs'] = ['' if not indelutils.has_indels(info['indelfos'][iseq]) else info['indel_reversed_seqs'][iseq] for iseq in range(len(info['unique_ids']))]  # if no indels, it's the same as 'input_seqs', so set indel_reversed_seqs to empty strings

# ----------------------------------------------------------------------------------------
def get_line_for_output(headers, info, glfo=None):
    """ Reverse the action of process_input_line() """
    # NOTE only used by (deprecated) csv writer now
    outfo = {}
    transfer_indel_info(info, outfo)
    for key in headers:
        if key in ['seqs', 'indelfo', 'indelfos']:
            continue

        if key not in special_indel_columns_for_output:  # these keys were already added to <outfo> (and are not in <info>), but everyone else needs to be transferred from <info> to <outfo>
            if key in info:
                if key in io_column_configs['lists-of-lists']:
                    outfo[key] = copy.deepcopy(info[key])
                else:
                    outfo[key] = info[key]
            else:
                add_extra_column(key, info, outfo, glfo=glfo, definitely_add_all_columns_for_csv=True)

        str_fcn = str
        if key in io_column_configs['floats']:
            str_fcn = repr  # keeps it from losing precision (we only care because we want it to match expectation if we read it back in)

        if key in io_column_configs['lists-of-lists']:
            if '_per_gene_support' in key:
                outfo[key] = list(outfo[key].items())
            for isl in range(len(outfo[key])):
                outfo[key][isl] = ':'.join([str_fcn(s) for s in outfo[key][isl]])
            outfo[key] = ';'.join(outfo[key])
        elif key in io_column_configs['lists']:
            outfo[key] = ':'.join([str_fcn(v) for v in outfo[key]])
        else:
            outfo[key] = str_fcn(outfo[key])

        # if key == 'tree':  # if this is commented, the newick tree strings get written with trailing newline, which looks weird when you look at the .csv by hand, but otherwise works just fine
        #     outfo[key] = outfo[key].strip()

    return outfo

# ----------------------------------------------------------------------------------------
def merge_simulation_files(outfname, file_list, headers, cleanup=True, n_total_expected=None, n_per_proc_expected=None, use_pyyaml=False, dont_write_git_info=False):
    if getsuffix(outfname) == '.csv':  # old way
        n_event_list, n_seq_list = merge_csvs(outfname, file_list, old_simulation=True, cleanup=True)
    elif getsuffix(outfname) == '.yaml':  # new way
        n_event_list, n_seq_list = merge_yamls(outfname, file_list, headers, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info, cleanup=True)
    else:
        raise Exception('unhandled annotation file suffix %s' % args.outfname)

    print('   read %d event%s with %d seqs from %d %s files' % (sum(n_event_list), plural(len(n_event_list)), sum(n_seq_list), len(file_list), getsuffix(outfname)))
    if n_total_expected is not None:
        if isinstance(n_per_proc_expected, list):  # different number for each proc
            if n_event_list != n_per_proc_expected:
                raise Exception('expected events per proc (%s), different from those read from the files (%s)' % (' '.join(str(n) for n in n_per_proc_expected), ' '.join([str(n) for n in n_event_list])))
        else:  # all procs the same number
            if n_event_list.count(n_per_proc_expected) != len(n_event_list):
                raise Exception('expected %d events per proc, but read: %s' % (n_per_proc_expected, ' '.join([str(n) for n in n_event_list])))
        if n_total_expected != sum(n_event_list):
            print('  %s expected %d total events but read %d (per-file couts: %s)' % (color('yellow', 'warning'), n_total_expected, sum(n_event_list), ' '.join([str(n) for n in n_event_list])))

# ----------------------------------------------------------------------------------------
def merge_csvs(outfname, csv_list, cleanup=False, old_simulation=False):
    header = None
    outfo = []
    if old_simulation:
        n_event_list, n_seq_list = [], []
    for infname in csv_list:
        if getsuffix(infname) not in ['.csv', '.tsv']:
            raise Exception('unhandled suffix, expected .csv or .tsv: %s' % infname)
        delimiter = ',' if getsuffix(infname) == '.csv' else '\t'
        with open(infname) as sub_outfile:
            reader = csv.DictReader(sub_outfile, delimiter=str(delimiter))
            header = reader.fieldnames
            if old_simulation:
                n_event_list.append(0)
                n_seq_list.append(0)
                last_reco_id = None
            for line in reader:
                outfo.append(line)
                if old_simulation:
                    n_seq_list[-1] += 1
                    if last_reco_id is None or line['reco_id'] != last_reco_id:
                        last_reco_id = line['reco_id']
                        n_event_list[-1] += 1
        if cleanup:
            os.remove(infname)
            os.rmdir(os.path.dirname(infname))

    mkdir(outfname, isfile=True)
    with open(outfname, csv_wmode()) as outfile:
        writer = csv.DictWriter(outfile, header, delimiter=str(delimiter))
        writer.writeheader()
        for line in outfo:
            writer.writerow(line)

    if old_simulation:
        return n_event_list, n_seq_list

# ----------------------------------------------------------------------------------------
def update_gene_names_in_line(line, name_mapping):
    """
    Update gene names in an annotation line using the provided name mapping.
    This is used when merging germline info to keep annotations in sync with the merged glfo.

    Args:
        line: annotation line (dict) with gene name references
        name_mapping: dict of {region: {dropped_name: retained_name}}

    Returns:
        None (modifies line in place)
    """
    for region in regions:
        if len(name_mapping[region]) == 0:  # skip regions with no name changes
            continue

        # Update main gene assignment field
        gene_key = region + '_gene'
        if line[gene_key] != '' and line[gene_key] in name_mapping[region]:
            line[gene_key] = name_mapping[region][line[gene_key]]

        # Update per-gene support dictionary (has gene names as keys)
        support_key = region + '_per_gene_support'
        if support_key in line and line[support_key] is not None:
            for old_name, new_name in name_mapping[region].items():
                if old_name in line[support_key]:
                    old_value = line[support_key].pop(old_name)
                    # If new name already exists, keep the one with higher support
                    if new_name in line[support_key]:
                        line[support_key][new_name] = max(line[support_key][new_name], old_value)
                    else:
                        line[support_key][new_name] = old_value

# ----------------------------------------------------------------------------------------
def merge_yamls(outfname, yaml_list, headers, cleanup=False, use_pyyaml=False, dont_write_git_info=False, remove_duplicates=False, return_merged_objects=False, debug=False):
    from . import glutils
    merged_annotation_list, merged_keys = [], set()
    merged_cpath, merged_glfo = None, None
    n_event_list, n_seq_list = [], []
    for infname in yaml_list:
        glfo, annotation_list, cpath = read_yaml_output(infname, dont_add_implicit_info=True)
        if debug:
            print('        %d sequences in %d clusters from %s' % (sum(len(l['unique_ids']) for l in annotation_list), len(annotation_list), infname))
        if remove_duplicates:  # NOTE this doesn't catch duplicates *within* each subfile, but atm I'm only worried about the case where they appear at most once in each subfile, so oh well
            annotation_list = [l for l in annotation_list if ':'.join(l['unique_ids']) not in merged_keys]
            for ptn in cpath.partitions:
                ptn = [c for c in ptn if ':'.join(c) not in merged_keys]
        n_event_list.append(len(annotation_list))
        n_seq_list.append(sum(len(l['unique_ids']) for l in annotation_list))
        if merged_cpath is None:
            merged_cpath = cpath
        else:
            assert len(cpath.partitions) == len(merged_cpath.partitions)
            assert cpath.i_best == merged_cpath.i_best  # not sure what to do otherwise (and a.t.m. i'm only using this  to merge simulation files, which only ever have one partition)
            for ip in range(len(cpath.partitions)):
                merged_cpath.partitions[ip] += cpath.partitions[ip]  # NOTE this assumes there's no overlap between files, e.g. if it's simulation and the files are totally separate
                merged_cpath.logprobs[ip] += cpath.logprobs[ip]  # they'll be 0 for simulation, but may as well handle it
                # NOTE i think i don't need to mess with these, but not totally sure: self.n_procs, self.ccfs, self.we_have_a_ccf
        if merged_glfo is None:
            merged_glfo = glfo
        elif glfo is not None:  # fall through if <glfo> is None
            merged_glfo, name_mapping = glutils.get_merged_glfo(glfo, merged_glfo)
            # Update gene names in annotations to match the merged glfo (only if there are names to remap)
            if any(len(name_mapping[r]) > 0 for r in regions):
                for line in annotation_list:
                    update_gene_names_in_line(line, name_mapping)
        merged_annotation_list += annotation_list
        merged_keys |= set(':'.join(l['unique_ids']) for l in annotation_list)
        if cleanup:
            os.remove(infname)
            os.rmdir(os.path.dirname(infname))

    if getsuffix(outfname) != '.yaml':
        raise Exception('wrong function for %s' % outfname)
    outdir = '.' if os.path.dirname(outfname) == '' else os.path.dirname(outfname)
    mkdir(outdir)

    write_annotations(outfname, merged_glfo, merged_annotation_list, headers, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info, partition_lines=merged_cpath.get_partition_lines())

    if debug:
        print('      read %d total seqs from %d yaml files' % (sum(n_seq_list), len(yaml_list)))

    if return_merged_objects:
        return merged_glfo, merged_annotation_list, merged_cpath
    else:
        return n_event_list, n_seq_list

# ----------------------------------------------------------------------------------------
# merge parameter dirs corresponding to <n_subsets> subsets in <basedir> with str <substr>-<isub> (only works with paired dir structure)
# some things are handled nicelycorrectly, others more hackily
def merge_parameter_dirs(merged_odir, subdfn, n_subsets, include_hmm_cache_files=False, ig_or_tr='ig'):
    from . import glutils, paircluster
    # ----------------------------------------------------------------------------------------
    print('    merging parameters from %d subdirs (e.g. %s) to %s' % (n_subsets, subdfn(0), merged_odir))
    for ltmp in sub_loci(ig_or_tr):
        if os.path.exists('%s/parameters/%s/hmm/germline-sets' % (merged_odir, ltmp)):  # just looks for one of the last thing we would've written
            print('       %s: subset-merged input exists, not rewriting' % locstr(ltmp))
            continue
        mkdir('%s/parameters/%s' % (merged_odir, ltmp))
        def swfn(dname): return '%s/parameters/%s/sw-cache.yaml'%(dname, ltmp)
        sub_swfs = [swfn(subdfn(i)) for i in range(n_subsets) if os.path.exists(swfn(subdfn(i)))]
        if len(sub_swfs) == 0:
            print('       %s: no sw cache files, skipping' % locstr(ltmp))
            continue
        merge_yamls(swfn(merged_odir), sub_swfs, sw_cache_headers, remove_duplicates=True)
        mean_mut_fns = ['%s/parameters/%s/hmm/all-mean-mute-freqs.csv'%(subdfn(i), ltmp) for i in range(n_subsets)]
        mean_mut_fns = [f for f in mean_mut_fns if os.path.exists(f)]
        makelink('%s/parameters/%s/hmm' % (fpath(merged_odir), ltmp), fpath(mean_mut_fns[0]), 'all-mean-mute-freqs.csv')  # NOTE just links to one subset's mut distribution, which should be fine
        merged_glfo, merged_gene_counts = None, {r : defaultdict(int) for r in regions}
        def gpfn(dname, l, r): return '%s/parameters/%s/hmm/%s_gene-probs.csv' % (dname, l, r)
        for isub in range(n_subsets):
            for hfn in glob.glob('%s/parameters/%s/hmm/hmms/*.yaml' % (subdfn(isub), ltmp)):  # these will get overwritten if they're in multiple dirs, which should be fine
                makelink('%s/parameters/%s/hmm/hmms' % (fpath(merged_odir), ltmp), fpath(hfn), os.path.basename(hfn))
            sub_glfo = glutils.read_glfo('%s/parameters/%s/hmm/germline-sets' % (subdfn(isub), ltmp), ltmp, dont_crash=True)
            name_mapping = None
            if merged_glfo is None:
                merged_glfo = sub_glfo
            elif sub_glfo is not None:  # fall through if <sub_glfo> is None
                merged_glfo, name_mapping = glutils.get_merged_glfo(merged_glfo, sub_glfo)
            for treg in regions:
                if not os.path.exists(gpfn(subdfn(isub), ltmp, treg)):
                    continue
                for tline in csvlines(gpfn(subdfn(isub), ltmp, treg)):
                    gene_name = tline['%s_gene'%treg]
                    # Update gene name if it was remapped during glfo merge
                    if name_mapping is not None and gene_name in name_mapping[treg]:
                        gene_name = name_mapping[treg][gene_name]
                    merged_gene_counts[treg][gene_name] += int(tline['count'])
        if merged_glfo is None:  # none of them exists
            continue
        glutils.write_glfo('%s/parameters/%s/hmm/germline-sets' % (merged_odir, ltmp), merged_glfo)
        for treg in regions:
            with open(gpfn(merged_odir, ltmp, treg), 'w') as gfile:
                writer = csv.DictWriter(gfile, ['%s_gene'%treg, 'count'])
                writer.writeheader()
                for gene, count in merged_gene_counts[treg].items():
                    writer.writerow({'%s_gene'%treg : gene, 'count' : count})
        if include_hmm_cache_files:  # these aren't parameters, but don't want to change the name, either, oh well
            subfns = ['%s/single-chain/persistent-cache-%s.csv'%(subdfn(i), ltmp) for i in range(n_subsets)]
            merge_csvs('%s/single-chain/persistent-cache-%s.csv'% (merged_odir, ltmp), subfns)

# ----------------------------------------------------------------------------------------
def get_nodelist_from_slurm_shorthand(nodestr, known_nodes, debug=False):
    if debug:
        print('    getting nodelist from \'%s\'' % nodestr)

    if '[' not in nodestr and ']' not in nodestr:  # single node (don't really need this, but maybe it's a little faster)
        return [nodestr]

    # first find the indices at which there're square braces
    nodes = []
    bracketfo = []
    ilastcomma = -1  # if the first one has brackets, the "effective" comma is at -1
    thisnodestr = ''
    for ich in range(len(nodestr)):
        ch = nodestr[ich]
        if ch == ',':
            ilastcomma = ich
        if ch == '[':
            bracketfo.append({'comma' : ilastcomma, 'ibrackets' : [ich, None]})
            thisnodestr = ''
            if debug:
                print('      start bracket    %s' % nodestr[ilastcomma + 1 : ich + 1])
        elif ch == ']':
            assert bracketfo[-1]['ibrackets'][1] is None
            bracketfo[-1]['ibrackets'][1] = ich
            thisnodestr = ''
            if debug:
                print('      end bracket      %s' % nodestr[bracketfo[-1]['ibrackets'][0] : bracketfo[-1]['ibrackets'][1] + 1])

        # if we're not within a bracket info
        if len(bracketfo) == 0 or bracketfo[-1]['ibrackets'][1] is not None:
            thisnodestr += ch
            thisnodestr = thisnodestr.strip(',[]')
            if len(thisnodestr) > 1 and ch == ',':  # if we just got to a comma, and there's something worth looking at in <thisnodestr>
                nodes.append(thisnodestr)
                thisnodestr = ''
                if debug:
                    print('      add no-bracket   %s' % nodes[-1])

    if debug:
        if len(nodes) > 0:
            print('    %d bracketless nodes: %s' % (len(nodes), ' '.join(nodes)))
        if len(bracketfo) > 0:
            print('      brackets:')

    # the expand the ranges in the brackets
    for bfo in bracketfo:
        ibp = bfo['ibrackets']
        original_str = nodestr[ibp[0] + 1 : ibp[1]]  # NOTE excludes the actual bracket characters
        bracketnodes = []
        for subnodestr in original_str.split(','):
            if '-' in subnodestr:  # range of node numbers
                startstoplist = [int(i) for i in subnodestr.split('-')]
                if len(startstoplist) != 2:
                    raise Exception('wtf %s' % subnodestr)
                istart, istop = startstoplist
                bracketnodes += list(range(istart, istop + 1))
            else: # single node
                bracketnodes.append(int(subnodestr))
        namestr = nodestr[bfo['comma'] + 1 : ibp[0]]  # the texty bit of the name (i.e. without the numbers)
        if debug:
            print('        %s: \'%s\' --> %s' % (namestr, original_str, ' '.join([str(n) for n in bracketnodes])))
        bracketnodes = [namestr + str(i) for i in bracketnodes]
        nodes += bracketnodes

    unknown_nodes = set(nodes) - set(known_nodes)
    if len(unknown_nodes) > 0:
        print('    %s unknown nodes parsed from \'%s\': %s' % (color('yellow', 'warning'), nodestr, ' '.join(unknown_nodes)))

    if debug:
        print('    %d final nodes: %s' % (len(nodes), ' '.join(nodes)))

    return nodes

# ----------------------------------------------------------------------------------------
def get_available_node_core_list(batch_config_fname, debug=False):  # for when you're running the whole thing within one slurm allocation, i.e. with  % salloc --nodes N ./bin/partis [...]
    if debug:
        print('')
        print('  figuring out slurm config')

    our_nodes = []

    if os.getenv('SLURM_NODELIST') is None:  # not within a slurm allocation
        if debug:
            print('  not inside a slurm allocation')
        return None

    # first get info on all nodes from config file
    nodefo = {}  # node : (that node's specifications in config file)
    with open(batch_config_fname) as bfile:
        for line in bfile:
            linefo = line.strip().split()
            node = None
            if len(linefo) > 0 and linefo[0].find('NodeName=') == 0:  # node config line
                for strfo in linefo:
                    tokenval = strfo.split('=')
                    if len(tokenval) != 2:
                        raise Exception('couldn\'t parse %s into \'=\'-separated key-val pairs' % strfo)
                    key, val = tokenval
                    if ',' in val:
                        val = val.split(',')
                    if key == 'NodeName':
                        node = val
                        # if node not in our_nodes:  # damn, doesn't work
                        #     continue
                        nodefo[node] = {}
                    if node is None or node not in nodefo:
                        raise Exception('first key wasn\'t NodeName')
                    nodefo[node][key] = val
    # multiply sockets times cores/socket
    for node, info in nodefo.items():
        if 'Sockets' not in info or 'CoresPerSocket' not in info:
            raise Exception('missing keys in: %s' % ' '.join(list(info.keys())))
        info['nproc'] = int(info['Sockets']) * int(info['CoresPerSocket'])
    if debug:
        print('    info for %d nodes in %s' % (len(nodefo), batch_config_fname))

    our_nodes = get_nodelist_from_slurm_shorthand(os.getenv('SLURM_NODELIST'), known_nodes=list(nodefo.keys()))
    if len(our_nodes) == 0:
        return []

    if debug:
        print('    current allocation includes %d nodes' % len(our_nodes))

    # then info on all current allocations
    quefo = {}  # node : (number of tasks allocated to that node, including ours)
    squeue_str = subprocess.check_output(['squeue', '--format', '%.18i %.2t %.6D %R'], universal_newlines=True)
    headers = ['JOBID', 'ST',  'NODES', 'NODELIST(REASON)']
    for line in squeue_str.split('\n'):
        linefo = line.strip().split()
        if len(linefo) == 0:
            continue
        if linefo[0] == 'JOBID':
            assert linefo == headers
        if linefo[headers.index('ST')] != 'R':  # skip jobs that aren't running
            continue
        nodes = get_nodelist_from_slurm_shorthand(linefo[headers.index('NODELIST(REASON)')], known_nodes=list(nodefo.keys()))
        for node in nodes:
            if node not in our_nodes:
                continue
            if node not in nodefo:
                print('  %s node %s in squeue output but not in config file %s' % (color('yellow', 'warning'), node, batch_config_fname))
                continue
            if node not in quefo:
                quefo[node] = 0
            quefo[node] += 1  # NOTE ideally this would be the number of cores slurm gave this task, rather than 1, but I can't figure out how to that info (and the docs make it sound like it might not be possible to)

    if debug:
        print('    %d "total tasks" allocated among nodes in our current allocation' % sum(quefo.values()))

    # and finally, decide how many procs we can send to each node
    corelist = []
    for node in our_nodes:
        if node not in nodefo:
            raise Exception('node %s in our allocation not in config file %s' % (node, batch_config_fname))
        if node not in quefo:
            raise Exception('node %s in our allocation not in squeue output' % node)
        n_cores_we_can_use = nodefo[node]['nproc'] - quefo[node] + 1  # add one to account for the fact that quefo[node] includes our allocation
        if n_cores_we_can_use == 0:
            print('  %s huh, n_cores_we_can_use is zero' % color('yellow', 'warning'))
            n_cores_we_can_use = 1
        elif n_cores_we_can_use < 0:
            print('  %s more tasks allocated to %s than available cores: %d - %d = %d (setting n_cores_we_can_use to 1 because, uh, not sure what else to do)' % (color('yellow', 'warning'), node, nodefo[node]['nproc'], quefo[node], nodefo[node]['nproc'] - quefo[node]))
            n_cores_we_can_use = 1
        corelist += [node for _ in range(n_cores_we_can_use)]

    corelist = sorted(corelist)  # easier to read if it's alphabetical

    if debug:
        print('    %d available cores:' % len(corelist))
        for node in set(corelist):
            print('        %d  %s' % (corelist.count(node), node))

    if len(corelist) == 0:
        return None

    return corelist

# ----------------------------------------------------------------------------------------
def set_slurm_nodelist(cmdfos, batch_config_fname=None, debug=False):
    # get info about any existing slurm allocation
    corelist = None
    n_procs = len(cmdfos)
    if not os.path.exists(batch_config_fname):
        print('  %s specified --batch-config-fname %s doesn\'t exist' % (color('yellow', 'warning'), batch_config_fname))
    else:
        corelist = get_available_node_core_list(batch_config_fname)  # list of nodes within our current allocation (empty if there isn't one), with each node present once for each core that we've been allocated on that node
        if corelist is not None and len(corelist) < n_procs:
            if 1.5 * len(corelist) < n_procs:
                print('  %s many fewer cores %d than processes %d' % (color('yellow', 'warning'), len(corelist), n_procs))
                print('      corelist: %s' % ' '.join(corelist))
            while len(corelist) < n_procs:
                corelist += sorted(set(corelist))  # add each node once each time through

    if corelist is None:
        return

    if debug:
        print('    %d final cores for %d procs' % (len(corelist), n_procs))
        print('        iproc     node')
        for iproc in range(n_procs):
            print('          %-3d    %s' % (iproc, corelist[iproc]))
    assert len(corelist) >= n_procs
    for iproc in range(n_procs):  # it's kind of weird to keep batch_system and batch_options as keyword args while putting nodelist into the cmdfos, but they're just different enough that it makes sense (we're only using nodelist if we're inside an existing slurm allocation)
        cmdfos[iproc]['nodelist'] = [corelist[iproc]]  # the downside to setting each proc's node list here is that each proc is stuck on that node for each restart (well, unless we decide to change it when we restart it)

# ----------------------------------------------------------------------------------------
def check_cmd(cmd, options='', return_bool=False):  # check for existence of <cmd> (this exists just because check_call() throws a 'no such file or directory' error, and people never figure out that that means the command isn't found)
    try:
        subprocess.check_call([cmd] + options, stdout=open('/dev/null'))
        if return_bool:  # return True if cmd exists + succeeds
            return True
    except OSError:
        if return_bool:
            return False
        else:
            raise Exception('command \'%s\' not found in path (maybe not installed?)' % cmd)

# ----------------------------------------------------------------------------------------
def run_r(cmdlines, workdir, dryrun=False, print_time=None, extra_str='', logfname=None, return_out_err=False, remove_cmdfile=False, debug=False):  # <print_time> is a string which, if set, is printed along with/labeling the time
    if dryrun:
        debug = True
    remove_workdir = False
    if workdir == 'auto':
        workdir = choose_random_subdir('/tmp/%s' % os.getenv('USER'))
        os.makedirs(workdir)
        remove_workdir = True  # NOTE can't use workdir=='auto' since we reset it
    if not os.path.exists(workdir):
        raise Exception('workdir %s doesn\'t exist' % workdir)
    check_cmd('R', options=['--slave', '--version'])
    cmdfname = workdir + '/run.r'
    if debug:
        print('      r cmd lines:')
        print(pad_lines('\n'.join(cmdlines)))
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write('\n'.join(cmdlines) + '\n')
    retval = simplerun('R --slave -f %s' % cmdfname, return_out_err=return_out_err, logfname=logfname, print_time=print_time, extra_str=extra_str, dryrun=dryrun, debug=debug)
    if remove_cmdfile or remove_workdir:
        os.remove(cmdfname)  # different sort of <cmdfname> to that in simplerun()
    if remove_workdir:
        os.rmdir(workdir)
    if return_out_err:
        outstr, errstr = retval
        return outstr, errstr

# ----------------------------------------------------------------------------------------
def mamba_cmds(env, only_prep=False):
    cmds = ['export PYTHONNOUSERSITE=1']  # make sure it doesn't use packages from local environment
    cmds += ['eval "$(micromamba shell hook --shell bash)"']
    if only_prep:
        return cmds
    cmds += ['micromamba activate %s'%env]
    return cmds

# ----------------------------------------------------------------------------------------
def run_ete_script(sub_cmd, return_for_cmdfos=False, dryrun=False, extra_str='', debug=True):
    prof_cmds = '' #' -m cProfile -s tottime -o prof.out'
    # Use QT_QPA_PLATFORM=offscreen for PyQt5 headless rendering (no need for xvfb-run)
    cmd = 'QT_QPA_PLATFORM=offscreen python3%s %s' % (prof_cmds, sub_cmd)
    if debug or dryrun:
        print('%s%s %s' % (extra_str, color('red', 'run'), cmd))
    if return_for_cmdfos:
        return cmd
    else:
        simplerun(cmd, shell=True, dryrun=dryrun, debug=False)

# ----------------------------------------------------------------------------------------
def write_cmd_file(cmd_str, cmdfname, dryrun=False, debug=False):
    mkdir(cmdfname, isfile=True)
    with open(cmdfname, 'w') as cmdfile:
        cmdfile.write(cmd_str)
    subprocess.check_call(['chmod', '+x', cmdfname])
    if debug or dryrun:
        print('      cmd lines in %s:' % cmdfname)
        print(pad_lines('\n'.join(cmd_str.split('\n'))))
        sys.stdout.flush()

# ----------------------------------------------------------------------------------------
def simplerun(cmd_str, shell=False, cmdfname=None, dryrun=False, return_out_err=False, print_time=None, extra_str='', logfname=None, debug=True):
    # ----------------------------------------------------------------------------------------
    def run_subp(cmd_str, shell, fout=None, ferr=None):
        # try:  # this doesn't seeem to do anything as is, but there were some cases where it might help so i'm leaving it here, but commented
        subprocess.check_call(cmd_str if shell else cmd_str.split(), env=os.environ, shell=shell, stdout=fout, stderr=ferr)  # maybe should add executable='/bin/bash'?
        # except subprocess.CalledProcessError as err:
        #     raise Exception() #err)
    # ----------------------------------------------------------------------------------------
    if cmdfname is not None:
        write_cmd_file(cmd_str, cmdfname, dryrun=dryrun, debug=debug)
        cmd_str = cmdfname

    if debug or dryrun:
        print('%s%s %s' % (extra_str, color('red', 'run'), cmd_str))
        if logfname is not None:
            print('        %slog: %s' % (extra_str, logfname))
        sys.stdout.flush()
    if dryrun:
        return '', '' if return_out_err else None
    if print_time is not None:
        start = time.time()

    if return_out_err:
        with tempfile.TemporaryFile() as fout, tempfile.TemporaryFile() as ferr:
            run_subp(cmd_str, shell, fout=fout, ferr=ferr)  # maybe should add executable='/bin/bash'?
            fout.seek(0)
            ferr.seek(0)
            outstr = ''.join(l.decode() for l in fout.readlines())
            errstr = ''.join(l.decode() for l in ferr.readlines())
    else:
        if logfname is not None:  # write cmd_str to logfname, then redirect stdout to it as well
            mkdir(logfname, isfile=True)
            subprocess.check_call('echo %s >%s'%(cmd_str, logfname), shell=True)
            cmd_str = '%s >>%s' % (cmd_str, logfname)
            shell = True
        run_subp(cmd_str, shell)

    if cmdfname is not None:
        os.remove(cmdfname)
    if print_time is not None:
        print('      %s time: %.1f' % (print_time, time.time() - start))
        sys.stdout.flush()

    if return_out_err:
        return outstr, errstr

# ----------------------------------------------------------------------------------------
# return fraction of total system memory that this process is using (as always with memory things, this is an approximation)
# NOTE this is exactly the same number as the stuff in ham/src/bcrutils.cc, so *don't* copy that here and convert to python (again)
def memory_usage_fraction(extra_str='', debug=False):
    if platform.system() != 'Linux':
        print('\n  note: utils.memory_usage_fraction() needs testing on platform \'%s\' to make sure unit conversions don\'t need changing' % platform.system())
    current_usage = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)  # kb
    total = float(psutil.virtual_memory().total) / 1000.  # returns bytes, then convert to kb
    if debug:
        print('  %susing %.0f / %.0f MB = %.4f' % (extra_str, current_usage / 1000, total / 1000, current_usage / total))
    return current_usage / total

# ----------------------------------------------------------------------------------------
def auto_n_procs():  # for running on the local machine
    n_procs = multiprocessing.cpu_count()
    # UPDATE eh whatever just use 'em all
    # if n_procs > 10: # if it's a huge server, we probably shouldn't use all the cores
    #     n_procs = int(float(n_procs) / 2.)
    return n_procs

# ----------------------------------------------------------------------------------------
def limit_procs(base_cmd_str, n_max_procs=None, sleep_time=1, procs=None, debug=False):  # <sleep_time> is seconds
    if base_cmd_str is None:
        if procs is None:
            raise Exception('<base_cmd_str> should be a (string) fragment of the command that will show up in ps, e.g. \'bin/partis\'')
        else:
            def n_running_jobs():
                return [p.poll() for p in procs].count(None)
    else:
        def n_running_jobs():
            tcmd = 'ps auxw | grep %s | grep -v grep | grep -v defunct | grep -v emacs | wc -l' % base_cmd_str  # ignoring emacs lines is hackey, but not sure what else to do
            return int(subprocess.check_output(tcmd, shell=True, universal_newlines=True))
    if n_max_procs is None:
        n_max_procs = auto_n_procs()
    n_jobs = n_running_jobs()
    while n_jobs >= n_max_procs:
        if debug:
            print('%d (>=%d) running jobs' % (n_jobs, n_max_procs))
        time.sleep(sleep_time)
        n_jobs = n_running_jobs()

# ----------------------------------------------------------------------------------------
def run_proc_functions(procs, n_procs=None, debug=False):  # <procs> is a list of multiprocessing.Process objects
    if n_procs is None:
        n_procs = auto_n_procs()
    if debug:
        print('    running %d proc fcns with %d procs' % (len(procs), n_procs))
        sys.stdout.flush()
    while True:
        while len(procs) > 0 and len(multiprocessing.active_children()) < n_procs:
            procs[0].start()
            procs.pop(0)
        if len(multiprocessing.active_children()) == 0 and len(procs) == 0:
            break

# ----------------------------------------------------------------------------------------
def get_batch_system_str(batch_system, cmdfo, fout, ferr, batch_options):
    prestr = ''

    if batch_system == 'slurm':
        prestr += 'srun --export=ALL --nodes 1 --ntasks 1'  # --exclude=data/gizmod.txt'  # --export=ALL seems to be necessary so XDG_RUNTIME_DIR gets passed to the nodes
        if 'threads' in cmdfo:
            prestr += ' --cpus-per-task %d' % cmdfo['threads']
        if 'nodelist' in cmdfo:
            prestr += ' --nodelist ' + ','.join(cmdfo['nodelist'])
    elif batch_system == 'sge':
        prestr += 'qsub -sync y -b y -V -o ' + fout + ' -e ' + ferr
        fout = None
        ferr = None
    else:
        assert False

    if batch_options is not None:
        prestr += ' ' + batch_options

    return prestr, fout, ferr

# ----------------------------------------------------------------------------------------
def cycle_log_files(logfname, debug=False):  # move any existing log file to .0, .1, etc.
    if not os.path.exists(logfname):  # nothing to do
        return
    itmp = 0
    while os.path.exists('%s.%d'%(logfname, itmp)):
        itmp += 1
    if debug:
        print('  %s --> %s' % (logfname, '%s.%d'%(logfname, itmp)))
    os.rename(logfname, '%s.%d'%(logfname, itmp))

# ----------------------------------------------------------------------------------------
def run_cmd(cmdfo, batch_system=None, batch_options=None, shell=False):
    cstr = cmdfo['cmd_str']  # don't want to modify the str in <cmdfo>
    fout = cmdfo['logdir'] + '/out'
    ferr = cmdfo['logdir'] + '/err'
    cycle_log_files(fout)
    cycle_log_files(ferr)

    if batch_system is not None:
        prestr, fout, ferr = get_batch_system_str(batch_system, cmdfo, fout, ferr, batch_options)
        cstr = prestr + ' ' + cstr

    mkdir(cmdfo['logdir'])

    proc = subprocess.Popen(cstr if shell else cstr.split(),
                            stdout=None if fout is None else open(fout, 'w'),
                            stderr=None if ferr is None else open(ferr, 'w'),
                            env=cmdfo.get('env'), shell=shell, universal_newlines=True,
                            executable='/bin/bash' if shell else None)  # adding executable= very late, not sure if it'll break something somewhere
    return proc

# ----------------------------------------------------------------------------------------
# <cmdfos> list of dicts, each dict specifies how to run one process, entries:
cmdfo_required_keys = [
    'cmd_str',  #  actual command to run
    'outfname',  # output file resulting from 'cmd_str'. Used to determine if command completed successfully
    'workdir',  # either this or 'logdir' must be set. If <clean_on_success> is set, this directory is removed when finished. Also serves as default for 'logdir'. (ok this is kind of messy, this probably used to be used for more things)
]
cmdfo_defaults = {  # None means by default it's absent
    'logdir' : 'workdir',  # where to write stdout/stderr (written as the files 'out' and 'err'). If not set, they go in 'workdir' (I don't like this default, but there's way too much code that depends on it to change it now)
    'workfnames' : None,  # if <clean_on_success> is set, this list of files is deleted before removing 'workdir' upon successful completion
    'dbgfo' : None,  # dict to store info from bcrham stdout about how many viterbi/forward/etc. calculations it performed
    'env' : None,  # if set, passed to the env= keyword arg in Popen
    'nodelist' : None,  # list of slurm nodes to allow; do not use, only set automatically
    'threads' : None,  # slurm cpus per task
}

# ----------------------------------------------------------------------------------------
# notes:
#  - set sleep to False if your commands are going to run really really really quickly
#  - unlike everywhere else, <debug> is not a boolean, and is either None (swallow out, print err)), 'print' (print out and err), 'write' (write out and err to file called 'log' in logdir), or 'write:<log file name>' (same as 'write', but you set your own base name)
#  - if both <n_max_procs> and <proc_limit_str> are set, it uses limit_procs() (i.e. a ps call) to count the total number of <proc_limit_str> running on the machine; whereas if only <n_max_procs> is set, it counts only subprocesses that it is itself running
#  - debug: can be None (stdout mostly gets ignored), 'print' (printed), 'write' (written to file 'log' in logdir), or 'write:<logfname>' (same, but use <logfname>)
def run_cmds(cmdfos, shell=False, n_max_tries=None, clean_on_success=False, batch_system=None, batch_options=None, batch_config_fname=None,
             debug=None, ignore_stderr=False, sleep=True, n_max_procs=None, proc_limit_str=None, allow_failure=False):
    if len(cmdfos) == 0:
        raise Exception('zero length cmdfos')
    if n_max_tries is None:
        n_max_tries = 1 if batch_system is None else 3
    per_proc_sleep_time = 0.01 / max(1, len(cmdfos))

    # check cmdfos and set defaults
    if len(set(cmdfo_required_keys) - set(cmdfos[0])) > 0:
        raise Exception('missing required keys in cmdfos: %s' % ' '.join(set(cmdfo_required_keys) - set(cmdfos[0])))
    if len(set(cmdfos[0]) - set(cmdfo_required_keys) - set(cmdfo_defaults)) > 0:
        raise Exception('unexpected key in cmdfos: %s' % ' '.join(set(cmdfos[0]) - set(cmdfo_required_keys) - set(cmdfo_defaults)))
    for iproc in range(len(cmdfos)):  # ok this is way overcomplicated now that I'm no longer adding None ones by default, oh well
        for ckey, dval in cmdfo_defaults.items():
            if ckey not in cmdfos[iproc] and dval is not None:
                cmdfos[iproc][ckey] = cmdfos[iproc][dval] if dval in cmdfos[iproc] else None  # first bit is only used for using workdir as the logdir default a.t.m.

    if batch_system == 'slurm' and batch_config_fname is not None:
        set_slurm_nodelist(cmdfos, batch_config_fname)

    procs, n_tries_list = [], []
    for iproc in range(len(cmdfos)):
        procs += [run_cmd(cmdfos[iproc], batch_system=batch_system, batch_options=batch_options, shell=shell)]
        n_tries_list.append(1)
        if sleep:
            time.sleep(per_proc_sleep_time)
        if n_max_procs is not None:
            limit_procs(proc_limit_str, n_max_procs, procs=procs)  # NOTE now that I've added the <procs> arg, I should remove all the places where I'm using the old cmd str method (I mean, it works fine, but it's hackier/laggier, and in cases where several different parent procs are running a log of the same-named subprocs on the same machine, the old way will be wrong [i.e. limit_procs was originally intended as a global machine-wide limit, whereas in this fcn we usually call it wanting to set a specific number of subproces for this process])

    dbgstrs = ['' for _ in procs]
    while procs.count(None) != len(procs):  # we set each proc to None when it finishes
        for iproc in range(len(cmdfos)):
            if procs[iproc] is None:  # already finished
                continue
            if procs[iproc].poll() is not None:  # it just finished
                status, dbgstrs[iproc] = finish_process(iproc, procs, n_tries_list[iproc], cmdfos[iproc], n_max_tries, dbgfo=cmdfos[iproc].get('dbgfo'), batch_system=batch_system, debug=debug, ignore_stderr=ignore_stderr, clean_on_success=clean_on_success, allow_failure=allow_failure)
                if status == 'restart':
                    print(dbgstrs[iproc])
                    procs[iproc] = run_cmd(cmdfos[iproc], batch_system=batch_system, batch_options=batch_options, shell=shell)
                    n_tries_list[iproc] += 1
        sys.stdout.flush()
        if sleep:
            time.sleep(per_proc_sleep_time)
    for dstr in dbgstrs:
        if dstr != '':
            print(dstr)

# ----------------------------------------------------------------------------------------
def pad_lines(linestr, padwidth=8):
    lines = [padwidth * ' ' + l for l in linestr.split('\n')]
    return '\n'.join(lines)

# ----------------------------------------------------------------------------------------
def get_slurm_node(errfname):
    if not os.path.exists(errfname):
        return None

    jobid = None
    try:
        jobid = subprocess.check_output(['head', '-n1', errfname], universal_newlines=True).split()[2]
    except (subprocess.CalledProcessError, IndexError) as err:
        print(err)
        print('      couldn\'t get jobid from err file %s with contents:' % errfname)
        subprocess.check_call(['cat', errfname])
        return None

    assert jobid is not None
    nodelist = None
    try:
        nodelist = subprocess.check_output(['squeue', '--job', jobid, '--states=all', '--format', '%N'], universal_newlines=True).split()[1]
    except (subprocess.CalledProcessError, IndexError) as err:
        print(err)
        print('      couldn\'t get node list from jobid \'%s\'' % jobid)
        return None

    assert nodelist is not None
    if ',' in nodelist:  # I think this is what it looks like if there's more than one, but I'm not checking
        raise Exception('multiple nodes in nodelist: \'%s\'' % nodelist)

    return nodelist

# ----------------------------------------------------------------------------------------
# deal with a process once it's finished (i.e. check if it failed, and tell the calling fcn to restart it if so)
def finish_process(iproc, procs, n_tried, cmdfo, n_max_tries, dbgfo=None, batch_system=None, debug=None, ignore_stderr=False, clean_on_success=False, allow_failure=False):
    # ----------------------------------------------------------------------------------------
    def logfname(ltype):
        return cmdfo['logdir'] + '/' + ltype
    # ----------------------------------------------------------------------------------------
    def getlogstrs(logtypes=None):  # not actually using stdout at all, but maybe I should?
        if logtypes is None:
            logtypes = ['out', 'err']
        returnstr = []
        for ltype in logtypes:
            if os.path.exists(logfname(ltype)) and os.stat(logfname(ltype)).st_size > 0:
                returnstr += ['        %s           %s' % (color('red', 'std%s:'%ltype), logfname(ltype))]
                returnstr += [pad_lines(subprocess.check_output(['cat', logfname(ltype)], universal_newlines=True), padwidth=12)]
        return '\n'.join(returnstr)
    # ----------------------------------------------------------------------------------------
    procs[iproc].communicate()
    outfname = cmdfo['outfname']
    rtn_strs = []

    # success
    if procs[iproc].returncode == 0:
        if not os.path.exists(outfname):
            rtn_strs.append('      proc %d succeeded but its output isn\'t there, so sleeping for a bit: %s' % (iproc, outfname))  # give a networked file system some time to catch up
            time.sleep(0.5)
        if os.path.exists(outfname):
            oestr = process_out_err(cmdfo['logdir'], extra_str='' if len(procs) == 1 else str(iproc), dbgfo=dbgfo, cmd_str=cmdfo['cmd_str'], debug=debug, ignore_stderr=ignore_stderr)
            if debug == 'print':
                rtn_strs.append(oestr)
            procs[iproc] = None  # job succeeded
            if clean_on_success:  # this is newer than the rest of the fcn, so it's only actually used in one place, but it'd be nice if other places started using it eventually
                if cmdfo.get('workfnames') is not None:
                    for fn in [f for f in cmdfo['workfnames'] if os.path.exists(f)]:
                        os.remove(fn)
                if os.path.isdir(cmdfo['workdir']):
                    os.rmdir(cmdfo['workdir'])
            return 'ok', '\n'.join(rtn_strs)

    # handle failure
    status_str = '    proc %d try %d' % (iproc, n_tried)
    if procs[iproc].returncode == 0 and not os.path.exists(outfname):  # don't really need both the clauses
        status_str += ' succeeded but output is missing: %s' % outfname
    else:
        status_str += ' failed with exit code %d, %s: %s' % (procs[iproc].returncode, 'but output exists' if os.path.exists(outfname) else 'and output is missing',  outfname)
    if batch_system == 'slurm':  # cmdfo['cmd_str'].split()[0] == 'srun' and
        if 'nodelist' in cmdfo:  # if we're doing everything from within an existing slurm allocation
            nodelist = cmdfo['nodelist']
        else:  # if, on the other hand, each process made its own allocation on the fly
            nodelist = get_slurm_node(cmdfo['logdir'] + '/err')
        if nodelist is not None:
            status_str += '     failed on node %s' % nodelist
        # try:
        #     print '        sshing to %s' % nodelist
        #     outstr = check_output('ssh -o StrictHostKeyChecking=no ' + nodelist + ' ps -eo pcpu,pmem,rss,cputime:12,stime:7,user,args:100 --sort pmem | tail', shell=True, universal_newlines=True)
        #     print pad_lines(outstr, padwidth=12)
        # except subprocess.CalledProcessError as err:
        #     print '        failed to ssh:'
        #     print err
    rtn_strs.append(status_str)
    if os.path.exists(outfname + '.progress'):  # glomerator.cc is the only one that uses this at the moment
        rtn_strs.append('        progress file (%s):' % (outfname + '.progress'))
        rtn_strs.append(pad_lines(subprocess.check_output(['cat', outfname + '.progress'], universal_newlines=True), padwidth=12))
    if n_tried < n_max_tries:
        rtn_strs.append(getlogstrs(['err']))
        rtn_strs.append('      restarting proc %d' % iproc)
        return 'restart', '\n'.join(rtn_strs)
    else:
        failstr = 'exceeded max number of tries (%d >= %d) for subprocess with command:\n        %s\n' % (n_tried, n_max_tries, cmdfo['cmd_str'])
        failstr += getlogstrs()  # used to try to only print err/out as needed, but it was too easy to miss useful info
        if allow_failure:
            rtn_strs.append('      %s\n      not raising exception for failed process' % failstr)
            procs[iproc] = None  # let it keep running any other processes
        else:
            raise Exception(failstr)

    return 'failed', '\n'.join(rtn_strs)

# ----------------------------------------------------------------------------------------
def process_out_err(logdir, extra_str='', dbgfo=None, cmd_str=None, debug=None, ignore_stderr=False):
    # NOTE something in this chain seems to block or truncate or some such nonsense if you make it too big
    err_strs_to_ignore = [
        'stty: standard input: Inappropriate ioctl for device',
        'queued and waiting for resources',
        'has been allocated resources',
        'srun: Required node not available (down, drained or reserved)',
        'GSL_RNG_TYPE=',
        'GSL_RNG_SEED=',
        '[ig_align] Read',
        '[ig_align] Aligned',
    ]
    def read_and_delete_file(fname):
        fstr = ''
        if os.stat(fname).st_size > 0:  # NOTE if <fname> doesn't exist, it probably means you have more than one process writing to the same log file
            ftmp = open(fname, encoding='ISO-8859-1')  # this encoding seems necessary for bppseqgen output, and doesn't seem to break anything else, at least yet?
            fstr = ''.join(ftmp.readlines())
            ftmp.close()
        os.remove(fname)
        return fstr
    def skip_err_line(line):
        if len(line.strip()) == 0:
            return True
        for tstr in err_strs_to_ignore:
            if tstr in line:
                return True
        return False

    logstrs = {tstr : read_and_delete_file(logdir + '/' + tstr) for tstr in ['out', 'err']}

    err_str = []
    for line in logstrs['err'].split('\n'):
        if skip_err_line(line):
            continue
        err_str += [line]
    err_str = '\n'.join(err_str)

    if 'bcrham' in cmd_str:
        for line in logstrs['out'].split('\n'):  # print debug info related to --n-final-clusters/--min-largest-cluster-size force merging
            if 'force' in line:
                print('    %s %s' % (color('yellow', 'force info:'), line))

    if dbgfo is not None:  # keep track of how many vtb and fwd calculations the process made
        for header, variables in bcrham_dbgstrs['partition'].items():  # 'partition' is all of them, 'annotate' is a subset
            dbgfo[header] = {var : None for var in variables}
            theselines = [ln for ln in logstrs['out'].split('\n') if header + ':' in ln]
            if len(theselines) == 0:
                continue
            if len(theselines) > 1:
                raise Exception('too many lines with dbgfo for \'%s\' in:\nstdout:\n%s\nstderr:\n%s' % (header, logstrs['out'], logstrs['err']))
            words = theselines[0].split()
            for var in variables:  # convention: value corresponding to the string <var> is the word immediately vollowing <var>
                if var in words:
                    if words.count(var) > 1:
                        raise Exception('found multiple instances of variable \'%s\' in line \'%s\'' % (var, theselines[0]))
                    dbgfo[header][var] = float(words[words.index(var) + 1])

    rtn_strs = []
    if debug is None:
        if not ignore_stderr and len(err_str) > 0:
            print(err_str)
    elif len(err_str) + len(logstrs['out']) > 0:
        if debug == 'print':
            if extra_str != '':
                tmpcolor = 'red_bkg' if len(err_str + logstrs['out']) != len_excluding_colors(err_str + logstrs['out']) else None  # if there's color in the out/err strs, make the 'proc 0' str colored as well
                rtn_strs.append('      --> %s' % color(tmpcolor, 'proc %s' % extra_str))
            rtn_strs.append(err_str + logstrs['out'])
        elif 'write' in debug:
            if debug == 'write':
                logfile = logdir + '/log'
            else:
                assert debug[:6] == 'write:'
                logfile = logdir + '/' + debug.replace('write:', '')
            cycle_log_files(logfile)
            with open(logfile, 'w') as dbgfile:
                if cmd_str is not None:
                    dbgfile.write('%s %s\n' % (color('red', 'run'), cmd_str))  # NOTE duplicates code in datascripts/run.py
                dbgfile.write(err_str + logstrs['out'])
        else:
            assert False

        return '\n'.join(rtn_strs)

# ----------------------------------------------------------------------------------------
def summarize_bcrham_dbgstrs(dbgfos, action):
    def defval(dbgcat):
        if dbgcat in bcrham_dbgstr_types[action]['sum']:
            return 0.
        elif dbgcat in bcrham_dbgstr_types[action]['same']:
            return None
        elif dbgcat in bcrham_dbgstr_types[action]['min-max']:
            return []
        else:
            assert False

    cache_read_inconsistency = False
    summaryfo = {dbgcat : {vtype : defval(dbgcat) for vtype in tlist} for dbgcat, tlist in bcrham_dbgstrs[action].items()}  # fill summaryfo with default/initial values
    for procfo in dbgfos:
        for dbgcat in bcrham_dbgstr_types[action]['same']:  # loop over lines in output for which every process should have the same values (e.g. cache-read)
            for vtype in bcrham_dbgstrs[action][dbgcat]:  # loop over values in that line (e.g. logprobs and naive-seqs)
                if summaryfo[dbgcat][vtype] is None:  # first one
                    summaryfo[dbgcat][vtype] = procfo[dbgcat][vtype]
                if procfo[dbgcat][vtype] != summaryfo[dbgcat][vtype]:  # make sure all subsequent ones are the same
                    cache_read_inconsistency = True
                    print('        %s bcrham procs had different \'%s\' \'%s\' info: %d vs %d' % (color('red', 'warning'), vtype, dbgcat, procfo[dbgcat][vtype], summaryfo[dbgcat][vtype]))
        for dbgcat in bcrham_dbgstr_types[action]['sum']:  # lines for which we want to add up the values
            for vtype in bcrham_dbgstrs[action][dbgcat]:
                if procfo[dbgcat][vtype] is None:  # can't seem to replicate this, but it happened once
                    print('  %s none type dbg info read from subprocess (maybe the batch system made the subprocess print out something extra so it didn\'t parse correctly?)' % color('yellow', 'warning'))
                else:
                    summaryfo[dbgcat][vtype] += procfo[dbgcat][vtype]
        for dbgcat in bcrham_dbgstr_types[action]['min-max']:  # lines for which we want to keep track of the smallest and largest values (e.g. time required)
            for vtype in bcrham_dbgstrs[action][dbgcat]:
                summaryfo[dbgcat][vtype].append(procfo[dbgcat][vtype])

    for dbgcat in bcrham_dbgstr_types[action]['min-max']:
        for vtype in bcrham_dbgstrs[action][dbgcat]:
            summaryfo[dbgcat][vtype] = min(summaryfo[dbgcat][vtype]), max(summaryfo[dbgcat][vtype])

    if cache_read_inconsistency:
        raise Exception('inconsistent cache reading information across processes (see above), probably due to file system issues')

    return summaryfo

# ----------------------------------------------------------------------------------------
def find_first_non_ambiguous_base(seq):
    """ return index of first non-ambiguous base """
    for ib in range(len(seq)):
        if seq[ib] not in all_ambiguous_bases:
            return ib
    assert False  # whole sequence was ambiguous... probably shouldn't get here

# ----------------------------------------------------------------------------------------
def find_last_non_ambiguous_base_plus_one(seq):
    for ib in range(len(seq) - 1, -1, -1):  # count backwards from the end
        if seq[ib] not in all_ambiguous_bases:  # find first non-ambiguous base
            return ib + 1  # add one for easy slicing

    assert False  # whole sequence was ambiguous... probably shouldn't get here

# ----------------------------------------------------------------------------------------
def remove_ambiguous_ends(seq):
    """ remove ambiguous bases from the left and right ends of <seq> """
    i_seq_start = find_first_non_ambiguous_base(seq)
    i_seq_end = find_last_non_ambiguous_base_plus_one(seq)
    return seq[i_seq_start : i_seq_end]

# ----------------------------------------------------------------------------------------
def split_clusters_by_cdr3(partition, sw_info, warn=False):  # NOTE maybe should combine this with some of the following fcns?
    new_partition = []
    all_cluster_splits = []
    for cluster in partition:
        cdr3_lengths = [sw_info[q]['cdr3_length'] for q in cluster]
        if len(set(cdr3_lengths)) == 1:
            new_partition.append(cluster)
        else:
            split_clusters = [list(group) for _, group in itertools.groupby(sorted(cluster, key=lambda q: sw_info[q]['cdr3_length']), key=lambda q: sw_info[q]['cdr3_length'])]  # TODO i think this should really be using group_seqs_by_value()
            for sclust in split_clusters:
                new_partition.append(sclust)
            all_cluster_splits.append((len(cluster), [len(c) for c in split_clusters]))
    if warn and len(all_cluster_splits) > 0:
        print('  %s split apart %d cluster%s that contained multiple cdr3 lengths (total clusters: %d --> %d)' % (color('yellow', 'warning'), len(all_cluster_splits), plural(len(all_cluster_splits)), len(partition), len(new_partition)))
        print('      cluster splits: %s' % ', '.join(('%3d --> %s'%(cl, ' '.join(str(l) for l in spls))) for cl, spls in all_cluster_splits))
    return new_partition

# ----------------------------------------------------------------------------------------
def split_partition_with_criterion(partition, criterion_fcn):  # this would probably be faster if I used the itertools stuff from collapse_naive_seqs()
    true_cluster_indices = [ic for ic in range(len(partition)) if criterion_fcn(partition[ic])]  # indices of clusters for which <criterion_fcn()> is true
    true_clusters = [partition[ic] for ic in true_cluster_indices]
    false_clusters = [partition[ic] for ic in range(len(partition)) if ic not in true_cluster_indices]
    return true_clusters, false_clusters

# ----------------------------------------------------------------------------------------
def group_seqs_by_value(queries, keyfunc, return_values=False):  # <queries> don't have to be related to seqs at all, only requirement is that the things in the iterable <queries> have to be valid arguments to <keyfunc()>
    vals, groups = zip(*[(val, sorted(list(group))) for val, group in itertools.groupby(sorted(queries, key=keyfunc), key=keyfunc)])  # NOTE sorted() around list(group) is so [in python 3] we'll get the same list order in reruns
    if return_values:
        return list(zip(*(vals, groups)))
    else:
        return list(groups)

# ----------------------------------------------------------------------------------------
def collapse_seqfos_with_identical_seqs(input_seqfos, keys_not_to_collapse=None, debug=False):
    # ----------------------------------------------------------------------------------------
    def kfcn(u):
        rstr = sfo_dict[u]['seq']
        if keys_not_to_collapse is not None:
            def kstr(v): return str(plotting.legend_groups.get(v, v))
            rstr = '-@-'.join([rstr] + [kstr(sfo_dict[u][k]) for k in keys_not_to_collapse])
        return rstr
    # ----------------------------------------------------------------------------------------
    if keys_not_to_collapse is not None:
        if any(k not in sfo for sfo in input_seqfos for k in keys_not_to_collapse):
            raise Exception('found at least one seqfo without one of the keys_not_to_collapse %s' % keys_not_to_collapse)
        from . import plotting  # UGH
    newfos = []
    sfo_dict = {s['name'] : s for s in input_seqfos}
    for uids in group_seqs_by_value(sorted(sfo_dict.keys()), kfcn):
        tseq = sfo_dict[uids[0]]['seq']
        extra_keys = set(k for u in uids for k in sfo_dict[u] if k not in ['name', 'seq', 'multiplicity'] + non_none([keys_not_to_collapse, []]))
        if len(extra_keys) > 0:  # should really just check if the values are the same in all the seqfos we're merging, and if they are then copy values
            print('    %s extra keys (beyond name, seq, multiplicity) found in seqfos to merge, new seqfos won\'t have them: %s' % (wrnstr(), ' '.join(extra_keys)))
        newfos.append({'name' : uids[0], 'seq' : tseq, 'multiplicity' : len(uids), 'duplicates' : uids[1:]})
    if debug:
        estr = '' if keys_not_to_collapse is None else ' (used also keys: %s)' % ' '.join(keys_not_to_collapse)
        print('    collapsed %d total seqs (with duplicates) into %d unique seqs before running tree inference%s' % (len(input_seqfos), len(newfos), estr))
        # print('         removed seqs: %s' % sorted(set(sfo_dict) - set(s['name'] for s in newfos)))
    return newfos

# ----------------------------------------------------------------------------------------
# takes sw info (or any dict keyed by uid with 'naive_seq' key in each uid's dict) and returns a partition of uids, where each cluster contains all uids with the same nave seq
def collapse_naive_seqs(swfo, queries=None, split_by_cdr3=False, debug=None):  # <split_by_cdr3> is only needed when we're getting synthetic sw info that's a mishmash of hmm and sw annotations
    start = time.time()
    if queries is None:
        queries = swfo['queries']  # don't modify this

    def keyfunc(q):  # while this is no longer happening before fwk insertion trimming (which was bad), it is still happening on N-padded sequences, which should be kept in mind
        if split_by_cdr3:
            return swfo[q]['cdr3_length'], swfo[q]['naive_seq']
        else:
            return swfo[q]['naive_seq']

    partition = group_seqs_by_value(queries, keyfunc)

    if debug:
        print('   collapsed %d queries into %d cluster%s with identical naive seqs (%.1f sec)' % (len(queries), len(partition), plural(len(partition)), time.time() - start))

    return partition

# ----------------------------------------------------------------------------------------
def collapse_naive_seqs_with_hashes(naive_seq_list, sw_info):  # this version is (atm) only used for naive vsearch clustering
    start = time.time()
    naive_seq_map = {}  # X[cdr3][hash(naive_seq)] : naive_seq
    naive_seq_hashes = {}  # X[cdr3][hash(naive_seq)] : [ustr1, ustr2, ustr3...]  # NOTE didn't used to be also subset by [cdr3], but it seems that they can have different cdr3 but same naive seq, which screws up untranslation
    for ustr, naive_seq in naive_seq_list:
        hashstr = uidhashstr(naive_seq)
        c3len = sw_info[ustr]['cdr3_length']
        if c3len not in naive_seq_map:
            naive_seq_map[c3len], naive_seq_hashes[c3len] = {}, {}
        if hashstr not in naive_seq_map[c3len]:
            naive_seq_map[c3len][hashstr] = naive_seq  # i.e. vsearch gets a hash of the naive seq (which maps to a list of ustrs with that naive sequence) instead of the ustr
            naive_seq_hashes[c3len][hashstr] = []  # first sequence that has this naive
        naive_seq_hashes[c3len][hashstr].append(ustr)
    print('        collapsed %d sequences (%d initial naive seqs/clusters) into %d unique naive sequences over %d cdr3 lengths (%.1f sec)' % (sum(len(u.split(':')) for u, _ in naive_seq_list), len(naive_seq_list), sum(len(d) for d in naive_seq_hashes.values()), len(naive_seq_hashes), time.time() - start))
    return naive_seq_map, naive_seq_hashes

# ----------------------------------------------------------------------------------------
def get_partition_from_annotation_list(annotation_list):
    return [copy.deepcopy(l['unique_ids']) for l in annotation_list]

# ----------------------------------------------------------------------------------------
def get_partition_from_reco_info(reco_info, ids=None):
    # Two modes:
    #  - if <ids> is None, it returns the actual, complete, true partition.
    #  - if <ids> is set, it groups them into the clusters dictated by the true partition in/implied by <reco_info> NOTE these are not, in general, complete clusters
    if ids is None:
        ids = list(reco_info.keys())
    def keyfunc(q):
        return reco_info[q]['reco_id']
    return [list(group) for _, group in itertools.groupby(sorted(ids, key=keyfunc), key=keyfunc)]  # sort 'em beforehand so all the ones with the same reco id are consecutive (if there are non-consecutive ones with the same reco id, it means there were independent events with the same rearrangment parameters)

# ----------------------------------------------------------------------------------------
def get_partition_from_str(partition_str):
    """ NOTE there's code in some other places that do the same thing """
    clusters = partition_str.split(';')
    partition = [cl.split(':') for cl in clusters]
    return partition

# ----------------------------------------------------------------------------------------
def get_str_from_partition(partition):
    """ NOTE there's code in some other places that do the same thing """
    clusters = [':'.join(cl) for cl in partition]
    partition_str = ';'.join(clusters)
    return partition_str

# ----------------------------------------------------------------------------------------
def get_clone_id(cluster):  # get unique hash id for <cluster>
    return uidhashstr(':'.join(cluster))

# ----------------------------------------------------------------------------------------
def get_cluster_ids(uids, partition):  # return map from each uid in <uids> to the index of its cluster in <partition>
    clids = {uid : [] for uid in uids}  # almost always list of length one with index (in <partition>) of the uid's cluster
    for iclust in range(len(partition)):
        for uid in partition[iclust]:
            if iclust not in clids[uid]:  # in case there's duplicates (from seed unique id)
                clids[uid].append(iclust)
    return clids

# ----------------------------------------------------------------------------------------
# return a new list of partitions that has no duplicate uids (choice as to which cluster gets to keep a duplicate id is entirely random [well, it's the first one that has it, so not uniform random, but you can't specify it])
def get_deduplicated_partitions(partitions, antn_list=None, glfo=None, debug=False):  # not using this atm since i wrote it for use in clusterpath, but then ended up not needing it UPDATE now using it during paired clustering resolution, but maybe only temporarily
    # ----------------------------------------------------------------------------------------
    def try_to_remove_annotation(cluster, new_cluster):
        ilines = [i for i, l in enumerate(antn_list) if l['unique_ids']==cluster]
        if len(ilines) == 0:
            print('  %s antn_list was set in get_deduplicated_partitions(), but there wasn\'t an annotation for cluster %s' % (color('yellow', 'warning'), cluster))  # to check for overlapping (but not identical) ones: print [l for l in antn_list if len(set(l['unique_ids']) & cluster) > 0]
        else:
            if len(new_cluster) == 0:  # removed all the uids
                antn_list.pop(ilines[0])
            else:
                tline = antn_list[ilines[0]]
                iseqs_to_keep = [i for i, u in enumerate(tline['unique_ids']) if u in new_cluster]
                restrict_to_iseqs(tline, iseqs_to_keep, glfo)
                if debug:
                    print('          successfully removed %d/%d seqs from corresponding annotation' % (len(cluster) - len(new_cluster), len(cluster)))
    # ----------------------------------------------------------------------------------------
    if debug:
        print('    deduplicating %d partition%s' % (len(partitions), plural(len(partitions))))
    new_partitions = [[] for _ in partitions]
    for ipart in range(len(partitions)):
        if debug:
            unique_uids, total_uids = len(set(u for c in partitions[ipart] for u in c)), sum(len(c) for c in partitions[ipart])
            print('      ipart %d with %d clusters: %d unique vs %d total uids (expect to remove %d)' % (ipart, len(partitions[ipart]), unique_uids, total_uids, total_uids - unique_uids))
            duplicated_uids = set()
            dbg_strs = []
        previously_encountered_uids = set()
        n_removed = 0
        for cluster in partitions[ipart]:
            new_uids = sorted(set(cluster) - previously_encountered_uids)  # remove any uids that were in previous clusters
            new_cluster = [u for u in cluster if u in new_uids]  # then make sure the order stays the same
            previously_encountered_uids |= set(new_cluster)
            if len(new_cluster) > 0:
                new_partitions[ipart].append(new_cluster)
            if len(new_cluster) < len(cluster):  # if we're actually removing any uids
                if antn_list is not None:
                    try_to_remove_annotation(cluster, new_cluster)
                if debug:
                    dbg_strs.append('%d/%d' % (len(cluster) - len(new_cluster), len(cluster)))
                    duplicated_uids |= set(cluster) - set(new_cluster)
                    n_removed += len(cluster) - len(new_cluster)
        if debug:
            dbg_strs = sorted(dbg_strs, key=lambda x: int(x.split('/')[1]), reverse=True)  # sort by denominator, i.e. original cluster size
            print('      %d uids appeared more than once (and %d total were removed)%s' % (len(duplicated_uids), n_removed, (':  ' + ' '.join(sorted(duplicated_uids))) if len(duplicated_uids) < 10 else ''))
            n_singletons = dbg_strs.count('1/1')  # NOTE similarity to utils.cluster_size_str()
            print('           removed uids/from clusters with size: %s (+%d singletons)' % ('  '.join(s for s in dbg_strs if s!='1/1'), n_singletons))
    return new_partitions

# ----------------------------------------------------------------------------------------
# return a new list of partitions where any clusters that share uids have been merged (including, if set, fixing antn_list to match)
def merge_overlapping_clusters(partitions, antn_list=None, glfo=None, dbgstr='', debug=False):  # similar to deduplicating fcn above
    # ----------------------------------------------------------------------------------------
    def try_to_update_annotations(cluster_i, cluster_j, newclust):
        if newclust not in [cluster_i, cluster_j]:  # if needed, get the new annotation
            if any(':'.join(c) not in antn_dict for c in [cluster_i, cluster_j]):
                print('      %s missing annotation for cluster[s] when trying to update annotations, so skipping: %s' % (wrnstr(), ', '.join([':'.join(c) for c in [cluster_i, cluster_j] if ':'.join(c) not in antn_dict])))
                return
            antn_i, antn_j = antn_dict[':'.join(cluster_i)], antn_dict[':'.join(cluster_j)]  # if it crashes here, the annotations aren't sync'd with the partition
            antn_dict[':'.join(newclust)] = combine_events(glfo, [antn_i, antn_j], extra_str='                        ', debug=debug)
        for tclust in [cluster_i, cluster_j]:
            if tclust != newclust and ':'.join(tclust) in antn_dict:  # maybe should warn if they're missing, but it's probably ok, there's lots of reasons clusters and annotations aren't synced
                del antn_dict[':'.join(tclust)]
    # ----------------------------------------------------------------------------------------
    def process_partition(newptn):
        n_changes = 0  # note that this the number of changes made, which in general can be much larger than the number of clusters
        for iclust in range(len(newptn)):
            for jclust in range(iclust + 1, len(newptn)):
                cluster_i, cluster_j = newptn[iclust], newptn[jclust]
                if cluster_i is None or cluster_j is None:  # when we merge two clusters, the merged cluster goes where one was, and the other is replaced by None
                    continue
                overlap_ids = sorted(set(cluster_i) & set(cluster_j))
                if len(overlap_ids) == 0:
                    continue
                newptn[iclust] = cluster_i + [u for u in cluster_j if u not in overlap_ids]  # keep order in i^th cluster, remove common ids from j^th cluster (NOTE this corresponds to the ordering in combine_events())
                if debug:
                    print('        %3d %3d   %s   %s   %4d     %4d   %s   %s' % (iclust, jclust,
                                                                                 color('red' if cluster_i==newptn[iclust] else None, str(len(cluster_i)), width=4),
                                                                                 color('red' if cluster_j==newptn[iclust] else None, str(len(cluster_j)), width=4),
                                                                                 len(newptn[iclust]), len(overlap_ids),
                                                                                 ':'.join(color('red' if u in overlap_ids else None, u) for u in cluster_i),
                                                                                 ':'.join(color('red' if u in overlap_ids else None, u) for u in cluster_j)))
                assert set(newptn[iclust]) == set(cluster_i) | set(cluster_j)
                if antn_list is not None:
                    try_to_update_annotations(cluster_i, cluster_j, newptn[iclust])
                newptn[jclust] = None
                n_changes += 1
        return n_changes
    # ----------------------------------------------------------------------------------------
    from . import clusterpath
    if antn_list is not None:
        antn_dict = get_annotation_dict(antn_list, ignore_duplicates=True, quiet=True)  # we expect lots of duplicates so don't need it to print the duplicate warning
    if debug:
        print('    %smerging overlapping clusters in %d partition%s' % (dbgstr, len(partitions), plural(len(partitions))))
    new_partitions = []
    for ipart in range(len(partitions)):
        if debug:
            unique_uids, total_uids = len(set(u for c in partitions[ipart] for u in c)), sum(len(c) for c in partitions[ipart])
            print('      ipart %d with %d clusters: %d unique vs %d total uids' % (ipart, len(partitions[ipart]), unique_uids, total_uids))
        newptn = [copy.copy(c) for c in partitions[ipart]]
        if debug:
            clusterpath.ptnprint(partitions[ipart], abbreviate=False)
            print('               (%s)             (%s)' % (color('red', 'same as new cluster'), color('red', 'uids in common')))
            print('          i   j   len(i) len(j) len(new)  common')
        n_changes = None
        while n_changes is None or n_changes > 0:
            if debug and n_changes is not None:
                print('      overlap merging: re-processing partition since previous n_changes %d > 0' % n_changes)
            n_changes = process_partition(newptn)
        newptn = [c for c in newptn if c is not None]
        new_partitions.append(newptn)
        if debug:
            clusterpath.ptnprint(newptn, abbreviate=False)
        print('    merged overlapping annotations in partition (N clusters %d --> %d, total cluster size %d --> %d, N unique seqs %d --> %d)' % (len(partitions[ipart]), len(newptn), sum(len(c) for c in partitions[ipart]), sum(len(c) for c in newptn), len(set(u for c in partitions[ipart] for u in c)), len(set(u for c in newptn for u in c))))
    # return new_partitions, [antn_dict[':'.join(c)] for p in new_partitions for c in p if ':'.join(c) in antn_dict]  # maybe it's worth putting the annotations in some order?
    if antn_list is not None:
        check_concordance_of_cpath_and_annotations(clusterpath.ClusterPath(partition=new_partitions[0]), antn_dict.values(), antn_dict)
    return new_partitions, list(antn_dict.values()) if antn_list is not None else None

# ----------------------------------------------------------------------------------------
def build_dummy_reco_info(true_partition):
    def tkey(c): return ':'.join(c)
    chashes = {tkey(tc) : uidhashstr(tkey(tc)) for tc in true_partition}
    return {u : {'reco_id' : chashes[tkey(tc)]} for tc in true_partition for u in tc}

# ----------------------------------------------------------------------------------------
def per_seq_correct_cluster_fractions(partition, true_partition, reco_info=None, seed_unique_id=None, dbg_str='', inf_label='inferred', true_label='true', debug=False):
    if seed_unique_id is None:
        check_intersection_and_complement(partition, true_partition, a_label=inf_label, b_label=true_label)
    if reco_info is None:  # build a dummy reco_info that just has reco ids
        reco_info = build_dummy_reco_info(true_partition)
    reco_ids = {uid : reco_info[uid]['reco_id'] for cluster in partition for uid in cluster}  # NOTE this iterates over the *inferred* partition, which maybe is important if seed id is set? (i don't remember atm)
    uids = set([uid for cluster in partition for uid in cluster])
    clids = get_cluster_ids(uids, partition)  # map of {uid : (index of cluster in <partition> in which that uid occurs)} (well, list of indices, in case there's duplicates)

    def get_clonal_fraction(uid, inferred_cluster):
        """ Return the fraction of seqs in <uid>'s inferred cluster which are really clonal. """
        n_clonal = 0
        for tmpid in inferred_cluster:  # NOTE this includes the case where tmpid equal to uid
            if reco_ids[tmpid] == reco_ids[uid]:
                n_clonal += 1
        return float(n_clonal) / len(inferred_cluster)

    def get_fraction_present(inferred_cluster, true_cluster):
        """ Return the fraction of the true clonemates in <true_cluster> that appear in <inferred_cluster>. """
        return len(set(inferred_cluster) & set(true_cluster)) / float(len(true_cluster))

    mean_clonal_fraction, mean_fraction_present = 0., 0.
    n_uids = 0
    for true_cluster in true_partition:
        if seed_unique_id is not None and seed_unique_id not in true_cluster:
            continue
        for uid in true_cluster:
            if seed_unique_id is not None and uid != seed_unique_id:
                continue
            if len(clids[uid]) != 1:  # this seems to only happen for earlier partitions (more than one proc) when seed_unique_id is set, since we pass seed_unique_id to all the subprocs. I.e. it's expected in these cases, and the ccfs don't make sense when a uid is in more than one cluster, since it's no longer a partition, so just return None, None
                if debug:
                    print('  %s found %s in multiple clusters while calculating ccfs (returning None, None)' % (color('red', 'warning'), uid))
                return None, None
            inferred_cluster = partition[clids[uid][0]]  # we only look at the first cluster in which it appears
            mean_clonal_fraction += get_clonal_fraction(uid, inferred_cluster)
            mean_fraction_present += get_fraction_present(inferred_cluster, true_cluster)
            n_uids += 1

    if n_uids > 1e6:
        raise Exception('you should start worrying about numerical precision if you\'re going to run on this many queries')

    if debug:
        print('    %scorrect cluster fractions:' % dbg_str)
        print('       clusters uids')
        print('         %3d    %3d  %s' % (len(true_partition), sum(len(c) for c in true_partition), true_label))
        print('         %3d    %3d  %s' % (len(partition), sum(len(c) for c in partition), inf_label))
        print('      purity:       %.1f / %d = %.3f' % (mean_clonal_fraction, n_uids, mean_clonal_fraction / n_uids))
        print('      completeness: %.1f / %d = %.3f' % (mean_fraction_present, n_uids, mean_fraction_present / n_uids))
    return mean_clonal_fraction / n_uids, mean_fraction_present / n_uids

# ----------------------------------------------------------------------------------------
def per_family_correct_cluster_fractions(partition, true_partition, debug=False):
    # The new ccfs above are pretty similar, except they're per-sequence rather than per-cluster, so they don't get all scatterbrained and shit when a sample's only got a few clusters.
    # Also, you still get partial credit for how good your cluster is, it's not just all-or-nothing.
    raise Exception('deprecated! use per_seq_correct_cluster_fractions() above')

    def find_clusters_with_ids(ids, partition):
        """ find all clusters in <partition> that contain at least one of <ids> """
        clusters = []
        for cluster in partition:
            for uid in ids:
                if uid in cluster:
                    clusters.append(cluster)
                    break
        return clusters

    check_intersection_and_complement(partition, true_partition)

    n_under_merged, n_over_merged = 0, 0
    for trueclust in true_partition:
        if debug:
            print('')
            print('   true %s' % (len(trueclust) if len(trueclust) > 15 else trueclust))
        infclusters = find_clusters_with_ids(trueclust, partition)  # list of inferred clusters that contain any ids from the true cluster
        if debug and len(infclusters) > 1:
            print('  infclusters %s' % infclusters)
        assert len(infclusters) > 0
        under_merged = len(infclusters) > 1  # ids in true cluster are not all in the same inferred cluster
        over_merged = False  # at least one inferred cluster with an id in true cluster also contains an id not in true cluster
        for iclust in infclusters:
            if debug:
                print('   inferred %s' % (len(iclust) if len(iclust) > 15 else iclust))
            for uid in iclust:
                if uid not in trueclust:
                    over_merged = True
                    break
            if over_merged:
                break
        if debug:
            print('  under %s   over %s' % (under_merged, over_merged))
        if under_merged:
            n_under_merged += 1
        if over_merged:
            n_over_merged += 1

    under_frac = float(n_under_merged) / len(true_partition)
    over_frac = float(n_over_merged) / len(true_partition)
    if debug:
        print('  under %.2f   over %.2f' % (under_frac, over_frac))
    return (1. - under_frac, 1. - over_frac)

# ----------------------------------------------------------------------------------------
# from definitions in mobille paper
# problems with these pairwise metrics:
#  - they scale quadratically, which (subjectively) I think is counter intuitive to how our brains work
#    - because of this, large clusters matter much more than small clusters, e.g. a two-sequence true cluster gets "counted" twice, while a 10-sequence one gets counted 100 times
#  - they're both 0 (or nan) on the singleton partition
#  - they generally just kind of ignore singletons, e.g. precision is supposed to measure overmerging, but it gives you no credit for not over merging singletons
def pairwise_cluster_metrics(mtstr, inf_ptn, tru_ptn, debug=False):
    debug = True
    # ----------------------------------------------------------------------------------------
    def id_dict(ptn):
        reco_info = build_dummy_reco_info(ptn)  # not actually reco info unless it's the true partition
        return {uid : reco_info[uid]['reco_id'] for cluster in ptn for uid in cluster}  # speed optimization
    # ----------------------------------------------------------------------------------------
    # inf_ptn, tru_ptn = [['c'], ['a', 'b', 'e'], ['d', 'g', 'f']], [['a', 'b', 'c'], ['d', 'e', 'f', 'g']]  # example from paper, should be pairwise: (0.666666, 0.44444444, ?), closeness: (0.857, 0.6, ?) (see note below, they calculate it wrong)
    # inf_ptn, tru_ptn = [['a'], ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']], [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n']]
    check_intersection_and_complement(inf_ptn, tru_ptn, a_label='true', b_label='inferred')
    if mtstr == 'pairwise':
        tp, fp, fn, n_tot = 0, 0, 0, 0
        tru_ids, inf_ids = [id_dict(ptn) for ptn in [tru_ptn, inf_ptn]]
        for u1, u2 in itertools.combinations(set(u for c in tru_ptn for u in c), 2):
            is_tru_clonal, is_inf_clonal = [tids[u1] == tids[u2] for tids in [tru_ids, inf_ids]]
            n_tot += 1
            if is_tru_clonal and is_inf_clonal:
                tp += 1
            elif is_tru_clonal:
                fn += 1
            elif is_inf_clonal:
                fp += 1
            else:  # singletons
                pass
    elif mtstr == 'closeness':
        tp, fp, fn, n_tot = set(), set(), set(), set()
        if debug:
            print('    infcl   trucl      tp   fp     fn')
        # for infcl, trucl in itertools.product(inf_ptn, tru_ptn):  # NOTE this is *not* what they mean, see next line
        for infcl in inf_ptn:  # NOTE this is apparently what they mean by "we first identified the best correspondence between inferred clonal lineages and correct clonal assignments" but note you'd get a *different* answer if you looped over tru_ptn
            trucl = sorted(tru_ptn, key=lambda c: len(set(c) & set(infcl)), reverse=True)[0]
            infset, truset = set(infcl), set(trucl)
            if len(infset & truset) == 0:
                continue
            tp |= infset & truset  # OMFG their example figure is wrong, it (correctly) shows 6 entries, but then when they calculate recall they switch it to 5
            fp |= infset - truset
            fn |= truset - infset
            n_tot |= truset
            if debug:
                print('  %20s   %20s   %20s   %20s   %20s' % (infcl, trucl, infset & truset, infset - truset, truset - infset))
        tp, fp, fn = [len(s) for s in [tp, fp, fn]]
    else:
        assert False
    precis = tp / float(tp + fp) if tp + fp > 0 else float('nan')
    recall = tp / float(tp + fn) if tp + fn > 0 else float('nan')  # same as sensitivity
    if debug:
        print('    pairwise clustering metrics:')
        print('        precision: tp / (tp + fp) = %d / (%d + %d) = %.2f' % (tp, tp, fp, precis))
        print('      sens/recall: tp / (tp + fn) = %d / (%d + %d) = %.2f' % (tp, tp, fn, recall))
    return {'precision' : precis, 'recall' : recall, 'f1' : 2 * precis * recall / float(precis + recall) if precis + recall > 0 else float('nan')}  # same as scipy.stats.hmean([precis, recall])

# ----------------------------------------------------------------------------------------
# return (# of within-cluster seq pairs) / (total # of seq pairs, i.e. n*(n-1)/2), i.e. if a "collision" is that two seqs are in a cluster together, this counts the number of actual collided sequence pairs, over the total number of possible collisions
def collision_fraction(partition):
# TODO this is only actually right on singleton true partitions.
# I think what i really want is two numbers:
#   - "naive" collision fraction: fraction of naive rearrangement pairs that collide (measures density/inverse diversity of naive repertoire)
#   - "mature" collision fraction: fraction of mature/mutated seq pairs from *different* families that collide (also adds in extra collisions from shm)
    def n_combos(n):
        return n * (n - 1) // 2
    n_collisions = sum(n_combos(len(c)) for c in partition)
    cfrac = float(n_collisions) / n_combos(sum(len(c) for c in partition))
    return cfrac

# ----------------------------------------------------------------------------------------
# NOTE n_biggest_clusters can be a pair of numbers
def partition_similarity_matrix(partition_a, partition_b, n_biggest_clusters, iscn_denominator='max', a_label='a', b_label='b', debug=False):
    # iscn_denominator: denominator to divide intersection size of each pair of clusters ('min': min of the two sizes, 'mean': mean of the two sizes)
    """ Return matrix whose ij^th entry is the size of the intersection between <partition_a>'s i^th biggest cluster and <partition_b>'s j^th biggest, divided by the mean/min/max size of the two clusters """
    # ----------------------------------------------------------------------------------------
    def sort_within_clusters(part):  # sort each cluster's uids alphabetically (so we can sort clusters alphabetically below in order to get a more consistent sorting between partitions, although maybe it would make more sense to sort afterwards by intersection size)
        for iclust in range(len(part)):
            part[iclust] = sorted(part[iclust])
    # ----------------------------------------------------------------------------------------
    def norm_factor(clust_a, clust_b):
        if iscn_denominator == 'min':
            return min(len(clust_a), len(clust_b))
        elif iscn_denominator == 'max':
            return max(len(clust_a), len(clust_b))
        elif iscn_denominator == 'mean':
            return 0.5 * (len(clust_a) + len(clust_b))  # mean size of the two clusters
        else:
            assert False
    # ----------------------------------------------------------------------------------------
    def getsclusts(ptn, n_big):
        sort_within_clusters(ptn)  # have to do this before taking N biggest so that when there's lots of ties in size we're more likely to get the same clusters for both methods
        srtclusts = sorted(sorted(ptn), key=len, reverse=True)
        n_to_use = min(n_big, len(srtclusts))
        plot_clusts = srtclusts[ : n_to_use]
        skip_clusts = srtclusts[n_to_use : ]
        return plot_clusts, skip_clusts
    # ----------------------------------------------------------------------------------------
    def incr_clust(csize, n_common, ifrac):
        sub_szs.append(csize)
        sub_ovlps.append(n_common)
        sub_fracs.append(ifrac)
    # ----------------------------------------------------------------------------------------
    check_intersection_and_complement(partition_a, partition_b, a_label=a_label, b_label=b_label, only_warn=True)
    n_big_a, n_big_b = n_biggest_clusters if hasattr(n_biggest_clusters, '__iter__') else (n_biggest_clusters, n_biggest_clusters)
    a_clusters, a_skip_clusts = getsclusts(partition_a, n_big_a)
    b_clusters, b_skip_clusts = getsclusts(partition_b, n_big_b)

    smatrix = [[float('nan') for _ in b_clusters] for _ in a_clusters]
    dszs, dovlps, dfracs = [], [], []
    overlap_matrix = [[set(clust_a) & set(clust_b) for clust_b in b_clusters] for clust_a in a_clusters]
    missing_b_clusts = list(range(len(b_clusters)))  # need to look for any clusters in b that overlapped only with small a clusters (i.e. won't otherwise show up in the debug printout)
    for ia, clust_a in enumerate(a_clusters):
        sub_szs, sub_ovlps, sub_fracs, n_found = [], [], [], 0
        for ib, clust_b in enumerate(b_clusters):  # first look just in the precalculated overlap matrix
            n_common = len(overlap_matrix[ia][ib])
            if n_common == 0:
                continue
            ifrac = float(n_common) / norm_factor(clust_a, clust_b)
            smatrix[ia][ib] = float('nan') if ifrac==0 else ifrac  # nan gives you transparent/empty color
            if debug:
                n_found += n_common
                incr_clust(len(clust_b), n_common, ifrac)
                if ib in missing_b_clusts:
                    missing_b_clusts.remove(ib)
        if debug and n_found < len(clust_a):  # if we didn't find any with overlap, look in smaller clusters (I think doing this when debug isn't set would be annoyingly slow. This fcn is meant to plot a summary of the largest clusters on potentially large repertoires, if you want really comprehensive info, just set debug)
            for bclst in b_skip_clusts:
                ncm = len(set(clust_a) & set(bclst))
                if ncm > 0:
                    n_found += ncm
                    incr_clust(len(bclst), ncm, float(ncm) / norm_factor(clust_a, bclst))
                sub_szs = sorted(sub_szs, reverse=True)
                sub_ovlps = sorted(sub_ovlps, reverse=True)
            if n_found != len(clust_a):
                print('  %s couldn\'t find %d / %d uids' % (wrnstr(), len(clust_a) - n_found, len(clust_a)))
        if debug:
            dszs.append(sub_szs)
            dovlps.append(sub_ovlps)
            dfracs.append(sub_fracs)
    if debug:
        max_dsz, max_ovl = max(len(s) for s in dszs), max(len(o) for o in dovlps)
        print('    cluster sizes (%s):' % color('blue', 'used'))
        for tlab, aclusts, eclusts in zip([a_label, b_label], [partition_a, partition_b], [a_clusters, b_clusters]):
            print('    %12s: %s' % (tlab, cluster_size_str(aclusts, clusters_to_emph=eclusts)))
        print('    %s  %s -->' % (a_label, b_label))
        print('     size   %s    %s      overlap / %s size (*100)' % (wfmt('sizes', 2*max_dsz, jfmt='-'), wfmt('overlaps', 2*max_ovl, jfmt='-'), iscn_denominator.replace('min', 'smaller').replace('max', 'larger')))
        cl = len(str(max(len(a_clusters[0]), len(b_clusters[0]))))
        for ia, clust_a in enumerate(a_clusters):
            szstr, ovlstr = ' '.join('%d'%s for s in dszs[ia]), ' '.join('%d'%o for o in dovlps[ia])
            print(('    %4d      %-s    %-s    %s') % (len(clust_a), wfmt(szstr, 2*max_dsz, jfmt='-'), wfmt(ovlstr, 2*max_ovl, jfmt='-'), ' '.join('%3.0f'%(100*f) for f in dfracs[ia])))
        if len(missing_b_clusts) > 0:
            print('  %s %d %s clusters only overlapped with small %s clusters, i.e. won\'t appear in debug output above (if there\'s large clusters here this means you should probably increase the number of clusters):' % (wrnstr(), len(missing_b_clusts), b_label, a_label))
            print('     %s' % cluster_size_str([b_clusters[i] for i in missing_b_clusts]))

    a_cluster_lengths, b_cluster_lengths = [len(c) for c in a_clusters], [len(c) for c in b_clusters]  # nice to return these here so you don't have to re-sort in the calling fcn
    return a_cluster_lengths, b_cluster_lengths, smatrix

# ----------------------------------------------------------------------------------------
def find_uid_in_partition(uid, partition):
    iclust, found = None, False
    for iclust in range(len(partition)):
        if uid in partition[iclust]:
            found = True
            break
    if found:
        return iclust
    else:
        raise Exception('couldn\'t find %s in %s\n' % (uid, partition))

# ----------------------------------------------------------------------------------------
def check_intersection_and_complement(part_a, part_b, only_warn=False, a_label='a', b_label='b', print_uids=False, debug=False):
    """ make sure two partitions have identical uid lists """
    uids_a, uids_b = [ptn_ids(p) for p in [part_a, part_b]]
    a_and_b = uids_a & uids_b
    a_not_b = uids_a - uids_b
    b_not_a = uids_b - uids_a
    if len(a_not_b) > 0 or len(b_not_a) > 0:  # NOTE this should probably also warn/pring if either of 'em has duplicate uids on their own
        failstr = '\'%s\' partition (%d total) and \'%s\' partition (%d total) don\'t have the same uids:   only %s %d    only %s %d    common %d' % (a_label, sum(len(c) for c in part_a), b_label, sum(len(c) for c in part_b), a_label, len(a_not_b), b_label, len(b_not_a), len(a_and_b))
        if print_uids:
            failstr += '\n' + '\n'.join(['    only %s: %s' % (a_label, ' '.join(a_not_b)),
                                         '    only %s: %s' % (b_label, ' '.join(b_not_a)),
                                         '    common: %s' % ' '.join(a_and_b)])
        if only_warn:
            print('  %s %s' % (color('yellow', 'warning'), failstr))
            if debug:
                for tl, tids in [(a_label, a_not_b), (b_label, b_not_a)]:
                    print('     only %s: %s' % (tl, ' '.join(tids)))
        else:
            raise Exception(failstr)
    else:
        if debug:
            print('  intersection and complement both ok')
    return a_and_b, a_not_b, b_not_a

# ----------------------------------------------------------------------------------------
def ptn_ids(partition):  # return set of uids in <partition> (adding this very late, so could probably be used in a lot of places it isn't)
    return set(u for c in partition for u in c)

# ----------------------------------------------------------------------------------------
def get_cluster_list_for_sklearn(part_a, part_b):
    # convert from partition format {cl_1 : [seq_a, seq_b], cl_2 : [seq_c]} to [cl_1, cl_1, cl_2]
    # NOTE this will be really slow for larger partitions

    # first make sure that <part_a> has every uid in <part_b> (the converse is checked below)
    for jclust in range(len(part_b)):
        for uid in part_b[jclust]:
            find_uid_in_partition(uid, part_a)  # raises exception if not found

    # then make the cluster lists
    clusts_a, clusts_b = [], []
    for iclust in range(len(part_a)):
        for uid in part_a[iclust]:
            clusts_a.append(iclust)
            clusts_b.append(find_uid_in_partition(uid, part_b))

    return clusts_a, clusts_b

# ----------------------------------------------------------------------------------------
def adjusted_mutual_information(partition_a, partition_b):
    return -1.  # not using it any more, and it's really slow
    # clusts_a, clusts_b = get_cluster_list_for_sklearn(partition_a, partition_b)
    # return sklearn.metrics.cluster.adjusted_mutual_info_score(clusts_a, clusts_b)

# ----------------------------------------------------------------------------------------
def add_missing_uids_to_partition(partition_with_missing_uids, ref_partition, allowed_uids=None, miss_label='', ref_label='', warn=False, fail_frac=None, dbgstr='', debug=False):
    """ return a copy of <partition_with_missing_uids> which has had any uids which were missing inserted as singletons (i.e. uids which were in <ref_partition>) """
    ref_ids, pmiss_ids = set(u for c in ref_partition for u in c), set(u for c in partition_with_missing_uids for u in c)
    missing_ids = ref_ids - pmiss_ids
    if allowed_uids is not None:  # if <allowed_uids> is set, we're only allowed to add those uids
        missing_ids &= set(allowed_uids)
    partition_with_uids_added = copy.deepcopy(partition_with_missing_uids) + [[u] for u in missing_ids]
    if len(missing_ids) > 0 and (debug or warn or fail_frac):
        if debug or warn:
            print('    %sadded %d / %d missing ids as singletons to %s partition, to match ids in %s partition%s' % (wrnstr()+' ' if warn else '', len(missing_ids), sum(len(c) for c in ref_partition), miss_label, ref_label, dbgstr)) #, ' '.join(missing_ids) if len(missing_ids) < 30 else '')
        miss_frac = len(missing_ids) / float(sum(len(c) for c in ref_partition))
        if fail_frac is not None and miss_frac > fail_frac:
            raise Exception('missing %.3f (more than %.3f, increase --max-ccf-fail-frac to avoid this error)' % (miss_frac, fail_frac))
    return partition_with_uids_added

# ----------------------------------------------------------------------------------------
def remove_missing_uids_from_partition(ref_partition, partition_with_missing_uids, ref_label='ref', miss_label='miss', fail_frac=None, dbgstr='', warn=False):
    """ return a copy of <ref_partition> which has had any uids which do not occur in <partition_with_missing_uids> removed """
    ref_ids, pmiss_ids = set(u for c in ref_partition for u in c), set(u for c in partition_with_missing_uids for u in c)
    missing_ids = ref_ids - pmiss_ids
    ref_partition_with_uids_removed = [[u for u in c] for c in ref_partition]
    for iclust, cluster in enumerate(ref_partition_with_uids_removed):
        ref_partition_with_uids_removed[iclust] = [u for u in cluster if u not in missing_ids]
    ref_partition_with_uids_removed = [c for c in ref_partition_with_uids_removed if len(c) > 0]  # remove any empty clusters
    n_remaining = sum([len(c) for c in ref_partition_with_uids_removed])
    if len(missing_ids) > 0:
        wstr = (color('yellow', 'warning')+' ') if (n_remaining==0 or warn) else ''
        miss_frac = len(missing_ids) / float(sum([len(c) for c in ref_partition]))
        print('    %sremoved %d/%d uids (leaving %d) from %s partition that were not in %s partition%s' % (wstr, len(missing_ids), sum([len(c) for c in ref_partition]), n_remaining, ref_label, miss_label, dbgstr)) #, ' '.join(missing_ids))
        if fail_frac is not None and miss_frac > fail_frac:
            raise Exception('missing %.3f (more than %.3f)' % (miss_frac, fail_frac))
    return ref_partition_with_uids_removed

# ----------------------------------------------------------------------------------------
def generate_incorrect_partition(true_partition, misassign_fraction, error_type, debug=False):
    """
    Generate an incorrect partition from <true_partition>.
    We accomplish this by removing <n_misassigned> seqs at random from their proper cluster, and putting each in either a
    cluster chosen at random from the non-proper clusters (<error_type> 'reassign') or in its own partition (<error_type> 'singleton').
    """
    # true_partition = [['b'], ['a', 'c', 'e', 'f'], ['d', 'g'], ['h', 'j']]
    # debug = True

    new_partition = copy.deepcopy(true_partition)
    nseqs = sum([len(c) for c in true_partition])
    n_misassigned = int(misassign_fraction * nseqs)
    if debug:
        print('  misassigning %d / %d seqs (should be clsoe to %.3f)' % (n_misassigned, nseqs, misassign_fraction))
        print('  before', new_partition)

    uids = [uid for cluster in true_partition for uid in cluster]
    for _ in range(n_misassigned):
        uid = uids[random.randint(0, len(uids) - 1)]  # choose a uid to misassign (note that randint() is inclusive)
        uids.remove(uid)
        iclust = find_uid_in_partition(uid, new_partition)
        new_partition[iclust].remove(uid)  # remove it
        if [] in new_partition:
            new_partition.remove([])
        if error_type == 'singletons':  # put the sequence in a cluster by itself
            new_partition.append([uid, ])
            if debug:
                print('    %s: %d --> singleton' % (uid, iclust))
        elif error_type == 'reassign':  # choose a different cluster to add it to
            inewclust = iclust
            while inewclust == iclust:  # hm, this won't work if there's only one cluster in the partition. Oh, well, that probably won't happen
                inewclust = random.randint(0, len(new_partition) - 1)
            new_partition[inewclust].append(uid)
            if debug:
                print('    %s: %d --> %d' % (uid, iclust, inewclust))
        else:
            raise Exception('%s not among %s' % (error_type, 'singletons, reassign'))
    if debug:
        print('  after', new_partition)
    return new_partition

# ----------------------------------------------------------------------------------------
def subset_files(uids, fnames, outdir, uid_header='Sequence ID', delimiter='\t', debug=False):
    """ rewrite csv files <fnames> to <outdir>, removing lines with uid not in <uids> """
    for fname in fnames:
        with open(fname) as infile:
            reader = csv.DictReader(infile, delimiter=str(delimiter))
            with open(outdir + '/' + os.path.basename(fname), csv_wmode()) as outfile:
                writer = csv.DictWriter(outfile, reader.fieldnames, delimiter=str(delimiter))
                writer.writeheader()
                for line in reader:
                    if line[uid_header] in uids:
                        writer.writerow(line)

# ----------------------------------------------------------------------------------------
def csv_to_fasta(infname, outfname=None, name_column='unique_ids', seq_column='input_seqs', n_max_lines=None, overwrite=True, remove_duplicates=False, debug=True):

    if not os.path.exists(infname):
        raise Exception('input file %s d.n.e.' % infname)
    if outfname is None:
        assert '.csv' in infname
        outfname = infname.replace('.csv', '.fa')
    if os.path.exists(outfname):
        if overwrite:
            if debug:
                print('  csv --> fasta: overwriting %s' % outfname)
        else:
            if debug:
                print('  csv --> fasta: leaving existing outfile %s' % outfname)
            return

    if '.csv' in infname:
        delimiter = ','
    elif '.tsv' in infname:
        delimiter = '\t'
    else:
        assert False

    uid_set = set()
    n_duplicate_ids = 0
    with open(infname) as infile:
        reader = csv.DictReader(infile, delimiter=str(delimiter))
        with open(outfname, 'w') as outfile:
            n_lines = 0
            for line in reader:
                if seq_column not in line:
                    raise Exception('specified <seq_column> \'%s\' not in line (keys in line: %s)' % (seq_column, ' '.join(list(line.keys()))))
                if name_column is not None:
                    if name_column not in line:
                        raise Exception('specified <name_column> \'%s\' not in line (keys in line: %s)' % (name_column, ' '.join(list(line.keys()))))
                    uid = line[name_column]
                else:
                    uid = uidhashstr(line[seq_column])
                if remove_duplicates:
                    if uid in uid_set:
                        n_duplicate_ids += 1
                        continue
                    uid_set.add(uid)
                n_lines += 1
                if n_max_lines is not None and n_lines > n_max_lines:
                    break
                outfile.write('>%s\n' % uid)
                outfile.write('%s\n' % line[seq_column])
    if debug and n_duplicate_ids > 0:
        print('   skipped %d / %d duplicate uids' % (n_duplicate_ids, len(uid_set)))

# ----------------------------------------------------------------------------------------
def print_heapy(extrastr, heap):
    'Partition of a set of 1511530 objects. Total size = 188854824 bytes.'
    heapstr = heap.__str__()
    total = None
    for line in heapstr.split('\n'):
        if 'Total size' in line:
            total = int(line.split()[10])
    if total is None:
        print('oops')
        print(heapstr)
        sys.exit()
    print('mem total %.3f MB    %s' % (float(total) / 1e6, extrastr))

# ----------------------------------------------------------------------------------------
def auto_slurm(n_procs):
    """ Return true if we want to force slurm usage, e.g. if there's more processes than cores """

    def slurm_exists():
        try:
            fnull = open(os.devnull, 'w')
            subprocess.check_output(['which', 'srun'], stderr=fnull, close_fds=True, universal_newlines=True)
            return True
        except subprocess.CalledProcessError:
            return False

    ncpu = multiprocessing.cpu_count()
    if n_procs > ncpu and slurm_exists():
        return True
    return False

# ----------------------------------------------------------------------------------------
def add_regional_alignments(glfo, aligned_gl_seqs, line, region, debug=False):
    if debug:
        print(' %s' % region)

    aligned_seqs = [None for _ in range(len(line['unique_ids']))]
    for iseq in range(len(line['seqs'])):
        qr_seq = line[region + '_qr_seqs'][iseq]
        gl_seq = line[region + '_gl_seq']
        aligned_gl_seq = aligned_gl_seqs[region][line[region + '_gene']]
        if len(qr_seq) != len(gl_seq):
            if debug:
                print('    qr %d and gl %d seqs different lengths for %s, setting invalid' % (len(qr_seq), len(gl_seq), ' '.join(line['unique_ids'])))
            line['invalid'] = True
            continue

        n_gaps = gap_len(aligned_gl_seq)
        if n_gaps == 0:
            if debug:
                print('   no gaps')
            aligned_seqs[iseq] = qr_seq
            continue

        if debug:
            print('   before alignment')
            print('      qr   ', qr_seq)
            print('      gl   ', gl_seq)
            print(' aligned gl', aligned_gl_seq)

        # add dots for 5p and 3p deletions
        qr_seq = gap_chars[0] * line[region + '_5p_del'] + qr_seq + gap_chars[0] * line[region + '_3p_del']
        gl_seq = gap_chars[0] * line[region + '_5p_del'] + gl_seq + gap_chars[0] * line[region + '_3p_del']

        if len(aligned_gl_seq) - n_gaps != len(gl_seq):
            if debug:
                print('    aligned germline seq without gaps (%d - %d = %d) not the same length as unaligned gl/qr seqs %d' % (len(aligned_gl_seq), n_gaps, len(aligned_gl_seq) - n_gaps, len(gl_seq)))
            line['invalid'] = True
            continue

        qr_seq = list(qr_seq)
        gl_seq = list(gl_seq)
        for ibase in range(len(aligned_gl_seq)):
            if aligned_gl_seq[ibase] in gap_chars:  # add gap to the qr and gl seq lists
                qr_seq.insert(ibase, gap_chars[0])
                gl_seq.insert(ibase, gap_chars[0])
            elif gl_seq[ibase] == aligned_gl_seq[ibase] or gl_seq[ibase] in gap_chars:  # latter is 5p or 3p deletion that we filled in above
                pass  # all is well
            else:  # all is not well, don't know why
                line['invalid'] = True
                break
        if line['invalid']:
            if debug:
                print('    unknown error during alignment process')
            continue
        qr_seq = ''.join(qr_seq)
        gl_seq = ''.join(gl_seq)

        if debug:
            print('   after alignment')
            print('      qr   ', qr_seq)
            print('      gl   ', gl_seq)
            print(' aligned gl', aligned_gl_seq)

        if len(qr_seq) != len(gl_seq) or len(qr_seq) != len(aligned_gl_seq):  # I don't think this is really possible as currently written
            if debug:
                print('    lengths qr %d gl %d and aligned gl %d not all the same after alignment' % (len(qr_seq), len(gl_seq), len(aligned_gl_seq)))
            line['invalid'] = True
            continue

        aligned_seqs[iseq] = qr_seq

    if line['invalid']:
        print('%s failed adding alignment info for %s' % (color('red', 'error'),' '.join(line['unique_ids'])))  # will print more than once if it doesn't fail on the last region
        aligned_seqs = [None for _ in range(len(line['seqs']))]

    line['aligned_' + region + '_seqs'] = aligned_seqs

# ----------------------------------------------------------------------------------------
def add_alignments(glfo, aligned_gl_seqs, line, debug=False):
    for region in regions:
        add_regional_alignments(glfo, aligned_gl_seqs, line, region, debug)

# ----------------------------------------------------------------------------------------
def intexterpolate(x1, y1, x2, y2, x):
    """ interpolate/extrapolate linearly based on two points in 2-space, returning y-value corresponding to <x> """
    if x1 == x2:
        print('  %s x1 equal to x2 in utils.intexterpolate()' % wrnstr())
        return 0  # arg, maybe this is right?
    m = (y2 - y1) / float(x2 - x1)
    b = 0.5 * (y1 + y2 - m*(x1 + x2))
    # if debug:
    #     for x in [x1, x2]:
    #         print '%f x + %f = %f' % (m, b, m*x + b)
    return m * x + b

# ----------------------------------------------------------------------------------------
def get_mean_mfreq(parameter_dir):
    from .hist import Hist
    mutehist = Hist(fname=parameter_dir + '/all-mean-mute-freqs.csv')
    return mutehist.get_mean(ignore_overflows=True)  # should I not ignore overflows here?

# ----------------------------------------------------------------------------------------
def get_naive_hamming_bounds(partition_method, parameter_dir=None, overall_mute_freq=None):  # parameterize the relationship between mutation frequency and naive sequence inaccuracy
    if parameter_dir is not None:
        assert overall_mute_freq is None
        mute_freq = get_mean_mfreq(parameter_dir)
    else:
        assert overall_mute_freq is not None
        mute_freq = overall_mute_freq

    # just use a line based on two points (mute_freq, threshold)
    x1, x2 = 0.05, 0.2  # 0.5x, 3x (for 10 leaves)

    if partition_method == 'naive-hamming':  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
        y1, y2 = 0.035, 0.06
        lo = intexterpolate(x1, y1, x2, y2, mute_freq)
        hi = lo
    elif partition_method == 'naive-vsearch':  # set lo and hi to the same thing, so we don't use log prob ratios, i.e. merge if less than this, don't merge if greater than this
        y1, y2 = 0.02, 0.05
        lo = intexterpolate(x1, y1, x2, y2, mute_freq)
        hi = lo
    elif partition_method == 'likelihood':  # these should almost never merge non-clonal sequences or split clonal ones, i.e. they're appropriate for naive hamming preclustering if you're going to then run the full likelihood (i.e. anything less than lo is almost certainly clonal, anything greater than hi is almost certainly not)
        y1, y2 = 0.015, 0.015  # would be nice to get better numbers for this
        lo = intexterpolate(x1, y1, x2, y2, mute_freq)  # ...and never merge 'em if it's bigger than this
        y1, y2 = 0.08, 0.15
        hi = intexterpolate(x1, y1, x2, y2, mute_freq)  # ...and never merge 'em if it's bigger than this
    else:
        assert False

    return lo, hi

# ----------------------------------------------------------------------------------------
def find_genes_that_have_hmms(parameter_dir):
    yamels = glob.glob(parameter_dir + '/hmms/*.yaml')
    if len(yamels) == 0:
        raise Exception('no yamels in %s' % parameter_dir + '/hmms')

    genes = []
    for yamel in yamels:
        gene = unsanitize_name(os.path.basename(yamel).replace('.yaml', ''))
        genes.append(gene)

    return genes

# ----------------------------------------------------------------------------------------
# takes the first uid from the first cluster in the best partition that meets the size criteria
#   - for paired, pass in the output dir
# <iseed>: skip this many clusters that meet the other criteria (then choose the next one)
def choose_seed_unique_id(infname, min_cluster_size, max_cluster_size, iseed=None, n_max_queries=-1, choose_random=False, paired=False, ig_or_tr='ig', debug=True):
    # ----------------------------------------------------------------------------------------
    def choose_sid(sclust):
        iseq = 0
        if not paired:
            return sclust[iseq], None  # don't think it matters which one
        antn_dict = get_annotation_dict(annotation_list)
        antn = antn_dict[':'.join(sclust)]
        for iseq, uid in enumerate(antn['unique_ids']):
            pids = antn['paired-uids'][iseq]
            if len(pids) == 0:
                continue
            correctly_paired_ids = [p for p in pids if is_correctly_paired(uid, p)]
            if len(correctly_paired_ids) == 0:
                continue
            return uid, get_single_entry(correctly_paired_ids)
        raise Exception('couldn\'t find a seed id in cluster with size %d (none had their correct paired id among their \'paired-uids\')')
    # ----------------------------------------------------------------------------------------
    if paired:
        from . import paircluster  # it barfs if i import at the top, and i know that means i'm doing something dumb, but whatever
        infname = paircluster.paired_fn(infname, heavy_locus(ig_or_tr), suffix='.yaml', actstr='auto')  # heavy chains seqs paired with either light chain
    _, annotation_list, cpath = read_output(infname, n_max_queries=n_max_queries, dont_add_implicit_info=True)
    if choose_random:
        assert iseed is None  # <iseed> if you want to specify a specific cluster, <choose_random> is if you want us to choose one at random from the whole file
        iseed = random.randint(0, len(cpath.best()) - 1)  # inclusive
        if debug:
            print('  chose random iseed: %d (of %d clusters)' % (iseed, len(cpath.best())))

    nth_seed = 0  # numbrer of clusters we've passed that meet the other criteria
    sid, pid, scluster = None, None, None
    for cluster in sorted(cpath.best(), key=len, reverse=True):
        if min_cluster_size is not None and len(cluster) < min_cluster_size:
            continue
        if max_cluster_size is not None and len(cluster) > max_cluster_size:
            continue
        if iseed is not None and int(iseed) > nth_seed:
            nth_seed += 1
            continue
        scluster = cluster
        sid, pid = choose_sid(cluster)
        break
    if sid is None:
        raise Exception('couldn\'t find seed in cluster between size %d and %d' % (min_cluster_size, max_cluster_size))

    if paired:
        plocus = pid.split('-')[-1]
        if plocus not in loci:
            raise Exception('couldn\'t get allowed locus (got \'%s\') from uid \'%s\'' % (plocus, pid))
        if debug:
            print('    chose seed uids %s %s with loci %s %s from heavy cluster with size %d%s%s%s (%s)' % (sid, pid, heavy_locus(ig_or_tr), plocus, len(cluster), ' (chosen at random)' if choose_random else '',
                                                                                                            '' if nth_seed==0 else ', after skipping %d/%d other cluster%s'%(nth_seed, len(cpath.best()), plural(nth_seed)),
                                                                                                            '' if min_cluster_size is None and max_cluster_size is None else ' (was asked for size in [%s, %s])'%(min_cluster_size, max_cluster_size),
                                                                                                            'chose from among %d heavy clusters with sizes %s' % (len(cpath.best()), cluster_size_str(cpath.best())),
                                                                                                            ))
        return ([sid, pid], [heavy_locus(ig_or_tr), plocus]), len(cluster)
    else:
        if debug:
            print('    chose seed uid %s from cluster with size %d%s%s%s (%s)' % (sid, len(cluster), ' (chosen at random)' if choose_random else '',
                                                                             '' if nth_seed==0 else ', after skipping %d/%d other cluster%s'%(nth_seed, len(cpath.best()), plural(nth_seed)),
                                                                                  '' if min_cluster_size is None and max_cluster_size is None else ' (was asked for size in [%s, %s])'%(min_cluster_size, max_cluster_size),
                                                                                  'chose from among %d clusters with sizes %s' % (len(cpath.best()), cluster_size_str(cpath.best())),
                                                                                  ))
        return sid, len(cluster)

# ----------------------------------------------------------------------------------------
# Takes two logd values and adds them together, i.e. takes (log a, log b) --> log a+b
# i.e. a *or* b
def add_in_log_space(first, second):
    if first == -float('inf'):
        return second
    elif second == -float('inf'):
        return first
    elif first > second:
        return first + math.log(1 + math.exp(second - first))
    else:
        return second + math.log(1 + math.exp(first - second))

# ----------------------------------------------------------------------------------------
def non_none(vlist):  # return the first non-None value in vlist (there are many, many places where i could go back and use this) [this avoids hard-to-read if/else statements that require writing the first val twice]
    for val in vlist:
        if val is not None:
            return val
    raise Exception('utils.non_none() called with all-None vlist: %s' % vlist)

# ----------------------------------------------------------------------------------------
def arglist_imatches(clist, argstr):
    assert argstr[:2] == '--'  # this is necessary since the matching assumes that argparse has ok'd the uniqueness of an abbreviated argument UPDATE now we've disable argparse prefix matching, but whatever
    return [i for i, c in enumerate(clist) if argstr==c]  # NOTE do *not* try to go back to matching just the start of the argument, in order to make that work you'd need to have access to the whole list of potential arguments in bin/partis, and you'd probably screw it up anyway

# ----------------------------------------------------------------------------------------
# NOTE this is not necessary any more since arglist_imatches() now only look for exact matches (since I've disable argparse prefix matchin), but I'm in the middle of too many things to remove it atm
def reduce_imatches(imatches, clist, argstr):  # restrict <imatches> to exact matches in an effort to get it down to one unique match
    imatches = [i for i in imatches if clist[i]==argstr]  # see if any of them are exact matches
    if len(imatches) > 1:
        raise Exception('multiple matches for argstr \'%s\' in cmd (this may be caused by not typing out the entirety of an arg string): %s' % (argstr, ' '.join(clist[i] for i in imatches)))
    return imatches

# ----------------------------------------------------------------------------------------
def arglist_index(clist, argstr):
    imatches = arglist_imatches(clist, argstr)
    if len(imatches) == 0:
        raise Exception('\'%s\' not found in cmd: %s' % (argstr, ' '.join(clist)))
    if len(imatches) > 1:
        imatches = reduce_imatches(imatches, clist, argstr)
    return get_single_entry(imatches)

# ----------------------------------------------------------------------------------------
def is_in_arglist(clist, argstr):  # accounts for argparse unique/partial matches (update: we no longer allow partial matches)
    return len(arglist_imatches(clist, argstr)) > 0

# ----------------------------------------------------------------------------------------
def has_arg_val(clist, argstr):  # return true if <argstr> has an accompanying value (i.e. is not a boolean)
    iarg = arglist_index(clist, argstr)
    return iarg < len(clist) - 1 and clist[iarg+1].find('--') != 0

# ----------------------------------------------------------------------------------------
def get_val_from_arglist(clist, argstr):
    imatch = arglist_index(clist, argstr)
    if imatch == len(clist) - 1:
        raise Exception('no argument for %s in %s' % (argstr, clist))
    val = clist[imatch + 1]
    if val[:2] == '--':
        raise Exception('no value for %s in %s (next word is %s)' % (argstr, clist, val))
    return val

# ----------------------------------------------------------------------------------------
def remove_from_arglist(clist, argstr, has_arg=False):
    if clist.count(None) > 0:
        raise Exception('None type value in clist %s' % clist)
    imatches = arglist_imatches(clist, argstr)
    if len(imatches) == 0:
        return clist
    if len(imatches) > 1:
        imatches = reduce_imatches(imatches, clist, argstr)
    iloc = imatches[0]
    # if clist[iloc] != argstr:
    #     print '  %s removing abbreviation \'%s\' from sys.argv rather than \'%s\'' % (color('yellow', 'warning'), clist[iloc], argstr)
    if has_arg:
        clist.pop(iloc + 1)
    clist.pop(iloc)
    return clist  # NOTE usually we don't use the return value (just modify it in memory), but at least in one place we're calling with a just-created list so it's nice to be able to use the return value

# ----------------------------------------------------------------------------------------
# replace the argument to <argstr> in <clist> with <replace_with>, or if <argstr> isn't there add it. If we need to add it and <insert_after> is set, add it after <insert_after>
def replace_in_arglist(clist, argstr, replace_with, insert_after=None, has_arg=False):
    if clist.count(None) > 0:
        raise Exception('None type value in clist %s' % clist)
    if not is_in_arglist(clist, argstr):
        if insert_after is None or insert_after not in clist:  # just append it
            clist.append(argstr)
            clist.append(replace_with)
        else:  # insert after the arg <insert_after>
            insert_in_arglist(clist, [argstr, replace_with], insert_after, has_arg=has_arg)
    else:
        clist[arglist_index(clist, argstr) + 1] = replace_with

# ----------------------------------------------------------------------------------------
def replace_in_argstr(cmd, argstr, replace_with, insert_after=None, has_arg=False):  # same as previous fcn, but on arg str (probably a bunch of places I could use this, but I'm adding it late)
    clist = cmd.split()
    replace_in_arglist(clist, argstr, replace_with, insert_after=insert_after, has_arg=has_arg)
    return ' '.join(clist)

# ----------------------------------------------------------------------------------------
# insert list <new_arg_strs> after <argstr> (unless <before> is set),  Use <has_arg> if <argstr> has an argument after which the insertion should occur
def insert_in_arglist(clist, new_arg_strs, argstr, has_arg=False, before=False):  # set <argstr> to None to put it at end (yeah it probably should've been a kwarg)
    i_insert = len(clist)
    if argstr is not None:
        i_insert = clist.index(argstr) + (0 if before else (2 if has_arg else 1))
    clist[i_insert : i_insert] = new_arg_strs

# ----------------------------------------------------------------------------------------
# return new arg list with args + values from both clist_1 and clist_2. Any args that are in both with different values will cause a crash
def merge_arg_lists(clist_1, clist_2, debug=False):
    if debug:
        print('    clist_1: %s' % ' '.join(clist_1))
        print('    clist_2: %s' % ' '.join(clist_2))
    argnames = [v for v in clist_1 if v.find('--')==0] + [v for v in clist_2 if v.find('--')==0 and v not in clist_1]  # would use a set, but then the order gets messed up
    new_clist = []
    for astr in argnames:
        alist = [astr]
        tclist = clist_1 if astr in clist_1 else clist_2
        if has_arg_val(tclist, astr):
            if astr in clist_1 and astr in clist_2:
                assert has_arg_val(clist_1, astr) and has_arg_val(clist_2, astr)
                if get_val_from_arglist(clist_1, astr) != get_val_from_arglist(clist_2, astr):  # if this arg has a value, make sure it's the same in both 
                    raise Exception('arg \'%s\' has different values in the two arg lists \'%s\' vs \'%s\'' % (astr, get_val_from_arglist(clist_1, astr), get_val_from_arglist(clist_2, astr)))
            alist.append(get_val_from_arglist(tclist, astr))
        new_clist += alist
    if debug:
        print('  new_clist: %s' % ' '.join(new_clist))
    return new_clist

# ----------------------------------------------------------------------------------------
def kbound_str(kbounds):
    return_str = []
    for region in ['v', 'd']:
        rkb = kbounds[region]
        return_str.append('k_%s %d' % (region, rkb['best']))
        if 'min' in rkb and 'max' in rkb:
            return_str.append(' [%s-%s)' % (str(rkb.get('min', ' ')), str(rkb.get('max', ' '))))
        return_str.append('  ')
    return ''.join(return_str).strip()

# ----------------------------------------------------------------------------------------
def jsdump(fname, jdata):
    with open(fname, 'w') as jfile:
        jfo = json.dumps(jdata)
        if sys.version_info.major < 3:
            jfo = jfo.decode()
        jfile.write(jfo)

# ----------------------------------------------------------------------------------------
def write_seqfos(fname, seqfos):  # NOTE basically just a copy of write_fasta(), except this writes to .yaml, and includes any extra info (beyond name and seq)
    mkdir(fname, isfile=True)
    jsdump(fname, seqfos)

# ----------------------------------------------------------------------------------------
def write_fasta(fname, seqfos, name_key='name', seq_key='seq'):  # should have written this a while ago -- there's tons of places where I could use this instead of writing it by hand, but I'm not going to hunt them all down now
    mkdir(fname, isfile=True)
    with open(fname, 'w') as seqfile:
        for sfo in seqfos:
            seqfile.write('>%s\n%s\n' % (sfo[name_key], sfo[seq_key]))

# ----------------------------------------------------------------------------------------
# NOTE replacement for (some cases of) read_fastx()
def read_seqfos(fname):  # queries=None, n_max_queries=-1, istartstop=None, ftype=None, n_random_queries=None): NOTE maybe add these args?
    if getsuffix(fname) not in ['.json', '.yaml']:
        raise Exception('unhandled suffix %s (should be .json or .yaml)' % getsuffix(fname))
    with open(fname) as sfile:
        try:  # ok this is kind of dumb but it's nice to be able to edit them by hand as yaml files
            seqfos = json.load(sfile)
        except ValueError:
            sfile.seek(0)
            seqfos = yaml.load(sfile, Loader=Loader)
    if 'germline-info' in seqfos:
        raise Exception('this is a standard yaml output file (not just a list of seq infos), so needs to be read with utils.read_yaml_output(): %s' % fname)
    return seqfos

# ----------------------------------------------------------------------------------------
# if <look_for_tuples> is set, look for uids that are actually string-converted python tuples, and add each entry in the tuple as a duplicate sequence. Can also pass in a list <tuple_info> if you need to do more with the info afterwards (this is to handle gctree writing fasta files with broken names; see usage also in datascripts/meta/taraki-XXX)
def read_fastx(fname, name_key='name', seq_key='seq', add_info=True, dont_split_infostrs=False, sanitize_uids=False, sanitize_seqs=False, queries=None, n_max_queries=-1, istartstop=None, ftype=None, n_random_queries=None, look_for_tuples=False, tuple_info=None, quiet=False):
    if ftype is None:
        suffix = getsuffix(fname)
        if suffix == '.fa' or suffix == '.fasta':
            ftype = 'fa'
        elif suffix == '.fq' or suffix == '.fastq':
            ftype = 'fq'
        else:
            raise Exception('unhandled file type: %s' % suffix)

    finfo = []
    iline = -1  # index of the query/seq that we're currently reading in the fasta
    n_fasta_queries = 0  # number of queries so far added to <finfo> (I guess I could just use len(finfo) at this point)
    missing_queries = set(queries) if queries is not None else None
    already_printed_forbidden_character_warning, already_printed_warn_char_warning = False, False
    with open(fname) as fastafile:
        startpos = None
        while True:
            if startpos is not None:  # rewind since the last time through we had to look to see when the next header line appeared
                fastafile.seek(startpos)
            headline = fastafile.readline()
            if not headline:
                break
            if headline.strip() == '':  # skip a blank line
                headline = fastafile.readline()

            if ftype == 'fa':
                if headline[0] != '>':
                    raise Exception('invalid fasta header line in %s:\n    %s' % (fname, headline))
                headline = headline.lstrip('>')

                seqlines = []
                nextline = fastafile.readline()
                while True:
                    if not nextline:
                        break
                    if nextline[0] == '>':
                        break
                    else:
                        startpos = fastafile.tell()  # i.e. very line that doesn't begin with '>' increments <startpos>
                    seqlines.append(nextline)
                    nextline = fastafile.readline()
                seqline = ''.join([l.strip() for l in seqlines]) if len(seqlines) > 0 else None
            elif ftype == 'fq':
                if headline[0] != '@':
                    raise Exception('invalid fastq header line in %s:\n    %s' % (fname, headline))
                headline = headline.lstrip('@')

                seqline = fastafile.readline()  # NOTE .fq with multi-line entries isn't supported, since delimiter characters are allowed to occur within the quality string
                plusline = fastafile.readline().strip()
                if plusline[0] != '+':
                    raise Exception('invalid fastq quality header in %s:\n    %s' % (fname, plusline))
                qualityline = fastafile.readline()
            else:
                raise Exception('unhandled ftype %s' % ftype)

            if not seqline:
                break

            iline += 1
            if istartstop is not None:
                if iline < istartstop[0]:
                    continue
                elif iline >= istartstop[1]:
                    continue

            if dont_split_infostrs:  # if this is set, we let the calling fcn handle all the infostr parsing (e.g. for imgt germline fasta files)
                infostrs = headline
                uid = infostrs
            else:  # but by default, we split by everything that could be a separator, which isn't really ideal, but we're reading way too many different kinds of fasta files at this point to change the default
                # NOTE commenting this since it fucking breaks on the imgt fastas, which use a different fucking format and i don't even remember whose stupid format this was for WTF PEOPLE
                # if ';' in headline and '=' in headline:  # HOLY SHIT PEOPLE DON"T PUT YOUR META INFO IN YOUR FASTA FILES
                #     infostrs = [s1.split('=') for s1 in headline.strip().split(';')]
                #     uid = infostrs[0][0]
                #     infostrs = dict(s for s in infostrs if len(s) == 2)
                # else:
                infostrs = [s3.strip() for s1 in headline.split(' ') for s2 in s1.split('\t') for s3 in s2.split('|')]  # NOTE the uid is left untranslated in here
                uid = infostrs[0]
            if sanitize_uids and any(fc in uid for fc in forbidden_characters):
                if not already_printed_forbidden_character_warning:
                    print('  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in forbidden_characters]), uid, uid.translate(forbidden_character_translations)))
                    already_printed_forbidden_character_warning = True
                uid = uid.translate(forbidden_character_translations)
            if sanitize_uids and any(wc in uid for wc in warn_chars):
                if not already_printed_warn_char_warning:
                    print('  %s: found a character that may cause problems if doing phylogenetic inference (one of %s) in sequence id \'%s\' (only printing this warning on first occurence).' % (color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in warn_chars]), uid))
                    already_printed_warn_char_warning = True

            if queries is not None:
                if uid not in queries:
                    continue
                missing_queries.remove(uid)

            seqfo = {name_key : uid, seq_key : seqline.strip().upper()}
            if add_info:
                seqfo['infostrs'] = infostrs
            if sanitize_seqs:
                seqfo[seq_key] = seqfo[seq_key].translate(ambig_translations).upper()
                if any(c not in alphabet for c in seqfo[seq_key]):
                    unexpected_chars = set([ch for ch in seqfo[seq_key] if ch not in alphabet])
                    raise Exception('unexpected character%s %s (not among %s) in input sequence with id %s:\n  %s' % (plural(len(unexpected_chars)), ', '.join([('\'%s\'' % ch) for ch in unexpected_chars]), alphabet, seqfo[name_key], seqfo[seq_key]))
            finfo.append(seqfo)

            n_fasta_queries += 1
            if n_max_queries > 0 and n_fasta_queries >= n_max_queries:
                break
            if queries is not None and len(missing_queries) == 0:
                break
    if n_max_queries > 0 and not quiet:
        print('    stopped after reading %d sequences from %s' % (n_max_queries, fname))
    if queries is not None:
        print('    only looked for %d specified sequences in %s' % (len(queries), fname))

    if n_random_queries is not None:
        if n_random_queries > len(finfo):
            print('  %s asked for n_random_queries %d from file with only %d entries, so just taking all of them (%s)' % (color('yellow', 'warning'), n_random_queries, len(finfo), fname))
            n_random_queries = len(finfo)
        finfo = numpy.random.choice(finfo, n_random_queries, replace=False)

    if look_for_tuples:  # this is because gctree writes broken fasta files with multiple uids in a line
        new_sfos, n_found = [], 0
        for sfo in finfo:
            headline = ' '.join(sfo['infostrs'])
            if any(c not in headline for c in '()'):
                new_sfos.append(sfo)
                continue
            istart, istop = [headline.find(c) for c in '()']
            try:
                nids = eval(headline[istart : istop + 1])
                if tuple_info is not None:
                    tuple_info.append(nids)
                n_found += 1
                print('         look_for_tuples: found %d uids in headline %s: %s' % (len(nids), headline, ' '.join(nids)))
                for uid in nids:
                    nsfo = copy.deepcopy(sfo)
                    nsfo['name'] = uid
                    new_sfos.append(nsfo)
            except:
                print('    %s failed parsing tuple from fasta line \'%s\' from %s' % (wrnstr(), headline, fname))
                new_sfos.append(sfo)
        if n_found > 0:
            print('      found %d seqfos with tuple headers (added %d net total seqs) in %s' % (n_found, len(new_sfos) - len(finfo), fname))
        finfo = new_sfos

    return finfo

# ----------------------------------------------------------------------------------------
def output_exists(args, outfname, outlabel=None, leave_zero_len=False, offset=None, todostr=None, debug=True):
    outlabel = '' if outlabel is None else ('%s ' % outlabel)
    if offset is None: offset = 22  # weird default setting method so we can call it also with the fcn below (without setting default value in two places)

    if not os.path.exists(outfname):
        return False

    if not leave_zero_len and os.stat(outfname).st_size == 0:
        if debug:
            print('%sdeleting zero length %s' % (offset * ' ', outfname))
        os.rmdir(outfname) if os.path.isdir(outfname) else os.remove(outfname)  # zero len dir means it's empty
        return False
    elif args.overwrite:
        if debug:
            print('%soverwriting %s%s' % (offset * ' ', outlabel, outfname))
        if os.path.isdir(outfname):
            raise Exception('output %s is a directory, rm it by hand' % outfname)
        else:
            os.remove(outfname)
        return False
    else:
        if debug:
            print('%s%soutput exists, %s (%s)' % (offset * ' ', outlabel, 'skipping' if todostr is None else todostr, outfname))
        return True

# ----------------------------------------------------------------------------------------
def all_outputs_exist(args, outfnames, outlabel=None, leave_zero_len=False, offset=None, dbgpathstr=None, debug=True):
    o_exist_list = [output_exists(args, ofn, outlabel=outlabel, leave_zero_len=leave_zero_len, offset=offset, debug=debug) for ofn in outfnames]
    n_exist = o_exist_list.count(True)
    if debug:
        ostr = '' if outlabel is None else ('%s ' % outlabel)
        print('    %s%d/%d outputs exist%s' % (ostr, n_exist, len(o_exist_list), '' if dbgpathstr is None else ' in %s'%dbgpathstr))
    return n_exist == len(o_exist_list)

# NOTE/reminder probably doesn't make sense to have this fcn, since output_exists() doesn't just test for existence, it also removes stuff, looks for zero length, etc.
# # ----------------------------------------------------------------------------------------
# def any_output_exists(args, outfnames, outlabel=None, offset=None, debug=True):
#     return any(output_exists(args, ofn, outlabel=outlabel, offset=offset, debug=debug) for ofn in outfnames)

# ----------------------------------------------------------------------------------------
# true if heavy and at least one light exist (note that this will still give the wrong answer if e.g. parameter writing fails between writing igk and igl, but i think it's better than the alternatives)
def lpair_outputs_exist(args, outfcn, outlabel=None, leave_zero_len=False, offset=None, dbgpathstr=None, todostr=None, ig_or_tr='ig', debug=True):
    existing_loci = []
    for lpr in locus_pairs[ig_or_tr]:  # note that we're use a locus pair here, but then only looking at the merged/non-pair parent dir (not the locus pair subdir)
        if all_outputs_exist(args, [outfcn(l) for l in lpr], outlabel=outlabel, leave_zero_len=leave_zero_len, offset=offset, dbgpathstr=dbgpathstr, debug=False):
            existing_loci += [l for l in lpr if l not in existing_loci]  # could use a set, but then i lose the order, so whatever
    if len(existing_loci) > 0:
        mstr = ''
        if len(existing_loci) < 3:
            missing_loci = [l for lp in locus_pairs[ig_or_tr] for l in lp if l not in existing_loci]
            mstr = ', %s %d %s: %s' % (color('yellow', 'missing'), len(missing_loci), 'locus' if len(missing_loci)==1 else 'loci', ' '.join(missing_loci))
        if debug:
            print('    %d / 3 locus outputs exist (%s%s), %s (e.g. %s)' % (len(existing_loci), ' '.join(existing_loci), mstr, 'skipping' if todostr is None else todostr, outfcn(existing_loci[0])))
        else:
            print('    %s' % mstr)
        return True
    return False

# ----------------------------------------------------------------------------------------
def getprefix(fname):  # basename before the dot
    if len(os.path.splitext(fname)) != 2:
        raise Exception('couldn\'t split %s into two pieces using dot' % fname)
    return os.path.splitext(fname)[0]

# ----------------------------------------------------------------------------------------
def getsuffix(fname):  # suffix, including the dot
    if len(os.path.splitext(fname)) != 2:
        raise Exception('couldn\'t split %s into two pieces using dot' % fname)
    return os.path.splitext(fname)[1]

# ----------------------------------------------------------------------------------------
def replace_suffix(fname, new_suffix):
    return fname.replace(getsuffix(fname), new_suffix)

# ----------------------------------------------------------------------------------------
def insert_before_suffix(insert_str, fname):
    if getsuffix(fname) == '':
        raise Exception('no suffix in %s' % fname)
    return fname.replace(getsuffix(fname), '%s%s' % (insert_str, getsuffix(fname)))

# ----------------------------------------------------------------------------------------
def read_vsearch_cluster_file(fname):
    id_clusters = {}
    with open(fname) as clusterfile:
        reader = csv.DictReader(clusterfile, fieldnames=['type', 'cluster_id', '3', '4', '5', '6', '7', 'crap', 'query', 'morecrap'], delimiter=str('\t'))
        for line in reader:
            if line['type'] == 'C':  # some lines are a cluster, and some are a query sequence. Skip the cluster ones.
                continue
            cluster_id = int(line['cluster_id'])
            if cluster_id not in id_clusters:
                id_clusters[cluster_id] = []
            uid = line['query']
            id_clusters[cluster_id].append(uid)
    partition = list(id_clusters.values())
    return partition

# ----------------------------------------------------------------------------------------
def read_vsearch_search_file(fname, userfields, seqdict, glfo, region, get_annotations=False, debug=False):
    from . import indelutils
    def get_mutation_info(query, matchfo, indelfo):
        tmpgl = glfo['seqs'][region][matchfo['gene']][matchfo['glbounds'][0] : matchfo['glbounds'][1]]
        if indelutils.has_indels(indelfo):
            tmpqr = indelfo['reversed_seq'][matchfo['qrbounds'][0] : matchfo['qrbounds'][1] - indelutils.net_length(indelfo)]
        else:
            tmpqr = seqdict[query][matchfo['qrbounds'][0] : matchfo['qrbounds'][1]]
        # color_mutants(tmpgl, tmpqr, print_result=True, align=True)
        return hamming_distance(tmpgl, tmpqr, return_len_excluding_ambig=True, return_mutated_positions=True)

    # first we add every match (i.e. gene) for each query
    query_info = {}
    with open(fname) as alnfile:
        reader = csv.DictReader(alnfile, fieldnames=userfields, delimiter=str('\t'))  # NOTE start/end positions are 1-indexed
        for line in reader:  # NOTE similarity to waterer.read_query()
            if line['query'] not in query_info:  # note that a surprisingly large number of targets give the same score, and you seem to get a little closer to what sw does if you sort alphabetically, but in the end it doesn't/shouldn't matter
                query_info[line['query']] = []
            query_info[line['query']].append({
                'ids' : int(line['ids']),
                'gene' : line['target'],
                'cigar' : line['caln'],
                'qrbounds' : (int(line['qilo']) - 1, int(line['qihi'])),
                'glbounds' : (int(line['tilo']) - 1, int(line['tihi'])),
            })

    # then we throw out all the matches (genes) that have id/score lower than the best one
    failed_queries = list(set(seqdict) - set(query_info))
    for query in list(query_info):
        if len(query_info[query]) == 0:
            print('%s zero vsearch matches for query %s' % (color('yellow', 'warning'), query))
            del query_info[query]  # uh... need to handle failures better than this
            failed_queries.append(query)
            continue
        query_info[query] = sorted(query_info[query], key=lambda d: d['ids'], reverse=True)  # sort the list of matches by decreasing score
        best_score = query_info[query][0]['ids']
        query_info[query] = [qinfo for qinfo in query_info[query] if qinfo['ids'] == best_score]  # keep all the matches that have the same score as the best match

    # then count up how many matches there were for each gene over all the sequences (to account for multiple matches with the same score, the count is not an integer)
    gene_counts = {}
    for query in query_info:
        counts_per_match = 1. / len(query_info[query])  # e.g. if there's four matches with the same score, give 'em each 0.25 counts
        for qinfo in query_info[query]:
            if qinfo['gene'] not in gene_counts:
                gene_counts[qinfo['gene']] = 0.
            gene_counts[qinfo['gene']] += counts_per_match

    annotations = OrderedDict()  # NOTE this is *not* a complete vdj annotation, it's just the info we have for an alignment to one region (presumably v)
    imatch = 0  # they all have the same score at this point, so just take the first one
    if get_annotations:  # it probably wouldn't really be much slower to just always do this
        if debug:
            print('      gene     score  mfreq  indels')
        for query in query_info:
            matchfo = query_info[query][imatch]
            v_indelfo = indelutils.get_indelfo_from_cigar(matchfo['cigar'], seqdict[query], matchfo['qrbounds'], glfo['seqs'][region][matchfo['gene']], matchfo['glbounds'], {'v' : matchfo['gene']}, vsearch_conventions=True, uid=query)
            n_mutations, len_excluding_ambig, mutated_positions = get_mutation_info(query, matchfo, v_indelfo)  # not sure this needs to just be the v_indelfo, but I'm leaving it that way for now
            combined_indelfo = indelutils.get_empty_indel()
            if indelutils.has_indels(v_indelfo):
                # huh, it actually works fine with zero-length d and j matches, so I don't need this any more
                # # arbitrarily give the d one base, and the j the rest of the sequence (I think they shouldn't affect anything as long as there's no d or j indels here)
                # if matchfo['qrbounds'][1] <= len(seqdict[query]) - 2:  # if we've got room to give 1 base to d and 1 base to j
                #     d_start = matchfo['qrbounds'][1]
                # else:  # ...but if we don't, then truncate the v match
                #     d_start = len(seqdict[query]) - 2
                # j_start = d_start + 1
                # tmpbounds = {
                #     'v' : (matchfo['qrbounds'][0], d_start),
                #     'd' : (d_start, j_start),
                #     'j' : (j_start, len(seqdict[query])),
                # }
                tmpbounds = {  # add zero-length d and j matches
                    'v' : matchfo['qrbounds'],
                    'd' : (matchfo['qrbounds'][1], matchfo['qrbounds'][1]),
                    'j' : (matchfo['qrbounds'][1], matchfo['qrbounds'][1]),
                }
                combined_indelfo = indelutils.combine_indels({'v' : v_indelfo}, seqdict[query], tmpbounds, debug=debug)
            annotations[query] = {
                region + '_gene' : matchfo['gene'],  # only works for v now, though
                'score' : matchfo['ids'],
                'n_' + region + '_mutations' : n_mutations,
                region + '_mut_freq' : float(n_mutations) / len_excluding_ambig,
                region + '_mutated_positions' : mutated_positions,
                'qrbounds' : {region : matchfo['qrbounds']},
                'glbounds' : {region : matchfo['glbounds']},
                'indelfo' : combined_indelfo,
            }
            if debug:
                atn = annotations[query]
                print('    %s  %3d   %5.2f  %s  %s' % (color_gene(matchfo['gene'], width=10), matchfo['ids'], atn[region+'_mut_freq'], color('red', 'indel') if indelutils.has_indels(atn['indelfo']) else '     ', query))

    return {'gene-counts' : gene_counts, 'annotations' : annotations, 'failures' : failed_queries}

# ----------------------------------------------------------------------------------------
# ok this sucks, but the original function below is used in a bunch of places that pass in a dict of seqs, and i don't want to mess with them since they're complicated, but now i also want it to be able to run on inputs with duplicate sequence ids
def run_vsearch_with_duplicate_uids(action, seqlist, workdir, threshold, **kwargs):
    # NOTE this is weird but necessary basically because in the normal partis workflows input always goes through seqfileopener, which renames duplicates, but now i want to be able to run this vsearch stuff without that (for bin/split-loci.py), which has .fa input files that i want to allow to have duplicates
    def get_trid(uid, seq):
        return '%s-DUP-%s' % (uid, uidhashstr(seq))
    returnfo = run_vsearch(action, {get_trid(s['name'], s['seq']) : s['seq'] for s in seqlist}, workdir, threshold, **kwargs)
    if set(returnfo) != set(['gene-counts', 'annotations', 'failures']):
        raise Exception('needs to be updated')
    annotation_list = []
    for sfo in seqlist:
        line = returnfo['annotations'].get(get_trid(sfo['name'], sfo['seq']), {'invalid' : True})  # these aren't really annotations, they just have a few random keys. I should've called them something else
        assert 'unique_ids' not in line  # make sure it wasn't added, since if it was we'd need to change it
        line['unique_ids'] = [sfo['name']]
        annotation_list.append(line)
    return annotation_list  # eh, probably don't need this:, returnfo['gene-counts']  # NOTE both annotation_list and <failures> can have duplicate uids in them

# ----------------------------------------------------------------------------------------
def get_platform_binstr():
    if platform.system() == 'Linux':
        return 'linux'
    elif platform.system() == 'Darwin':
        return 'macos'
    else:
        raise Exception('%s no vsearch binary in bin/ for platform \'%s\' (you can specify your own full vsearch path with --vsearch-binary)' % (color('red', 'error'), platform.system()))

# ----------------------------------------------------------------------------------------
# NOTE use the previous fcn if you expect duplicate uids
def run_vsearch(action, seqdict, workdir, threshold, match_mismatch='2:-4', gap_open=None, no_indels=False, minseqlength=None, consensus_fname=None, msa_fname=None, glfo=None, print_time=False, vsearch_binary=None, get_annotations=False, expect_failure=False, extra_str='  vsearch:', debug=False):
    from . import clusterpath
    from . import glutils
    # note: '2:-4' is the default vsearch match:mismatch, but I'm setting it here in case vsearch changes it in the future
    # single-pass, greedy, star-clustering algorithm with
    #  - add the target to the cluster if the pairwise identity with the centroid is higher than global threshold <--id>
    #  - pairwise identity definition <--iddef> defaults to: number of (matching columns) / (alignment length - terminal gaps)
    #  - the search process sorts sequences in decreasing order of number of k-mers in common
    #    - the search process stops after --maxaccept matches (default 1), and gives up after --maxreject non-matches (default 32)
    #    - If both are zero, it searches the whole database
    #    - I do not remember why I set both to zero. I just did a quick test, and on a few thousand sequences, it seems to be somewhat faster with the defaults, and a tiny bit less accurate.
    #      - UPDATE i seem to have gone back to default
    region = 'v'
    userfields = [  # all 1-indexed (note: only used for 'search')
        'query',
        'target',
        'qilo',  # first pos of query that aligns to target (1-indexed), skipping initial gaps (e.g. 1 if first pos aligns, 4 if fourth pos aligns but first three don't)
        'qihi',  # last pos of same (1-indexed)
        'tilo',  # same, but pos of target that aligns to query (1-indexed)
        'tihi',  # last pos of same (1-indexed)
        'ids',
        'caln',  # cigar string
    ]
    expected_success_fraction = 0.75  # if we get alignments for fewer than this, print a warning cause something's probably wrong

    start = time.time()
    prep_dir(workdir)
    infname = workdir + '/input.fa'

    # write input
    with open(infname, 'w') as fastafile:
        for name, seq in seqdict.items():
            fastafile.write('>' + name + '\n' + seq + '\n')

    # figure out which vsearch binary to use
    if vsearch_binary is None:
        vsearch_binary = '%s/bin/vsearch-2.4.3-%s-x86_64' % (get_partis_dir(), get_platform_binstr())

    # build command
    cmd = vsearch_binary
    cmd += ' --id ' + str(1. - threshold)  # reject if identity lower than this
    match, mismatch = [int(m) for m in match_mismatch.split(':')]
    assert mismatch < 0  # if you give it a positive one it doesn't complain, so presumably it's actually using that positive  (at least for v identification it only makes a small difference, but vsearch's default is negative)
    cmd += ' --match %d'  % match  # default 2
    cmd += ' --mismatch %d' % mismatch  # default -4
    # cmd += ' --gapext %dI/%dE' % (2, 1)  # default: (2 internal)/(1 terminal)
    # it would be nice to clean this up
    if gap_open is None:
        gap_open = 1000 if no_indels else 50
    cmd += ' --gapopen %dI/%dE' % (gap_open, 2)  # default: (20 internal)/(2 terminal)
    if minseqlength is not None:
        cmd += ' --minseqlength %d' % minseqlength
    if action == 'cluster':
        outfname = workdir + '/vsearch-clusters.txt'
        cmd += ' --cluster_fast ' + infname
        cmd += ' --uc ' + outfname
        if consensus_fname is not None:  # workdir cleanup below will fail if you put it in this workdir
            cmd += ' --consout ' + consensus_fname  # note: can also output a file with msa and consensus
        if msa_fname is not None:  # workdir cleanup below will fail if you put it in this workdir
            cmd += ' --msaout ' + msa_fname
        # cmd += ' --maxaccept 0 --maxreject 0'  # see note above
    elif action == 'search':
        outfname = workdir + '/aln-info.tsv'
        dbdir = workdir + '/' + glutils.glfo_dir
        glutils.write_glfo(dbdir, glfo)
        cmd += ' --usearch_global ' + infname
        cmd += ' --maxaccepts 5'  # it's sorted by number of k-mers in common, so this needs to be large enough that we'll almost definitely get the exact best gene match
        cmd += ' --db ' + glutils.get_fname(dbdir, glfo['locus'], region)
        cmd += ' --userfields %s --userout %s' % ('+'.join(userfields), outfname)  # note that --sizeout --dbmatched <fname> adds up all the matches from --maxaccepts, i.e. it's not what we want
    else:
        assert False
    # cmd += ' --threads ' + str(n_procs)  # a.t.m. just let vsearch use all the cores (it'd be nice to be able to control it a little, but it's hard to keep it separate from the number of slurm procs, and we pretty much always want it to be different to that)
    cmd += ' --quiet'
    cmdfos = [{
        'cmd_str' : cmd,
        'outfname' : outfname,
        'workdir' : workdir,
        # 'threads' : n_procs},  # NOTE that this does something very different (adjusts slurm command) to the line above ^ (which talks to vsearch)
    }]

    # run
    run_cmds(cmdfos)

    # read output
    n_seqs = len(seqdict)
    if action == 'cluster':
        returnfo = read_vsearch_cluster_file(outfname)
    elif action == 'search':
        returnfo = read_vsearch_search_file(outfname, userfields, seqdict, glfo, region, get_annotations=get_annotations, debug=debug)
        glutils.remove_glfo_files(dbdir, glfo['locus'])
        succ_frac = sum(returnfo['gene-counts'].values()) / float(n_seqs)
        if succ_frac < expected_success_fraction and not expect_failure:
            print('%s vsearch only managed to align %d / %d = %.3f of the input sequences (cmd below)   %s\n  %s' % (color('yellow', 'warning'), sum(returnfo['gene-counts'].values()), n_seqs, sum(returnfo['gene-counts'].values()) / float(n_seqs), reverse_complement_warning(), cmd))
    else:
        assert False
    os.remove(infname)
    os.remove(outfname)
    os.rmdir(workdir)

    if print_time:
        if action == 'search':
            # NOTE you don't want to remove these failures, since sw is much smarter about alignment than vsearch, i.e. some failures here are actually ok
            n_passed = int(round(sum(returnfo['gene-counts'].values())))
            tw = str(len(str(n_seqs)))  # width of format str from number of digits in N seqs
            print(('%s %'+tw+'d / %-'+tw+'d %s annotations (%d failed) with %d %s gene%s in %.1f sec') % (extra_str, n_passed, n_seqs, region, n_seqs - n_passed, len(returnfo['gene-counts']),
                                                                                                          region, plural(len(returnfo['gene-counts'])), time.time() - start))
        else:
            print('can\'t yet print time for clustering')

    return returnfo

# ----------------------------------------------------------------------------------------
def run_swarm(seqs, workdir, differences=1, n_procs=1):
    # groups together all sequence pairs that have <d> or fewer differences (--differences, default 1)
    #  - if d=1, uses algorithm of linear complexity (d=2 or greater uses quadratic algorithm)
    #  - --fastidious (only for d=1) extra pass to reduce the number of small OTUs

    prep_dir(workdir)
    infname = workdir + '/input.fa'
    outfname = workdir + '/clusters.txt'

    dummy_abundance = 1
    with open(infname, 'w') as fastafile:
        for name, seq in seqs.items():
            fastafile.write('>%s_%d\n%s\n' % (name, dummy_abundance, remove_ambiguous_ends(seq).replace('N', 'A')))

    cmd = '%s/bin/swarm-2.1.13-linux-x86_64 %s' % (get_partis_dir(), infname)
    cmd += ' --differences ' + str(differences)
    if differences == 1:
        cmd += ' --fastidious'
    cmd += ' --threads ' + str(n_procs)
    cmd += ' --output-file ' + outfname
    simplerun(cmd)

    partition = []
    with open(outfname) as outfile:
        for line in outfile:
            partition.append([uidstr.rstrip('_%s' % dummy_abundance) for uidstr in line.strip().split()])

    os.remove(infname)
    os.remove(outfname)
    os.rmdir(workdir)

    cp = clusterpath.ClusterPath()
    cp.add_partition(partition, logprob=0., n_procs=1)
    cp.print_partitions(abbreviate=True)

    return partition


# ----------------------------------------------------------------------------------------
def compare_vsearch_to_sw(sw_info, vs_info):
    from . import indelutils
    # pad vsearch indel info so it'll match the sw indel info (if the sw indel info is just copied from the vsearch info, and you're not going to use the vsearch info for anything after, there's no reason to do this)
    for query in sw_info['indels']:
        if query not in sw_info['queries']:
            continue
        if query not in vs_info['annotations']:
            continue
        if not indelutils.has_indels(vs_info['annotations'][query]['indelfo']):
            continue
        indelutils.pad_indelfo(vs_info['annotations'][query]['indelfo'], ambig_base * sw_info[query]['padlefts'][0], ambig_base * sw_info[query]['padrights'][0])

    from .hist import Hist
    hists = {
        # 'mfreq' : Hist(30, -0.1, 0.1),
        # 'n_muted' : Hist(20, -10, 10),
        # 'vs' : Hist(30, 0., 0.4),
        # 'sw' : Hist(30, 0., 0.4),
        'n_indels' : Hist(7, -3.5, 3.5),
        'pos' : Hist(21, -10.5, 10.5),
        'len' : Hist(11, -5, 5),
        # 'type' : Hist(2, -0.5, 1.5),
    }

    for uid in sw_info['queries']:
        if uid not in vs_info['annotations']:
            continue
        vsfo = vs_info['annotations'][uid]
        swfo = sw_info[uid]
        # swmfreq, swmutations = utils.get_mutation_rate_and_n_muted(swfo, iseq=0, restrict_to_region='v')
        # hists['mfreq'].fill(vsfo['v_mut_freq'] - swmfreq)
        # hists['n_muted'].fill(vsfo['n_v_mutations'] - swmutations)
        # hists['vs'].fill(vsfo['v_mut_freq'])
        # hists['sw'].fill(swmfreq)
        iindel = 0
        vsindels = vsfo['indelfo']['indels']
        swindels = swfo['indelfos'][0]['indels']
        hists['n_indels'].fill(len(vsindels) - len(swindels))
        if len(vsindels) == 0 or len(swindels) == 0:
            continue
        vsindels = sorted(vsindels, key=lambda d: (d['len'], d['pos']), reverse=True)
        swindels = sorted(swindels, key=lambda d: (d['len'], d['pos']), reverse=True)
        for name in [n for n in hists if n != 'n_indels']:
            hists[name].fill(vsindels[iindel][name] - swindels[iindel][name])

    from . import plotting
    for name, hist in hists.items():
        fig, ax = plotting.mpl_init()
        hist.mpl_plot(ax, hist, label=name, color=plotting.default_colors[list(hists.keys()).index(name)])
        # hist.write('_output/vsearch-test/%s/%s.csv' % (match_mismatch.replace(':', '-'), name))
        plotting.mpl_finish(ax, '.', name, xlabel='vs - sw')

# ----------------------------------------------------------------------------------------
def get_chimera_max_abs_diff(line, iseq, chunk_len=75, max_ambig_frac=0.1, debug=False):
    naive_seq, mature_seq = subset_iseq(line, iseq, restrict_to_region='v')  # self.info[uid]['naive_seq'], self.info[uid]['seqs'][0]

    if ambig_frac(naive_seq) > max_ambig_frac or ambig_frac(mature_seq) > max_ambig_frac:
        if debug:
            print('  too much ambig %.2f %.2f' % (ambig_frac(naive_seq), ambig_frac(mature_seq)))
        return None, 0.

    if debug:
        color_mutants(naive_seq, mature_seq, print_result=True)
        print(' '.join(['%3d' % s for s in isnps]))

    _, isnps = hamming_distance(naive_seq, mature_seq, return_mutated_positions=True)

    max_abs_diff, imax = 0., None
    for ipos in range(chunk_len, len(mature_seq) - chunk_len):
        if debug:
            print(ipos)

        muts_before = [isn for isn in isnps if isn >= ipos - chunk_len and isn < ipos]
        muts_after = [isn for isn in isnps if isn >= ipos and isn < ipos + chunk_len]
        mfreq_before = len(muts_before) / float(chunk_len)
        mfreq_after = len(muts_after) / float(chunk_len)

        if debug:
            print('    len(%s) / %d = %.3f' % (muts_before, chunk_len, mfreq_before))
            print('    len(%s) / %d = %.3f' % (muts_after, chunk_len, mfreq_after))

        abs_diff = abs(mfreq_before - mfreq_after)
        if imax is None or abs_diff > max_abs_diff:
            max_abs_diff = abs_diff
            imax = ipos

    return imax, max_abs_diff  # <imax> is break point

# ----------------------------------------------------------------------------------------
def get_version_info(debug=False):
    git_dir = '%s/.git' % get_partis_dir()
    if not os.path.exists(git_dir):
        try:
            from importlib.metadata import version
            return {'tag': version('partis-bcr'), 'commit': 'installed'}
        except ImportError:
            return {'tag': 'unknown', 'commit': 'installed'}
    vinfo = {}
    vinfo['commit'] = subprocess.check_output(['git', '--git-dir', git_dir, 'rev-parse', 'HEAD'], universal_newlines=True).strip()
    if debug:
        print('  commit: %s' % vinfo['commit'])
    cmd = 'git --git-dir %s describe --always --tags' % git_dir
    out, err = simplerun(cmd, return_out_err=True, debug=False)
    if '-' in out:
        if out.count('-') == 2:
            vinfo['tag'], vinfo['n_ahead_of_tag'], commit_hash_abbrev = out.strip().split('-')
            if debug:
                ahead_str = ''
                if int(vinfo['n_ahead_of_tag']) > 0:
                    ahead_str = '  (well, %d commits ahead of)' % int(vinfo['n_ahead_of_tag'])
                print('     tag: %s%s' % (vinfo['tag'], ahead_str))
        else:
            vinfo['tag'] = '?'
            print('    %s utils.get_version_info() couldn\'t figure out tag from \'%s\' output: %s' % (color('red', 'error'), cmd, out))
    else:
        vinfo['tag'] = out.strip()
        print('     tag: %s' % vinfo['tag'])

    return vinfo

# ----------------------------------------------------------------------------------------
def write_only_partition(fname, partition):
    from . import clusterpath
    cpath = clusterpath.ClusterPath(partition=partition)
    write_annotations(fname, None, [], None, partition_lines=cpath.get_partition_lines())

# ----------------------------------------------------------------------------------------
def write_empty_annotations(fname, locus):
    from . import glutils
    write_annotations(fname, glutils.get_empty_glfo(locus), [], annotation_headers)

# ----------------------------------------------------------------------------------------
def write_annotations(fname, glfo, annotation_list, headers, synth_single_seqs=False, failed_queries=None, partition_lines=None, use_pyyaml=False, dont_write_git_info=False):
    from . import clusterpath
    if os.path.exists(fname):
        os.remove(fname)
    elif not os.path.exists(os.path.dirname(os.path.abspath(fname))):
        os.makedirs(os.path.dirname(os.path.abspath(fname)))

    if getsuffix(fname) == '.csv':
        assert partition_lines is None
        write_csv_annotations(fname, headers, annotation_list, synth_single_seqs=synth_single_seqs, glfo=glfo, failed_queries=failed_queries)
    elif getsuffix(fname) == '.yaml':
        if partition_lines is None:
            partition_lines = clusterpath.ClusterPath(partition=get_partition_from_annotation_list(annotation_list)).get_partition_lines()
        write_yaml_output(fname, headers, glfo=glfo, annotation_list=annotation_list, synth_single_seqs=synth_single_seqs, failed_queries=failed_queries, partition_lines=partition_lines, use_pyyaml=use_pyyaml, dont_write_git_info=dont_write_git_info)
    else:
        raise Exception('unhandled file extension \'%s\' on %s' % (getsuffix(fname), fname))

# ----------------------------------------------------------------------------------------
def write_csv_annotations(fname, headers, annotation_list, synth_single_seqs=False, glfo=None, failed_queries=None):
    with open(fname, csv_wmode()) as csvfile:
        writer = csv.DictWriter(csvfile, headers)
        writer.writeheader()
        for line in annotation_list:
            if synth_single_seqs:
                for iseq in range(len(line['unique_ids'])):
                    outline = get_line_for_output(headers, synthesize_single_seq_line(line, iseq), glfo=glfo)
                    writer.writerow(outline)
            else:
                outline = get_line_for_output(headers, line, glfo=glfo)
                writer.writerow(outline)
        if failed_queries is not None:
            for failfo in failed_queries:
                assert len(failfo['unique_ids']) == 1
                writer.writerow({'unique_ids' : failfo['unique_ids'][0], 'invalid' : failfo['invalid'], 'input_seqs' : failfo['input_seqs'][0]})  # this is ugly, but the corresponding one in the yaml fcn is nice

# ----------------------------------------------------------------------------------------
def get_yamlfo_for_output(line, headers, glfo=None):
    yamlfo = {}
    if not line['invalid']:
        transfer_indel_info(line, yamlfo)
    for key in [k for k in headers if k not in yamlfo]:
        if key in line:
            yamlfo[key] = line[key]
        elif key in list(input_metafile_keys.values()):  # these are optional, so if they're missing, don't worry about it
            continue
        else:
            if line['invalid']:
                continue
            add_extra_column(key, line, yamlfo, glfo=glfo)
    return yamlfo

# ----------------------------------------------------------------------------------------
def write_yaml_output(fname, headers, glfo=None, annotation_list=None, synth_single_seqs=False, failed_queries=None, partition_lines=None, use_pyyaml=False, dont_write_git_info=False):
    # ----------------------------------------------------------------------------------------
    def check_ids():  # really just want to check that there's *some* overlap between the partitions and annotations. It's normal that there's annotations for only some uids in the partition, but if there is a partition, at least some of its uids should be in the annotations (usually the uids with annotations is a strict subset, but sometimes there might be annotations for uids from somewhere else)
        ptnids = set(u for p in partition_lines for c in p['partition'] for u in c)
        antnids = set(u for l in annotation_list for u in l['unique_ids'])  # note: doesn't include duplicates, but that's fine
        if len(ptnids) > 0 and len(antnids) > 0 and len(ptnids & antnids) == 0:
            print('  %s writing partitions (%d uids) and annotations (%d uids) with zero overlap to %s' % (color('yellow', 'warning'), len(ptnids), len(antnids), fname))
    # ----------------------------------------------------------------------------------------
    if annotation_list is None:
        annotation_list = []
    if partition_lines is None:
        partition_lines = []
    check_ids()

    version_info = {'partis-yaml' : 0.1, 'partis-git' : '' if dont_write_git_info else get_version_info()}
    yaml_annotations = [get_yamlfo_for_output(l, headers, glfo=glfo) for l in annotation_list]
    if failed_queries is not None:
        yaml_annotations += failed_queries
    yamldata = {'version-info' : version_info,
                'germline-info' : glfo,
                'partitions' : partition_lines,
                'events' : yaml_annotations}
    if use_pyyaml:  # slower, but easier to read by hand for debugging (use this instead of the json version to make more human-readable files)
        with open(fname, 'w') as yamlfile:
            yaml.dump(yamldata, yamlfile, width=400, Dumper=Dumper, default_flow_style=False, allow_unicode=False)  # set <allow_unicode> to false so the file isn't cluttered up with !!python.unicode stuff
    else:  # way tf faster than full yaml (only lost information is ordering in ordered dicts, but that's only per-gene support and germline info, neither of whose order we care much about)
        jsdump(fname, yamldata) #, sort_keys=True, indent=4)

# ----------------------------------------------------------------------------------------
def parse_yaml_annotations(glfo, yamlfo, n_max_queries, synth_single_seqs, dont_add_implicit_info):
    annotation_list = []
    n_queries_read = 0
    for line in yamlfo['events']:
        if not line['invalid']:
            transfer_indel_reversed_seqs(line)
            if 'all_matches' in line and isinstance(line['all_matches'], dict):  # it used to be per-family, but then I realized it should be per-sequence, so any old cache files lying around have it as per-family
                line['all_matches'] = [line['all_matches']]  # also, yes, it makes me VERY ANGRY that this needs to be here, but i just ran into a couple of these old files and otherwise they cause crashes
            if not dont_add_implicit_info:  # it's kind of slow, although most of the time you probably want all the extra info
                add_implicit_info(glfo, line)  # don't use the germline info in <yamlfo>, in case we decide we want to modify it in the calling fcn
        if synth_single_seqs and len(line['unique_ids']) > 1:
            for iseq in range(len(line['unique_ids'])):
                annotation_list.append(synthesize_single_seq_line(line, iseq))
        else:
            annotation_list.append(line)

        n_queries_read += len(line['unique_ids'])
        if n_max_queries > 0 and n_queries_read >= n_max_queries:
            break

    return annotation_list

# ----------------------------------------------------------------------------------------
def read_cpath(fname, n_max_queries=-1, seed_unique_id=None, skip_annotations=False):
    _, _, cpath = read_output(fname, n_max_queries=n_max_queries, seed_unique_id=seed_unique_id, dont_add_implicit_info=True, skip_annotations=skip_annotations) # can't skip annotations by default since older (simulation, i think) files will need them to make partitions
    return cpath

# ----------------------------------------------------------------------------------------
def read_output(fname, n_max_queries=-1, synth_single_seqs=False, dont_add_implicit_info=False, seed_unique_id=None, cpath=None, skip_annotations=False, glfo=None, glfo_dir=None, locus=None, skip_failed_queries=False, is_partition_file=False, debug=False):
    from . import clusterpath
    from . import glutils
    annotation_list = None

    if getsuffix(fname) == '.csv':
        cluster_annotation_fname = fname.replace('.csv', '-cluster-annotations.csv')
        if os.path.exists(cluster_annotation_fname) or is_partition_file:  # i.e. if <fname> is a partition file
            assert cpath is None   # see note in read_yaml_output()
            cpath = clusterpath.ClusterPath(fname=fname, seed_unique_id=seed_unique_id)  # NOTE I'm not sure if I really want to pass in the seed here -- it should be stored in the file -- but if it's in both places it should be the same. um, should.
            fname = cluster_annotation_fname  # kind of hackey, but oh well

        if not skip_annotations:
            if not dont_add_implicit_info and glfo is None:
                if glfo_dir is not None:
                    glfo = glutils.read_glfo(glfo_dir, locus)
                else:
                    raise Exception('glfo is None, but we were asked to add implicit info for an (old-style) csv output file')
            n_queries_read = 0
            annotation_list = []
            with open(fname) as csvfile:
                for line in csv.DictReader(csvfile):
                    process_input_line(line, skip_literal_eval=dont_add_implicit_info)  # NOTE kind of weird to equate implicit info adding and literal eval skipping... but in the end they're both mostly speed optimizations
                    if not dont_add_implicit_info:
                        add_implicit_info(glfo, line)
                    annotation_list.append(line)
                    n_queries_read += 1
                    if n_max_queries > 0 and n_queries_read >= n_max_queries:
                        break

    elif getsuffix(fname) == '.yaml':  # NOTE this replaces any <glfo> that was passed (well, only within the local name table of this fcn, unless the calling fcn replaces it themselves, since we return this glfo)
        glfo, annotation_list, cpath = read_yaml_output(fname, n_max_queries=n_max_queries, synth_single_seqs=synth_single_seqs,
                                                        dont_add_implicit_info=dont_add_implicit_info, seed_unique_id=seed_unique_id, cpath=cpath, skip_annotations=skip_annotations, debug=debug)
    else:
        raise Exception('unhandled file extension \'%s\' on %s' % (getsuffix(fname), fname))

    if cpath is not None and len(cpath.partitions) == 0 and not skip_annotations:  # old simulation files didn't write the partition separately, but we may as well get it
        cpath.add_partition(get_partition_from_annotation_list(annotation_list), -1., 1)

    if skip_failed_queries:
        failed_queries = set([u for l in annotation_list for u in l['unique_ids'] if l['invalid']])
        annotation_list = [l for l in annotation_list if not l['invalid']]
        removed_queries = set()
        for ip, partition in enumerate(cpath.partitions):
            for cluster in partition:
                for fq in set(cluster) & failed_queries:
                    cluster.remove(fq)
                    removed_queries.add(fq)
            cpath.partitions[ip] = [c for c in partition if len(c) > 0]  # remove any empty ones (typically we'll be removing single-sequence clusters, since they failed)
        if len(removed_queries) > 0:
            print('      removed %d failed queries when reading partition: %s' % (len(removed_queries), ' '.join(sorted(removed_queries))))

    if annotation_list is not None:
        for antn in annotation_list:
            add_per_seq_keys(antn)

    return glfo, annotation_list, cpath  # NOTE if you want a dict of annotations, use utils.get_annotation_dict() above

# ----------------------------------------------------------------------------------------
# NOTE if there's multiple files with conflicting info this will *not* detect it (unlike in seqfileopener, which will)
# NOTE also this'll only merge correctly input-metafo style yaml/json files, i.e. that're a simple dict at top level
def read_json_yamls(fnames):  # try to read <fname> as json (since it's faster), on exception fall back to yaml
    # ----------------------------------------------------------------------------------------
    def merge_yfos(tfo, yfo):
        for uid in tfo:
            if uid in yfo:
                yfo[uid].update(tfo[uid])
            else:
                yfo[uid] = tfo[uid]
    # ----------------------------------------------------------------------------------------
    yamlfo = None
    for fn in fnames:
        tmpfo = read_json_yaml(fn)
        if yamlfo is None:
            yamlfo = tmpfo
        else:
            merge_yfos(tmpfo, yamlfo)
    return yamlfo

# ----------------------------------------------------------------------------------------
def read_json_yaml(fname):  # try to read <fname> as json (since it's faster), on exception fall back to yaml
    with open(fname) as yamlfile:
        try:
            yamlfo = json.load(yamlfile)  # way tf faster than full yaml (only lost information is ordering in ordered dicts, but that's only per-gene support and germline info, neither of whose order we care much about)
        except ValueError:  # I wish i could think of a better way to do this, but I can't
            yamlfile.seek(0)
            yamlfo = yaml.load(yamlfile, Loader=Loader)  # use this instead of the json version to make more human-readable files
    return yamlfo

# ----------------------------------------------------------------------------------------
def read_yaml_output(fname, n_max_queries=-1, synth_single_seqs=False, dont_add_implicit_info=False, seed_unique_id=None, cpath=None, skip_annotations=False, debug=False):
    from . import clusterpath
    yamlfo = read_json_yaml(fname)
    if isinstance(yamlfo, list):
        raise Exception('read list of seqfos from file, instead of the expected standard yaml output with germline-info, annotations, and partitions. Run read_seqfos() instead: %s' % fname)
    if debug:
        print('  read yaml version %s from %s' % (yamlfo['version-info']['partis-yaml'], fname))

    glfo = yamlfo['germline-info']  # it would probably be good to run the glfo through the checks that glutils.read_glfo() does, but on the other hand since we're reading from our own yaml file, those have almost certainly already been done

    annotation_list = None
    if not skip_annotations:  # may not really be worthwhile, but oh well
        annotation_list = parse_yaml_annotations(glfo, yamlfo, n_max_queries, synth_single_seqs, dont_add_implicit_info)

    partition_lines = yamlfo['partitions']
    if cpath is None:   # allowing the caller to pass in <cpath> is kind of awkward, but it's used for backward compatibility in clusterpath.readfile()
        cpath = clusterpath.ClusterPath(seed_unique_id=seed_unique_id)  # NOTE I'm not sure if I really want to pass in the seed here -- it should be stored in the file -- but if it's in both places it should be the same. um, should.
    if len(partition_lines) > 0:  # *don't* combine this with the cluster path constructor, since then we won't modify the path passed in the arguments
        cpath.readlines(partition_lines)

    if annotation_list is not None:
        for antn in annotation_list:
            add_per_seq_keys(antn)

    return glfo, annotation_list, cpath  # NOTE if you want a dict of annotations, use utils.get_annotation_dict() above

# ----------------------------------------------------------------------------------------
def get_gene_counts_from_annotations(annotations, only_regions=None):
    gene_counts = {r : {} for r in (only_regions if only_regions is not None else regions)}
    for query, line in annotations.items():
        for tmpreg in gene_counts:
            gene = line[tmpreg + '_gene']
            if gene not in gene_counts[tmpreg]:
                gene_counts[tmpreg][gene] = 0.
            gene_counts[tmpreg][gene] += 1.  # vsearch info counts partial matches based of score, but I don't feel like dealing with that here at the moment
    return gene_counts

# ----------------------------------------------------------------------------------------
def parse_constant_regions(species, locus, annotation_list, workdir, aa_dbg=False, csv_outdir=None, debug=False):
    from . import glutils
    # ----------------------------------------------------------------------------------------
    def algncreg(tkey, n_min_seqs=5):
        if tkey == 'leader':  # for leaders, group together seqs that align to each V gene
            all_v_genes = set(l['v_gene'] for l in annotation_list)
            vg_antns = {}
            for vgene in all_v_genes:
                vg_antns[vgene] = [l for l in annotation_list if l['v_gene']==vgene]
            print('  found %d total V genes' % len(all_v_genes))
        else:  # For c_genes, do everyone at once
            vg_antns = {'IGHVx-x*x' : annotation_list}
        # leader_seq_infos = []  # final/new leader seqs
        writefos = []
        print('      gene    families  seqs')
        for ivg, (vgene, vgalist) in enumerate(sorted(list(vg_antns.items()), key=lambda q: sum(len(l['unique_ids']) for l in q[1]), reverse=True)):
            print('    %s %3d     %3d' % (color_gene(vgene, width=10), len(vgalist), sum(len(l['unique_ids']) for l in vgalist)))
            if ivg > 2 and sum(len(l['unique_ids']) for l in vgalist) < n_min_seqs:  # always do the first 3 v gene groups, but past that skip any vgene groups that have too few sequences
                print('            too few seqs')
                continue
            n_zero_len = len([s for l in vgalist for s in l['%s_seqs'%tkey] if len(s) == 0])
            if n_zero_len > 0:
                print('    ignoring %d seqs with zero length %ss' % (n_zero_len, tkey))
            qfos = [{'name' : u, 'seq' : s} for l in vgalist for u, s in zip(l['unique_ids'], l['%s_seqs'%tkey]) if len(s) > 0]  # it might be easier to do each annotation separately, but this way i can control parallelization better and there's less overhead
            if len(qfos) == 0:
                print('    no query infos')
                return
            matchfos, mstats = run_blastn(qfos, tgtfos[tkey], workdir, debug=debug) #, print_all_matches=True)  # , short=True
            if len(matchfos) == 0:
                print('    %s no %s matches from %d targets' % (wrnstr(), tkey, len(tgtfos[tkey])))
            qdict, tdict = [{s['name'] : s['seq'] for s in sfos} for sfos in [qfos, tgtfos[tkey]]]  # should really check for duplicates
            for itg, (tgt, mfos) in enumerate(mstats):  # loop over each target seq and the sequences for whom it was a best match
                print('      %s  %d' % (color('blue', tgt), len(mfos)))
                if itg > 2 and tkey=='c_gene' and len(mfos) < n_min_seqs:  # always print the first three target matches, but past that skip any targets that have to few sequences matched to them
                    print('          too small, skipping')
                    continue
                tmsfos = [{'name' : tgt, 'seq' : tdict[tgt]}] + [{'name' : m['query'], 'seq' : qdict[m['query']]} for m in mfos]
                msa_seqfos = align_many_seqs(tmsfos, no_gaps=True, extra_str='            ', debug=False)  # this really, really sucks to re-align, but it's only for debug and otherwise I'd have to write a parser for the stupid blast alignment output UPDATE well now i'm parsing the blast alignment a little more (although still not the gaps, whose locations may not be there?) so maybe coule remove this now
                cseq = cons_seq(aligned_seqfos=[s for s in msa_seqfos if s['name']!=tgt], extra_str='         ', debug=debug)
                # leader_seq_infos.append({'seq' : cseq, 'match-name' : tgt})
                writefos += [{'name' : s['name'], 'seq' : s['seq'], 'target' : tgt} for s in msa_seqfos + [{'name' : 'consensus-%s'%tgt, 'seq' : cseq, 'target' : tgt}]]
                if debug:
                    print(color_mutants(cseq, get_single_entry([s for s in msa_seqfos if s['name']==tgt])['seq'], seq_label='target: ', post_str=' %s'%tgt, extra_str='              '))  # can't put the target in the cons seq calculation, so have to print separately
                    if aa_dbg:
                        aasfos = [{'name' : s['name'], 'seq' : ltranslate(pad_nuc_seq(s['seq'].replace('-', ambig_base), side='left' if tkey=='leader' else 'right'))} for s in msa_seqfos]  # pad on left, since we assume the last three bases of 'leader_seqs' are in frame w/respect to V
                        aa_msa_seqfos = align_many_seqs(aasfos, no_gaps=True, extra_str='            ', aa=True, debug=False)
                        aa_cseq = cons_seq(aligned_seqfos=[s for s in aa_msa_seqfos if s['name']!=tgt], aa=True, extra_str='         ', debug=debug)
                        print(color_mutants(aa_cseq, get_single_entry([s for s in aa_msa_seqfos if s['name']==tgt])['seq'], amino_acid=True, seq_label='target: ', post_str=' %s'%tgt, extra_str='              '))  # can't put the target in the cons seq calculation, so have to print separately
        if csv_outdir is not None and len(writefos) > 0:
            cfn = '%s/%s-seq-alignments.csv'%(csv_outdir, tkey)
            print('      writing to %s' % cfn)
            with open(cfn, csv_wmode()) as cfile:
                writer = csv.DictWriter(cfile, fieldnames=list(writefos[0].keys()))
                writer.writeheader()
                for wfo in writefos:
                    writer.writerow(wfo)
    # ----------------------------------------------------------------------------------------
    def read_seqs(tkey, fname, tdbg=False):
        if tdbg:
            print('    reading %s seqs from %s' % (tkey, fname))
        tfos = read_fastx(fname)
        for tfo in tfos:
            if len(tfo['infostrs']) > 1:
                tfo['name'] = glutils.get_imgt_info(tfo['infostrs'], 'gene')
        for tk in ['name', 'seq']:  # eh, duplicate seqs are ok i think
            all_vals = [t[tk] for t in tfos]
            if any(all_vals.count(v) > 1 for v in set(all_vals)):
                dupfo = sorted([(v, all_vals.count(v)) for v in sorted(set(all_vals)) if all_vals.count(v) > 1], key=operator.itemgetter(1), reverse=True)
                print('    %s %d %s %ss appear more than once (%d total occurences) in %s' % (wrnstr(), len(dupfo), tkey, tk, sum(n for _, n in dupfo), fname))
                if tdbg:
                    print('      %s' % '\n      '.join('%3d  %s'%(n, v) for v, n in dupfo))
                if tk == 'name':
                    tkvals, new_tfos = [], []
                    for tfo in tfos:
                        if tfo[tk] in tkvals:
                            continue
                        new_tfos.append(tfo)
                        tkvals.append(tfo[tk])
                    print('      removed %d duplicate %ss' % (len(tfos) - len(new_tfos), tk))
                    tfos = new_tfos
        return tfos
    # ----------------------------------------------------------------------------------------
    lsrc = 'imgt'
    leaderfn, cgfn = 'data/germlines/leaders/%s/%s/%s/%sv.fa' % (lsrc, species, locus, locus), 'data/germlines/constant/%s/%sc.fa'%(species, locus)
    tgtfos = collections.OrderedDict()
    for tkey, fn in zip(['leader', 'c_gene'], [leaderfn, cgfn]):
        print(' %s seqs' % color('blue', tkey))
        tgtfos[tkey] = read_seqs(tkey, fn)
        algncreg(tkey)
