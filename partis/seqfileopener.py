from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import numpy
import bz2
import gzip
import copy
import os
import sys
import csv
from collections import OrderedDict
import random
import re

from . import utils
from io import open

delimit_info = {'.csv' : ',', '.tsv' : '\t'}

# ----------------------------------------------------------------------------------------
def read_input_metafo(input_metafnames, annotation_list, required_keys=None, n_warn_print=10, debug=False):
    for mfname in input_metafnames:
        read_single_input_metafo(mfname, annotation_list, required_keys=required_keys, n_warn_print=n_warn_print, debug=debug)

# ----------------------------------------------------------------------------------------
def read_single_input_metafo(input_metafname, annotation_list, required_keys=None, n_warn_print=10, debug=False):  # read input metafo from <input_metafname> and put in <annotation_list> (when we call this below, <annotation_list> is <input_info> (wait so at this point it sounds like this fcn and the next one should be merged [although that would be hard and dangerous, so i'm not doing it now)
    if not os.path.exists(input_metafname):
        print('  %s specified input meta file %s doesn\'t exist, so not reading anything' % (utils.wrnstr(), input_metafname))
        return
    # NOTE <annotation_list> doesn't need to be real annotations, it only uses the 'unique_ids' key
    metafo = utils.read_json_yaml(input_metafname)
    if any(isinstance(tkey, int) for tkey in metafo):  # would be better to check for not being a string, but that's harder, and this probably only happens for my simulation hash ids
        raise Exception('meta info keys need to be string (maybe just need to add \'\' around sequence ids in yaml file), but got: %s' % ' '.join(str(type(tk)) for tk in metafo if isinstance(tk, int)))
    metafile_keys = set(k for mfo in metafo.values() for k in mfo)
    if len(metafile_keys) == 0:  # zero length meta info file
        return
    utils.add_input_meta_keys(metafile_keys)
    if required_keys is not None and len(set(required_keys) - metafile_keys) > 0:
        raise Exception('required metafile key(s) (%s) not found in %s' % (', '.join(set(required_keys) - metafile_keys), input_metafname))

    usets = {c : {lk : set() for lk in utils.input_metafile_keys.values()} for c in ['no-info', 'added']}  # for counting each seq/uid
    llists = {c : {lk : [] for lk in utils.input_metafile_keys.values()} for c in ['no-info', 'added']}  # for counting each line
    n_modified, modified_keys = 0, set()
    for line in annotation_list:
        for input_key, line_key in utils.input_metafile_keys.items():
            mvals = [utils.input_metafile_defaults(line_key) for _ in line['unique_ids']]  # fill <mvals> with default values, then modify the values for seqs that have info
            n_seqs_with_info = 0  # have to keep track of this, rather than just counting the non-default-value ones (which we used to do), because sometimes the default is info in the file that we want to keep (e.g. for paired uids)
            for iseq, uid in enumerate(line['unique_ids']):
                if uid not in metafo or input_key not in metafo[uid]:  # skip this seq if it doesn't have info
                    usets['no-info'][line_key].add(uid)
                    continue
                mval = metafo[uid][input_key]
                if line_key in line and mval != line[line_key][iseq]:  # in general, the meta info shouldn't already be in the input file if you're also specifying a separate meta file
                    if n_modified < n_warn_print:
                        print('  %s replacing \'%s\'/\'%s\' value for \'%s\' with value from %s: %s --> %s' % (utils.color('yellow', 'warning'), input_key, line_key, uid, input_metafname, line[line_key][iseq], mval))
                    n_modified += 1
                    modified_keys.add(line_key)
                if input_key == 'multiplicity' and mval == 0:
                    raise Exception('input meta info value for \'multiplicity\' must be greater than 1 (since it includes this sequence), but got %d for \'%s\'' % (mval, uid))
                mvals[iseq] = mval
                usets['added'][line_key].add(uid)
                n_seqs_with_info += 1

            if n_seqs_with_info > 0:  # we used to add it even if they were all empty, but that means that you always get all the possible input meta keys, which is super messy (the downside of skipping them is some seqs can have them while others don't)
                line[line_key] = mvals
                llists['added'][line_key].append(line['unique_ids'])
            else:
                llists['no-info'][line_key].append(line['unique_ids'])

    if n_modified > 0:
        print('%s replaced input metafo for %d instances of key%s %s%s' % (utils.color('yellow', 'warning'), n_modified, utils.plural(modified_keys), ', '.join(modified_keys), (' (see above, only printed the first %d)'%n_warn_print) if n_modified > n_warn_print else ''))
    added_keys = set(k for k, us in usets['added'].items() if len(us) > 0)
    added_uids = set(u for us in usets['added'].values() for u in us)
    print('  --input-metafnames: %s %s%s' % ('no meta info added from' if len(added_uids)==0 else 'added meta info for %d sequences from'%len(added_uids), input_metafname, '' if len(added_keys)==0 else ': '+' '.join(sorted(added_keys))))
    if debug:
        print('  read_input_metafo(): add input metafo from meta file to annotations')
        print('                       uids      uids     lines  lines')
        print('                     no-info  with-info  no-info added')
        for ik, lk in utils.input_metafile_keys.items():
            if len(usets['added'][lk]) == 0:
                continue
            print('     %15s   %3d      %3d       %3d    %3d' % (lk, len(usets['no-info'][lk]), len(usets['added'][lk]), len(llists['no-info'][lk]), len(llists['added'][lk])))

# ----------------------------------------------------------------------------------------
def add_input_metafo(input_info, annotation_list, keys_not_to_overwrite=None, n_max_warn_print=10, overwrite_all=False, debug=False):  # transfer input metafo from <input_info> (i.e. what was in --input-metafnames) to <annotation_list>
    # NOTE this input meta info stuff is kind of nasty, just because there's so many ways that/steps at which we want to be able to specify it: from --input-metafnames, from <input_info>, from <sw_info>. If it's all consistent it's fine, and if it isn't consistent it'll print the warning, so should also be fine.
    # NOTE <keys_not_to_overwrite> should be the *line* key
    # overwrite_all: by default, for each meta key, this ignores <line>s in <annotation_list> that don't have any info for that meta key in <input_info>.
    #    If <overwrite_all> is set, this instead overwrites all <line>s for all meta keys, even if they have no info in <input_info>.
    #    Used e.g. when updating meta info in existing output files where the meta info has changed (if the meta info is consistent, there's no reason to look at keys/lines that have no info in <input_info>)
    # NOTE may need to add <overwrite_all> to some more calls (only added it to the one that I was sure of)
    # ----------------------------------------------------------------------------------------
    def incr(itp, line_key, uids):
        usets[itp][line_key]['seqs'] |= set(uids)  # it's a little weird to use a set, but the same uid is sometimes in different <line>s, and i guess maybe it'd be better to treat them separately, this way is probably simpler
        usets[itp][line_key]['clusters'] += 1
    # ----------------------------------------------------------------------------------------
    def inpfo(uid):  # ick, but sometimes a uid is missing from <input_info>, probably only possible when running 'update-meta-info' action
        return input_info.get(uid, {})
    # ----------------------------------------------------------------------------------------
    itypes = ['no-info', 'with-info', 'one-source', 'identical', 'not-overwritten', 'replaced']
    usets = {c : {lk : {'seqs' : set(), 'clusters' : 0} for lk in utils.input_metafile_keys.values()} for c in itypes}
    n_warn_printed = 0
    for line in annotation_list:
        uids = line['unique_ids']
        for line_key in sorted(utils.input_metafile_keys.values()):
            n_with_info = len([u for u in uids if line_key in inpfo(u)])  # NOTE there can be info in <line> but *not* in <input_info>, e.g. if we added multiplicity information in waterer.
            if n_with_info == 0:
                incr('no-info', line_key, uids)  # it's a little weird to use a set, but the same uid is sometimes in different <line>s, and i guess maybe it'd be better to treat them separately, this way is probably simpler
                if not overwrite_all:
                    continue
            incr('with-info', line_key, uids)  # it's a little weird to use a set, but the same uid is sometimes in different <line>s, and i guess maybe it'd be better to treat them separately, this way is probably simpler
            mvals = [inpfo(u)[line_key][0] if line_key in inpfo(u) else utils.input_metafile_defaults(line_key) for u in uids]

            if line_key not in line:  # info in <input_info>, but not in <line> so just copy it over
                line[line_key] = mvals
                incr('one-source', line_key, uids)
                continue

            if mvals == line[line_key]:  # both sources agree, so nothing to do
                incr('identical', line_key, uids)
                continue

            if keys_not_to_overwrite is not None and line_key in keys_not_to_overwrite:  # they're different, but we were explicitly told not to overwrite the info in <line> (e.g. in waterer we reset 'multiplicities' to account for duplicate sequences, then when this fcn gets called by partitiondriver we need to keep those reset values)
                if debug:
                    print('    not overwriting non-matching info for %s' % line_key)
                incr('not-overwritten', line_key, uids)
                continue

            # actually replace the info in <line>
            if n_warn_printed < n_max_warn_print:
                print('  %s replacing input metafo \'%s\' value for %d seq cluster \'%s\': %s --> %s' % (utils.color('yellow', 'warning'), line_key, len(uids), ':'.join(uids), line[line_key], mvals))
                n_warn_printed += 1
            line[line_key] = mvals
            incr('replaced', line_key, uids)

    if any(len(us['seqs']) > 0 for us in usets['replaced'].values()) > 0:
        modified_keys = [lk for lk, us in usets['replaced'].items() if len(us['seqs']) > 0]
        def tsum(itp, st): return sum((len if st=='seqs' else utils.pass_fcn)(us[st]) for us in usets[itp].values())
        print('%s replaced input metafo for %d total key%s in %d seqs (%d cluster%s): %s%s' % (utils.color('yellow', 'warning'), len(modified_keys), utils.plural(len(modified_keys)), tsum('replaced', 'seqs'), tsum('replaced', 'clusters'), utils.plural(tsum('replaced', 'clusters')),
                                                                                               ', '.join(modified_keys), (' (see above, only printed the first %d clusters)'%n_max_warn_print) if n_warn_printed >= n_max_warn_print else ''))
        debug = True
    if debug:
        print('  add_input_metafo(): transferring input metafo from <input_info> to <line> (counts in table are uids, not lines)')
        print('                            no-info  with-info  one-src ident. not-overwr. replaced')
        for ik, lk in utils.input_metafile_keys.items():
            if len(usets['with-info'][lk]['seqs']) == 0:
                continue
            tstr = ''.join('%9d'%len(usets[itp][lk]['seqs']) for itp in itypes)
            print('     %s%s' % (utils.wfmt(lk, max(len(k) for k in utils.input_metafile_keys.values())), tstr))

# ----------------------------------------------------------------------------------------
def post_process(input_info, reco_info, args, infname, found_seed, is_data, iline):
    if args is None:
        return

    if args.istartstop is not None:
        n_lines_in_file = iline + 1
        if n_lines_in_file < args.istartstop[1]:
            raise Exception('--istartstop upper bound %d larger than number of lines in file %d' % (args.istartstop[1], n_lines_in_file))
    if len(input_info) == 0:
        if args.queries is not None:
            raise Exception('didn\'t find the specified --queries (%s) in %s' % (str(args.queries), infname))
        if args.reco_ids is not None:
            raise Exception('didn\'t find the specified --reco-ids (%s) in %s' % (str(args.reco_ids), infname))
    if args.queries is not None:
        missing_queries = set(args.queries) - set(input_info)
        extra_queries = set(input_info) - set(args.queries)  # this is just checking for a bug in the code just above here...
        if len(missing_queries) > 0:  # used to crash on this, but with paired loci we'll never find the uids from the opposite chain
            print('  %s didn\'t find some of the specified --queries: %s' % (utils.wrnstr(), ' '.join(missing_queries)))
        if len(extra_queries) > 0:
            raise Exception('extracted uids %s that weren\'t specified with --queries' % ' '.join(extra_queries))
    if args.seed_unique_id is not None and not found_seed:
        raise Exception('couldn\'t find seed unique id %s in %s' % (args.seed_unique_id, infname))
    elif args.random_seed_seq:  # already checked (in bin/partis) that other seed args aren't set
        args.seed_unique_id = random.choice(list(input_info.keys()))
        print('    chose random seed unique id %s' % args.seed_unique_id)

    if args.n_random_queries is not None:
        included_queries = set()  # only for dbg printing
        uids_to_choose_from = list(input_info.keys())
        if args.seed_unique_id is not None:
            uids_to_choose_from.remove(args.seed_unique_id)
            included_queries.add(args.seed_unique_id)
        if args.queries_to_include is not None:
            for uid in args.queries_to_include:
                if args.seed_unique_id is not None and uid == args.seed_unique_id:
                    continue
                if uid not in uids_to_choose_from:
                    raise Exception('couldn\'t find requested query %s in %s' % (uid, infname))
                uids_to_choose_from.remove(uid)
                included_queries.add(uid)
        if args.n_random_queries >= len(input_info):
            print('  %s --n-random-queries %d >= number of queries read from %s (so just keeping everybody)' % (utils.color('yellow', 'warning'), args.n_random_queries, infname))
        else:
            uids_to_remove = numpy.random.choice(uids_to_choose_from, len(input_info) - args.n_random_queries, replace=False)
            for uid in uids_to_remove:
                del input_info[uid]
                if reco_info is not None:
                    del reco_info[uid]
            print('  --n-random-queries: keeping %d / %d sequences from input file (removed %d%s)' % (len(input_info), len(input_info) + len(uids_to_remove), len(uids_to_remove),
                                                                                                      (' and specifically kept %s' % ' '.join(included_queries)) if len(included_queries) > 0 else ''))

# ----------------------------------------------------------------------------------------
def get_seqfile_info(x, is_data=False):
    raise Exception('renamed and changed returned vals (see below)')

# ----------------------------------------------------------------------------------------
def read_sequence_file(infname, is_data, n_max_queries=-1, args=None, simglfo=None, quiet=False, more_input_info=None, dont_add_implicit_info=False):
    # NOTE renamed this from get_seqfile_info() since I'm changing the return values, but I don't want to update the calls everywhere (e.g. in compareutils)
    yaml_glfo = None
    suffix = utils.getsuffix(infname)
    if suffix in delimit_info:
        seqfile = open(infname)  # closes on function exit. no, this isn't the best way to do this
        reader = csv.DictReader(seqfile, delimiter=delimit_info[suffix])
        sanitize_uids = True
    elif suffix in ['.fa', '.fasta', '.fq', '.fastq', '.fastx']:
        add_info = args is not None and args.name_column is not None and 'fasta-info-index' in args.name_column
        reader = utils.read_fastx(infname, name_key='unique_ids', seq_key='input_seqs', add_info=add_info, sanitize_uids=True, n_max_queries=n_max_queries,  # NOTE don't use istarstop kw arg here, 'cause it fucks with the istartstop treatment in the loop below
                                  queries=(args.queries if (args is not None and not args.abbreviate) else None), sanitize_seqs=args.sanitize_input_seqs)  # NOTE also can't filter on args.queries here if we're also translating
        sanitize_uids = False  # already did it, so don't need to do it below (unless fasta-info-index is set, in which case the above will have sanitized the wrong thing...
    elif suffix == '.yaml':
        yaml_glfo, reader, _ = utils.read_yaml_output(infname, n_max_queries=n_max_queries, synth_single_seqs=True, dont_add_implicit_info=True)  # not really sure that long term I want to synthesize single seq lines, but for backwards compatibility it's nice a.t.m.
        if not is_data:
            simglfo = yaml_glfo  # doesn't replace the contents, of course, which is why we return it
        sanitize_uids = True
    else:
        raise Exception('unhandled file extension \'%s\' on file \'%s\'' % (suffix, infname))

    input_info = OrderedDict()
    reco_info = None
    if not is_data:
        reco_info = OrderedDict()
    n_duplicate_uids = 0
    printed_simu_mismatch_warning, already_printed_forbidden_character_warning = False, False
    n_queries_added = 0
    found_seed = False
    potential_names, used_names = None, None  # for abbreviating
    iname = None  # line number -- used as sequence id if there isn't a name column in the file
    iline = -1
    for line in reader:
        iline += 1
        if args is not None:
            if args.istartstop is not None:
                if iline < args.istartstop[0]:
                    continue
                if iline >= args.istartstop[1]:
                    break
            if args.name_column is not None:
                if 'infostrs' in line and args.name_column.split('-')[:3] == ['fasta', 'info', 'index']:
                    assert len(args.name_column.split('-')) == 4
                    line['unique_ids'] = line['infostrs'][int(args.name_column.split('-')[3])]
                else:
                    line['unique_ids'] = line[args.name_column]
                    del line[args.name_column]
            if args.seq_column is not None:
                line['input_seqs'] = line[args.seq_column]
                if args.seq_column != 'seqs':  # stupid god damn weird backwards compatibility edge case bullshit
                    del line[args.seq_column]
        if iname is None and 'unique_ids' not in line and 'unique_id' not in line:
            print('  %s: couldn\'t find a name (unique id) column, so using line number as the sequence label (you can set the name column with --name-column)' % (utils.color('yellow', 'warning')))
            iname = 0
        if iname is not None:
            line['unique_ids'] = '%09d' % iname
            iname += 1
        if 'input_seqs' not in line and 'seq' not in line:
            raise Exception('couldn\'t find a sequence column in %s (you can set this with --seq-column)' % infname)
        if suffix != '.yaml':
            utils.process_input_line(line)
        if len(line['unique_ids']) > 1:
            raise Exception('can\'t yet handle multi-seq csv input files')
        uid = line['unique_ids'][0]
        if sanitize_uids and any(fc in uid for fc in utils.forbidden_characters):
            if not already_printed_forbidden_character_warning:
                print('  %s: found a forbidden character (one of %s) in sequence id \'%s\'. This means we\'ll be replacing each of these forbidden characters with a single letter from their name (in this case %s). If this will cause problems you should replace the characters with something else beforehand.' % (utils.color('yellow', 'warning'), ' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), uid, uid.translate(utils.forbidden_character_translations)))
                already_printed_forbidden_character_warning = True
            uid = uid.translate(utils.forbidden_character_translations)
        if uid in input_info:
            # Crash on duplicates if requested (typically set by paired loci parent process for single-chain subprocesses)
            if args is not None and args.crash_on_duplicate_uids:
                raise Exception('Found duplicate UID \'%s\' in %s. Cannot handle duplicate UIDs since pairing info references the original names. Please remove duplicates from input files.' % (uid, infname))
            uid, n_duplicate_uids = utils.choose_non_dup_id(uid, n_duplicate_uids, input_info)
        inseq = line['input_seqs'][0]

        # # it would be nice to check here for forbidden characters (in addition to in the .fa code above), but it's hard because we won't have read the csv properly above if it has them
        # if any(fc in uid for fc in utils.forbidden_characters):
        #     raise Exception('found a forbidden character (one of %s) in sequence id \'%s\'' % (' '.join(["'" + fc + "'" for fc in utils.forbidden_characters]), uid))
        if args is not None:
            if args.abbreviate:  # note that this changes <uid>, but doesn't modify <line>
                uid, potential_names, used_names = utils.choose_new_uid(potential_names, used_names)
            if args.queries is not None and uid not in args.queries:
                continue
            if args.reco_ids is not None and line['reco_id'] not in args.reco_ids:
                continue
            if args.seed_unique_id is not None and uid == args.seed_unique_id:
                found_seed = True

        if uid in input_info:
            raise Exception('found uid \'%s\' twice in input file %s' % (uid, infname))

        if args.sanitize_input_seqs:
            inseq = inseq.translate(utils.ambig_translations).upper()
        if len(set(inseq) - set(utils.alphabet)) > 0:  # NOTE should really be integrated with sanitize_seqs arg in utils.read_fastx() UPDATE just added the clause above, which maybe is sufficient?
            unexpected_chars = set([ch for ch in inseq if ch not in utils.alphabet])
            raise Exception('unexpected character%s %s (not among %s) in input sequence with id %s (maybe should set --sanitize-input-seqs?):\n  %s' % (utils.plural(len(unexpected_chars)), ', '.join([('\'%s\'' % ch) for ch in unexpected_chars]), utils.alphabet, uid, inseq))

        # da business
        input_info[uid] = {'unique_ids' : [uid, ], 'seqs' : [inseq, ]}

        if not is_data:
            if 'v_gene' not in line and not dont_add_implicit_info:
                raise Exception('simulation info not found in %s' % infname)
            reco_info[uid] = line  # this used to be deepcopy'd, but it's really slow and i'm really pretty sure it's not necessary
            if uid != line['unique_ids'][0] and not printed_simu_mismatch_warning:
                print('     note: uid in simulation info %s doesn\'t match input file uid %s (latter was probably changed above). Simulation info will be internally consistent, but the key indexing that info in <reco_info> will be different, since it corresponds to the newly chosen uid above.' % (uid, line['unique_ids'][0]))
                printed_simu_mismatch_warning = True
            if simglfo is not None and not dont_add_implicit_info:
                utils.add_implicit_info(simglfo, reco_info[uid])
            for line_key in utils.input_metafile_keys.values():
                if line_key in reco_info[uid]:  # this is kind of weird to copy from sim info to input info, but it makes sense because affinity is really meta info (the only other place affinity could come from is --input-metafnames below). Where i'm defining meta info more or less as any input info besides name and sequence (i think the distinction is only really important because we want to support fastas, which can't [shouldn't!] handle anything else))
                    input_info[uid][line_key] = copy.deepcopy(reco_info[uid][line_key])  # note that the args.input_metafnames stuff below should print a warning if you've also specified that (which you shouldn't, if it's simulation)

        n_queries_added += 1
        if n_max_queries > 0 and n_queries_added >= n_max_queries:
            if not quiet:  # just adding <quiet>, and too lazy to decide what other print statements it should effect, this is the only one I care about right now
                print('  --n-max-queries: stopped after reading %d queries from input file' % len(input_info))
            break

    if n_duplicate_uids > 0:
        print('  %s renamed %d duplicate uids from %s' % (utils.color('yellow', 'warning'), n_duplicate_uids, infname))

    if more_input_info is not None:  # if you use this on simulation, the extra queries that aren't in <reco_info> may end up breaking something down the line (but I don't imagine this really getting used on simulation)
        if len(set(more_input_info) & set(input_info)) > 0:  # check for sequences in both places
            common_uids = set(more_input_info) & set(input_info)
            print('  note: found %d queries in both --infname and --queries-to-include-fname: %s' % (len(common_uids), ' '.join(common_uids)))  # not necessarily a problem, but you probably *shouldn't* have sequences floating around in two different files
            differing_seqs = [q for q in common_uids if more_input_info[q]['seqs'][0] != input_info[q]['seqs'][0]]
            if len(differing_seqs) > 0:  # if they have different sequences, though, that's a problem
                for q in differing_seqs:
                    print(q)
                    utils.color_mutants(input_info[q]['seqs'][0], more_input_info[q]['seqs'][0], align_if_necessary=True, print_result=True, ref_label='  --infname  ', seq_label='  --queries-to-include-fname  ')
                raise Exception('inconsistent sequences for %d of the queries in both --infname and --queries-to-include-fname (see preceding lines)' % len(differing_seqs))
        if args is not None and args.seed_unique_id is not None and args.seed_unique_id in more_input_info:
            found_seed = True
        input_info.update(more_input_info)
    if args is not None and args.input_metafnames is not None:
        read_input_metafo(args.input_metafnames, list(input_info.values()))
    post_process(input_info, reco_info, args, infname, found_seed, is_data, iline)

    if len(input_info) == 0:
        print('  %s didn\'t read any sequences from %s' % (utils.wrnstr(), infname))

    return input_info, reco_info, yaml_glfo
