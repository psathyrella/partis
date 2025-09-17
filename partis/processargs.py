from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import os
import random
import sys
import subprocess

from . import utils
from . import glutils
from . import treeutils

def get_dummy_outfname(workdir, locus=None):
    return '%s/XXX-dummy-simu%s.yaml' % (workdir, '-'+locus if locus is not None else '')

actions_not_requiring_input = ['simulate', 'view-output', 'merge-paired-partitions', 'view-annotations', 'view-partitions', 'view-cluster-annotations', 'plot-partitions', 'view-alternative-annotations', 'get-selection-metrics', 'get-linearham-info', 'write-fake-paired-annotations']

# ----------------------------------------------------------------------------------------
# split this out so we can call it from both bin/partis and bin/test-germline-inference.py
def process_gls_gen_args(args):  # well, also does stuff with non-gls-gen new allele args
    if args.locus is not None:  # if args.paired_loci is not set
        # ----------------------------------------------------------------------------------------
        def get_defaults(aname):
            return utils.get_arg_list(getattr(glutils, 'default_'+aname.replace('_sim', ''))[args.locus], key_list=utils.regions)  # note that this leaves them as strings
        # ----------------------------------------------------------------------------------------
        args.n_genes_per_region = utils.get_arg_list(args.n_genes_per_region, intify=True, key_list=utils.regions, default_vals=get_defaults('n_genes_per_region'))
        if args.n_genes_per_region is not None and not utils.has_d_gene(args.locus) and args.n_genes_per_region['d'] > 1:
            print('  --n-genes-per-region: reducing light chain d genes from %d to 1' % args.n_genes_per_region['d'])  # this is kind of hackey/annoying, but atm it's better than adding a separate arg for light chain genes
            args.n_genes_per_region['d'] = 1
        args.n_sim_alleles_per_gene = utils.get_arg_list(args.n_sim_alleles_per_gene, floatify=True, key_list=utils.regions, default_vals=get_defaults('n_sim_alleles_per_gene'))
    positions = {
        'snp' : utils.get_arg_list(args.snp_positions),
        'indel' : utils.get_arg_list(args.indel_positions),
    }
    numbers = {
        'snp' : utils.get_arg_list(args.nsnp_list, intify=True),
        'indel' : utils.get_arg_list(args.nindel_list, intify=True),
    }
    delattr(args, 'snp_positions')  # just to make sure you don't accidentally use them (should only use the new args.new_allele_info that gets created below)
    delattr(args, 'indel_positions')
    delattr(args, 'nsnp_list')
    delattr(args, 'nindel_list')

    n_new_alleles = None
    mtypes = ['snp', 'indel']
    for mtype in mtypes:
        if positions[mtype] is not None:  # if specific positions were specified on the command line
            positions[mtype] = [[int(p) for p in pos_str.split(',')] for pos_str in positions[mtype]]  # NOTE I think I could switch this to utils.get_arg_list() with list_of_lists=True
            if len(positions[mtype]) != len(args.sim_v_genes):  # we shouldn't be able to get here unless args has .sim_v_genes
                raise Exception('--%s-positions %s and --sim-v-genes %s not the same length (%d vs %d)' % (mtype, positions[mtype], args.sim_v_genes, len(positions[mtype]), len(args.sim_v_genes)))
        if numbers[mtype] is not None:
            if not args.generate_germline_set and len(numbers[mtype]) != len(args.sim_v_genes):  # we shouldn't be able to get here unless args has .sim_v_genes
                raise Exception('--n%s-list %s and --sim-v-genes %s not the same length (%d vs %d)' % (mtype, numbers[mtype], args.sim_v_genes, len(numbers[mtype]), len(args.sim_v_genes)))
            if positions[mtype] is not None:
                raise Exception('can\'t specify both --n%s-list and --%s-positions' % (mtype, mtype))
            positions[mtype] = [[None for _ in range(number)] for number in numbers[mtype]]  # the <None> tells glutils to choose a position at random
        if positions[mtype] is not None:
            if n_new_alleles is None:
                n_new_alleles = len(positions[mtype])
            if len(positions[mtype]) != n_new_alleles:
                raise Exception('mismatched number of new alleles for %s' % ' vs '.join(mtypes))
    if n_new_alleles is None:
        n_new_alleles = 0
    for mtype in mtypes:
        if positions[mtype] is None:  # if it wasn't specified at all, i.e. we don't want to generate any new alleles
            positions[mtype] = [[] for _ in range(n_new_alleles)]
    args.new_allele_info = [{'gene' : args.sim_v_genes[igene] if not args.generate_germline_set else None,  # we shouldn't be able to get here unless args has .sim_v_genes
                             'snp-positions' : positions['snp'][igene],
                             'indel-positions' : positions['indel'][igene]}
                            for igene in range(n_new_alleles)]

# ----------------------------------------------------------------------------------------
def get_workdir(batch_system):  # split this out so we can use it in datascripts (ok, then I ended up commenting it in datascripts, but maybe later I want to uncomment)
    basestr = os.getenv('USER', default='partis-work')
    if batch_system is not None and os.path.exists('/fh/fast/matsen_e'):
        workdir = utils.choose_random_subdir('/fh/fast/matsen_e/%s/_tmp/hmms' % basestr)
    else:
        workdir = utils.choose_random_subdir('/tmp/%s/hmms' % basestr)
        if batch_system is not None:
            print('  %s: using batch system %s with default --workdir (%s) -- if this dir isn\'t visible to your batch nodes, you\'ll need to set --workdir to something that is' % (utils.color('red', 'warning'), batch_system, workdir))
    return workdir

# ----------------------------------------------------------------------------------------
def process(args):
    if args.action == 'subset-partition':
        if args.infname is None and args.paired_indir is None:
            raise Exception('have to set either --infname or --paired-indir for \'subset-partition\'')
        if args.outfname is None and args.paired_outdir is None:
            raise Exception('have to set either --outfname or --paired-outdir for \'subset-partition\'')
        if not args.paired_loci:
            print('  note: turning on --paired-loci since \'subset-partition\' requires it (and turning on --keep-all-unpaired-seqs)')
            args.keep_all_unpaired_seqs = True
            sys.argv.append('--keep-all-unpaired-seqs')
            args.paired_loci = True
            sys.argv.append('--paired-loci')
            if args.paired_outdir is None:
                args.paired_outdir = utils.getprefix(args.outfname)
                args.outfname = None
                utils.remove_from_arglist(sys.argv, '--outfname', has_arg=True)
        if args.paired_outdir is None:
            raise Exception('have to set --paired-outdir (or --outfname) for \'subset-partition\'')
        if args.seed_unique_id is not None:  # note sure that it'd be much work, but there's much less reason to need it
            raise Exception('--seed-unique-id is not yet supported with \'subset-partition\'')

    if args.outfname is None and args.paired_outdir is None and args.action in utils.existing_output_actions:
        raise Exception('--outfname (or --paired-outdir, if using --paired-loci) required for %s' % args.action)

    if args.action == 'run-viterbi':
        print('  note: replacing deprecated action name \'run-viterbi\' with current name \'annotate\' (you don\'t need to change anything unless you want this warning message to go away)')
        args.action = 'annotate'
    if args.action == 'view-alternative-naive-seqs':
        print('  note: replacing deprecated action name \'view-alternative-naive-seqs\' with current name \'view-alternative-annotations\' (you don\'t need to change anything unless you want this warning message to go away)')
        args.action = 'view-alternative-annotations'
    if args.seed_seq is not None:
        raise Exception('--seed-seq is deprecated, use --seed-unique-id and --queries-to-include-fname')
    if args.seed is not None:
        print('  note: moving value from deprecated arg --seed to --random-seed (you don\'t need to change anything unless you want this warning message to go away)')
        args.random_seed = args.seed
        delattr(args, 'seed')
        assert '--seed' in sys.argv
        utils.replace_in_arglist(sys.argv, '--seed', str(args.random_seed))  # have to also modify argv in case we're using it e.g. in paired locus stuff
        sys.argv[utils.arglist_index(sys.argv, '--seed')] = '--random-seed'
    if args.input_metafname is not None:
        print('  note: moving value from deprecated arg --input-metafname to --input-metafnames [note plural] (you don\'t need to change anything unless you want this warning message to go away)')
        assert args.input_metafnames is None
        args.input_metafnames = args.input_metafname
        delattr(args, 'input_metafname')
        assert '--input-metafname' in sys.argv
        utils.replace_in_arglist(sys.argv, '--input-metafname', args.input_metafnames)  # have to also modify argv in case we're using it e.g. in paired locus stuff
        sys.argv[utils.arglist_index(sys.argv, '--input-metafname')] = '--input-metafnames'

    if args.action == 'merge-paired-partitions':
        assert args.paired_loci
    if args.paired_loci:
        args.locus = None
        if [args.infname, args.paired_indir].count(None) == 0:
            raise Exception('can\'t specify both --infname and --paired-indir')
        if args.paired_indir is not None:
            if args.guess_pairing_info:
                print('  %s --guess-pairing-info has no effect on input when --paired-indir is set (only when --infname is set)' % utils.wrnstr())
            for attrname in ['max', 'random']:
                if getattr(args, 'n_%s_queries'%attrname) != (-1 if attrname=='max' else None):  # ick ick ick i wish i could go back and make default for max be None
                    print('  %s --n-%s-queries can only try to choose properly paired seqs if --infname is set (since the args can then be passed to bin/split-loci.py); otherwise seqs are chosen from the different locus input files without regard to whether they\'re paired with each other.\n          But since --paired-indir is set we now have to assume the args were passed when that dir was generated (if they were, then all is well).' % (utils.wrnstr(), attrname))
        if args.outfname is not None:
            raise Exception('can\'t set --outfname if --paired-loci is set (use --paired-outdir)')
        if args.plotdir == 'paired-outdir':
            args.plotdir = args.paired_outdir
        if args.plotdir is None and args.action == 'plot-partitions':
            args.plotdir = args.paired_outdir
        if args.seed_unique_id is not None:
            args.seed_unique_id = utils.get_arg_list(args.seed_unique_id)
            args.seed_loci = utils.get_arg_list(args.seed_loci, choices=utils.loci)
            if len(args.seed_unique_id) != 2 or args.seed_loci is None or len(args.seed_loci) != 2:
                raise Exception('if --seed-unique-id and --paired-loci are set, both --seed-unique-id and --seed-loci must be set to colon-separated lists of length two')
            if utils.has_d_gene(args.seed_loci[1]) or not utils.has_d_gene(args.seed_loci[0]):
                raise Exception('--seed-loci has to have one heavy and one light locus, with the heavy one first (e.g. igh:igk) but got %s' % args.seed_loci)
        else:
            if args.seed_loci is not None:
                raise Exception('doesn\'t make sense to set --seed-loci without also setting --seed-unique-id')
        if args.persistent_cachefname is not None:
            assert args.persistent_cachefname == 'paired-outdir'

        args.light_chain_fractions = utils.get_arg_list(args.light_chain_fractions, key_val_pairs=True, floatify=True, choices=utils.light_loci(args.ig_or_tr))
        if args.light_chain_fractions is not None and not utils.is_normed(list(args.light_chain_fractions.values())):
            raise Exception('--light-chain-fractions %s don\'t add to 1: %f' % (args.light_chain_fractions, sum(args.light_chain_fractions.values())))
        if args.single_light_locus is not None:
            if args.single_light_locus not in utils.light_loci(args.ig_or_tr):
                raise Exception('--single-light-chain-locus must be among: %s (got \'%s\')' % (utils.light_loci(args.ig_or_tr), args.single_light_locus))
            for ltmp in utils.light_loci(args.ig_or_tr):
                args.light_chain_fractions[ltmp] = 1. if ltmp == args.single_light_locus else 0.
        args.droplet_id_indices = utils.get_arg_list(args.droplet_id_indices, intify=True)
        if [args.droplet_id_separators, args.droplet_id_indices].count(None) not in [0, 2]:
            raise Exception('if you set either --droplet-id-separators or --droplet-id-indicies you need to set both of them (guessing defaults is proving to be too dangerous)')
        if args.plot_annotation_performance:
            print('  %s ignoring --plot-annotation-performance for paired clustering since it\'s going to be a bit fiddly to implement' % utils.color('yellow', 'warning'))
    else:
        if args.paired_indir is not None:
            raise Exception('need to set --paired-loci if --paired-indir is set')
        args.ig_or_tr = args.locus[:2]  # this maybe/probably doesn't need to happen
    if not args.paired_loci and (args.paired_indir is not None or args.paired_outdir is not None):
        raise Exception('--paired-loci must be set if either --paired-indir or --paired-outdir is set')
    if args.reverse_negative_strands and not args.paired_loci:
        raise Exception('--reverse-negative-strands has no effect unless --paired-loci is set (maybe need to run bin/split-loci.py separately?)')

    args.only_genes = utils.get_arg_list(args.only_genes)
    if args.paired_loci and args.action == 'simulate' and args.only_genes is not None:
        for _, l_locus in utils.locus_pairs[args.ig_or_tr]:
            if len([g for g in args.only_genes if utils.get_locus(g)==l_locus]) == 0:
                args.light_chain_fractions[l_locus] = 0
                args.light_chain_fractions[utils.get_single_entry([l for l in args.light_chain_fractions if l!=l_locus])] = 1
                print('  note: no %s genes among --only-genes %s, so setting --light-chain-fractions accordingly (%s)' % (l_locus, ':'.join(args.only_genes), args.light_chain_fractions))
    args.queries = utils.get_arg_list(args.queries)
    args.queries_to_include = utils.get_arg_list(args.queries_to_include)

    args.meta_info_to_emphasize = utils.get_arg_list(args.meta_info_to_emphasize, key_val_pairs=True)
    args.meta_emph_formats = utils.get_arg_list(args.meta_emph_formats, key_val_pairs=True)
    args.meta_info_bkg_vals = utils.get_arg_list(args.meta_info_bkg_vals)
    utils.meta_emph_arg_process(args)

    args.reco_ids = utils.get_arg_list(args.reco_ids)
    args.istartstop = utils.get_arg_list(args.istartstop, intify=True)
    if args.istartstop is not None:
        if args.istartstop[0] >= args.istartstop[1] or args.istartstop[0] < 0:
            raise Exception('invalid --istartstop specification: %d %d' % (args.istartstop[0], args.istartstop[1]))
    args.n_max_per_region = utils.get_arg_list(args.n_max_per_region, intify=True)
    if len(args.n_max_per_region) != 3:
        raise Exception('n-max-per-region should be of the form \'x:y:z\', but I got ' + str(args.n_max_per_region))
    args.write_additional_cluster_annotations = utils.get_arg_list(args.write_additional_cluster_annotations, intify=True)
    if args.write_additional_cluster_annotations is not None and len(args.write_additional_cluster_annotations) != 2:
        raise Exception('--write-additional-cluster-annotations must be specified as two numbers \'m:n\', but I got %s' % args.write_additional_cluster_annotations)
    args.extra_annotation_columns = utils.get_arg_list(args.extra_annotation_columns, choices=utils.extra_annotation_headers)
    args.cluster_size_bins = utils.get_arg_list(args.cluster_size_bins, floatify=True)

    args.input_metafnames = utils.get_arg_list(args.input_metafnames)
    if args.input_metafnames is not None and len(args.input_metafnames) > 1 and os.path.exists(':'.join(args.input_metafnames)):  # hackishly try to handle single meta fname paths that have ':'s in them (this will break if you *also* try to have multiple input metfnames
        print('  %s: guessing that --input-metafnames is only length one despite \':\'s: %s' % (utils.color('yellow', 'warning'), ':'.join(args.input_metafnames)))
        args.input_metafnames = [':'.join(args.input_metafnames)]

    args.cluster_indices = utils.get_arg_list(args.cluster_indices, intify_with_ranges=True)

    pair_allow_fmt = args.paired_loci and args.allowed_cdr3_lengths is not None  # if it's set for paired simulation, it'll be the key/val syntax
    args.allowed_cdr3_lengths = utils.get_arg_list(args.allowed_cdr3_lengths, key_val_pairs=pair_allow_fmt, intify_with_ranges=not pair_allow_fmt)
    if args.allowed_cdr3_lengths is not None:
        if pair_allow_fmt:
            missing_loci = [l for l in utils.sub_loci(args.ig_or_tr) if l not in args.allowed_cdr3_lengths]
            if len(missing_loci) > 0:
                print('  note: missing %d loc%s (%s) from --allowed-cdr3-lengths, so %s lengths will not be restricted' % (len(missing_loci), 'us' if len(missing_loci)==1 else 'i', ' '.join(missing_loci), 'its' if len(missing_loci)==1 else 'their'))
        else:
            print('  --allowed-cdr3-lengths: restricting to %s' % args.allowed_cdr3_lengths)

    args.region_end_exclusions = {r : [args.region_end_exclusion_length if ('%s_%s' % (r, e)) in utils.real_erosions else 0 for e in ['5p', '3p']] for r in utils.regions}
    args.region_end_exclusion_length = None  # there isn't really a big reason to set it to None, but this makes clear that I should only be using the dict version

    if args.dont_remove_framework_insertions:
        raise Exception('--dont-remove-framework-insertions is deprecated: the columns \'leader_seqs\' and \'c_gene_seqs\' are now always available')

    args.typical_genes_per_region_per_subject = utils.get_arg_list(args.typical_genes_per_region_per_subject, intify=True)
    if len(args.typical_genes_per_region_per_subject) != len(utils.regions):
        raise Exception('wrong length for --typical-genes-per-region-per-subject, has to be three')
    tmpfrac, ntmp = args.min_allele_prevalence_fraction, args.typical_genes_per_region_per_subject
    args.min_allele_prevalence_fractions = {r : tmpfrac * ntmp[utils.regions.index('v')] / ntmp[utils.regions.index(r)] for r in utils.regions}
    delattr(args, 'min_allele_prevalence_fraction')  # delete the non-plural version
    delattr(args, 'typical_genes_per_region_per_subject')  # and we don't need this any more either

    if args.annotation_clustering and args.action != 'annotate':
            raise Exception('action must be \'annotate\' for annotation clustering')
    args.naive_hamming_bounds = utils.get_arg_list(args.naive_hamming_bounds, floatify=True)
    if args.fast:
        args.naive_vsearch = True
        delattr(args, 'fast')
    if args.naive_hamming_cluster and args.naive_vsearch:
        raise Exception('can\'t specify both --naive-hamming-cluster and --naive-vsearch')
    if args.synthetic_distance_based_partition:
        if not args.is_simu:
            raise Exception('--synthetic-distance-based-partition: have to set --is-simu (and be running on a simulation file)')
        if not args.dont_calculate_annotations:
            print('  --synthetic-distance-based-partition: turning on --dont-calculate-annotations, since things may crash otherwise (because true naive seqs can be out of sync with sw info)')
            args.dont_calculate_annotations = True
        if not args.naive_vsearch:
            print('  --synthetic-distance-based-partition: turning on --naive-vsearch')
            args.naive_vsearch = True
            args.naive_hamming_cluster = False
    if args.small_clusters_to_ignore is not None:
        if '-' in args.small_clusters_to_ignore:
            lo, hi = [int(cluster_size) for cluster_size in args.small_clusters_to_ignore.split('-')]
            args.small_clusters_to_ignore = list(range(lo, hi + 1))
        else:
            args.small_clusters_to_ignore = utils.get_arg_list(args.small_clusters_to_ignore, intify=True)
    if not args.paired_loci and args.seed_unique_id is not None:  # if --paired-loci is set, there will be two seed uids/seqs, which requires totally different handling, so do it above
        args.seed_unique_id = args.seed_unique_id.strip()  # protect against the space you may put in front of it if it's got an initial minus sign (better way is to use an equals sign)
        if args.queries is not None and args.seed_unique_id not in args.queries:
            raise Exception('seed uid %s not in --queries %s' % (args.seed_unique_id, ' '.join(args.queries)))
        if args.random_seed_seq:
            raise Exception('can\'t specify both --seed-unique-id and --random-seed-seq')

        if args.queries_to_include is None:  # make sure the seed is in --queries-to-include
            args.queries_to_include = [args.seed_unique_id]
        elif args.seed_unique_id not in args.queries_to_include:
            args.queries_to_include = [args.seed_unique_id] + args.queries_to_include  # may as well put it first, I guess (?)

    args.extra_print_keys = utils.get_arg_list(args.extra_print_keys)

    if args.sw_debug is None:  # if not explicitly set, set equal to regular debug
        args.sw_debug = args.debug

    if args.only_genes is not None:
        for gene in args.only_genes:  # make sure they're all at least valid ig genes
            utils.split_gene(gene)

    if args.print_git_commit or args.action == 'version':
        utils.get_version_info(debug=True)
        if args.action == 'version':
            sys.exit(0)

    args.is_data = not args.is_simu  # whole code base uses is_data, this is better than changing all of that

    if args.collapse_duplicate_sequences and not args.is_data:
        print('  %s collapsing duplicates on simulation, which is often not a good idea since it makes keeping track of performance harder (e.g. purity/completeness of partitions is harder to calculate)' % utils.color('red', 'warning'))

    if args.simultaneous_true_clonal_seqs:
        if args.is_data:
            raise Exception('can only pass true clonal families to multi-hmm together on simulation and with --is-simu set')
        if args.n_simultaneous_seqs is not None:
            raise Exception('can\'t specify both --n-simultaneous-seqs and --simultaneous-true-clonal-seqs')
        if args.all_seqs_simultaneous:
            raise Exception('can\'t specify both --all-seqs-simultaneous and --simultaneous-true-clonal-seqs')
    if args.n_simultaneous_seqs is not None and args.all_seqs_simultaneous:
        raise Exception('doesn\'t make sense to set both --n-simultaneous-seqs and --all-seqs-simultaneous.')

    if args.no_indels or args.all_seqs_simultaneous: # or args.simultaneous_true_clonal_seqs:
        print('  forcing --gap-open-penalty to %d to prevent indels, since --no-indels or --all-seqs-simultaneous were set (you can also adjust this penalty directly)' % args.no_indel_gap_open_penalty)  # for all_seqs_simultaneous, we run the msa indel stuff so don't also want sw indels
        args.gap_open_penalty = args.no_indel_gap_open_penalty

    if args.locus is not None and 'tr' in args.locus and args.mutation_multiplier is None:
        args.mutation_multiplier = 0.

    if args.workdir is None:  # set default here so we know whether it was set by hand or not
        args.workdir = get_workdir(args.batch_system)
    else:
        args.workdir = args.workdir.rstrip('/')
    if os.path.exists(args.workdir):
        raise Exception('workdir %s already exists' % args.workdir)

    if args.batch_system == 'sge' and args.batch_options is not None:
        if '-e' in args.batch_options or '-o' in args.batch_options:
            print('%s --batch-options contains \'-e\' or \'-o\', but we add these automatically since we need to be able to parse each job\'s stdout and stderr. You can control the directory under which they\'re written with --workdir (which is currently %s).' % (utils.color('red', 'warning'), args.workdir))

    if args.outfname is not None and not args.presto_output and not args.airr_output and not args.generate_trees:
        if utils.getsuffix(args.outfname) not in ['.csv', '.yaml']:
            raise Exception('unhandled --outfname suffix %s' % utils.getsuffix(args.outfname))
        if utils.getsuffix(args.outfname) != '.yaml':
            print('  %s --outfname uses deprecated file format %s. This will still mostly work ok, but the new default .yaml format doesn\'t have to do all the string conversions by hand (so is less buggy), and includes annotations, partitions, and germline info in the same file (so you don\'t get crashes or inconsistent results if you don\'t keep track of what germline info goes with what output file).' % (utils.color('yellow', 'note:'), utils.getsuffix(args.outfname)))
        if args.action in ['view-annotations', 'view-partitions'] and utils.getsuffix(args.outfname) == '.yaml':
            raise Exception('have to use \'view-output\' action to view .yaml output files')

    if args.presto_output:
        if args.outfname is None:
            raise Exception('have to set --outfname if --presto-output is set')
        if args.action == 'annotate' and utils.getsuffix(args.outfname) != '.tsv':
            raise Exception('--outfname suffix has to be .tsv for annotation with --presto-output (got %s)' % utils.getsuffix(args.outfname))
        if args.action == 'partition' and utils.getsuffix(args.outfname) not in ['.fa', '.fasta']:
            raise Exception('--outfname suffix has to be .fa or .fasta for partitioning with --presto-output (got %s)' % utils.getsuffix(args.outfname))
        if args.aligned_germline_fname is None:
            assert args.locus is not None
            args.aligned_germline_fname = '%s/%s/imgt-alignments/%s.fa' % (args.default_initial_germline_dir, args.species, args.locus)
        if not os.path.exists(args.aligned_germline_fname):
            raise Exception('--aligned-germline-fname %s doesn\'t exist, but we need it in order to write presto output' % args.aligned_germline_fname)
    if not args.paired_loci and args.airr_output and not args.generate_trees:
        if args.outfname is None:
            if args.action != 'cache-parameters':
                print('  note: no --outfname set')
        else:
            if utils.getsuffix(args.outfname) == '.tsv':
                print('  note: writing only airr .tsv to %s' % args.outfname)
            elif utils.getsuffix(args.outfname) in ['.yaml', '.csv']:
                print('  note: writing both partis %s to %s and airr .tsv to %s' % (utils.getsuffix(args.outfname), args.outfname, utils.replace_suffix(args.outfname, '.tsv')))
            else:
                raise Exception('--outfname suffix has to be either .tsv or .yaml if --airr-output is set (got %s)' % utils.getsuffix(args.outfname))
    if args.airr_input:
        args.seq_column = 'sequence'
        args.name_column = 'sequence_id'

    if args.cluster_annotation_fname is None and args.outfname is not None and utils.getsuffix(args.outfname) == '.csv':  # if it wasn't set on the command line (<outfname> _was_ set), _and_ if we were asked for a csv, then use the old file name format
        args.cluster_annotation_fname = utils.insert_before_suffix('-cluster-annotations', args.outfname)

    if args.calculate_alternative_annotations and args.outfname is None and args.paired_outdir is None:
        raise Exception('have to specify --outfname in order to calculate alternative annotations')
    if args.calculate_alternative_annotations and args.outfname is not None and utils.getsuffix(args.outfname) == '.csv':
        raise Exception('old-style .csv output files no longer supported for alternative annotations (use yaml)')
    if args.subcluster_annotation_size == 'None':  # i want it turned on by default, but also to be able to turn it off on the command line
        args.subcluster_annotation_size = None
    else:
        args.subcluster_annotation_size = int(args.subcluster_annotation_size)  # can't set it in add_argument(), sigh
    if args.subcluster_annotation_size is not None and args.write_additional_cluster_annotations is not None:
        raise Exception('can\'t set --write-additional-cluster-annotations if --subcluster-annotation-size is also set (you get duplicate annotations, which confuses and crashes things)')
    if args.subcluster_annotation_size is not None and args.mimic_data_read_length:
        raise Exception('can\'t run subcluster annotations if --mimic-data-read-length is set, so need to set --subcluster-annotation-size to \'None\' (although this could be implemented/fixed)')
    if args.action == 'view-alternative-annotations' and args.persistent_cachefname is None and not args.paired_loci:  # handle existing old-style output
        assert args.outfname is not None
        if os.path.exists(utils.getprefix(args.outfname) + '-hmm-cache.csv'):
            args.persistent_cachefname = utils.getprefix(args.outfname) + '-hmm-cache.csv'  # written by bcrham, so has to be csv, not yaml

    if args.min_largest_cluster_size is not None and args.n_final_clusters is not None:
        print('  note: both --min-largest-cluster-size and --n-final-clusters are set, which means we\'ll stop clustering when *either* of their criteria are satisfied (not both)')  # maybe it should be both, but whatever
    if args.min_largest_cluster_size is not None or args.n_final_clusters is not None:
        if args.seed_unique_id is not None and args.n_procs == 1:
            raise Exception('--n-procs must be set to greater than 1 if --seed-unique-id, and either --min-largest-cluster-size or --n-final-clusters, are set (so that a second clustering iteration is run after removing)')  # yes, this could also be fixed by making the algorithm that decides when to stop clustering smarter, but that would be hard
        args.n_partitions_to_write = 999  # need to make sure to get all the ones after the best partition, and this is a somewhat hackey way to do that
        print('  note: setting --n-partitions-to-write to 999 since either --min-largest-cluster-size or --n-final-clusters was set')  # bcrham argument parser barfs if you use sys.maxint

    args.existing_output_run_cfg = utils.get_arg_list(args.existing_output_run_cfg, choices=['single', 'paired', 'merged', 'no-fake-paired'])
    if args.existing_output_run_cfg is None:
        args.existing_output_run_cfg = []

    if args.action == 'infer-trees':  # don't need to make a new action, but calculate only the simplest smetric so it's faster
        args.action = 'get-selection-metrics'
        args.selection_metrics_to_calculate = 'lbi'
        args.selection_metric_plot_cfg = 'lb-scatter'
        sys.argv[sys.argv.index('infer-trees')] = 'get-selection-metrics'  # sys.argv gets used for arg manipulation in paired stuff
        print('  note: changing action \'infer-trees\' to \'get-selection-metrics\'')
    if args.action == 'get-selection-metrics' or args.get_selection_metrics:
        if args.paired_loci:
            if args.paired_outdir is None and args.selection_metric_fname is None:
                print('    %s calculating selection metrics, but neither --paired-outdir nor --selection-metric-fname were set, which means nothing will be written to disk' % utils.color('yellow', 'warning'))
            elif args.selection_metric_fname is None:
                if args.add_selection_metrics_to_outfname:
                    print('  note: --add-selection-metrics-to-outfname has no effect on final output files when --paired-loci is set since there isn\'t a unique output file (although if --existing-output-run-cfg is set they can get written to those files)')
                args.selection_metric_fname = treeutils.smetric_fname(args.paired_outdir)
        else:
            if args.outfname is None and args.selection_metric_fname is None:
                    print('    %s calculating selection metrics, but neither --outfname nor --selection-metric-fname were set, which means nothing will be written to disk' % utils.color('yellow', 'warning'))
            elif args.selection_metric_fname is None and args.action == 'get-selection-metrics' and not args.add_selection_metrics_to_outfname:
                args.selection_metric_fname = treeutils.smetric_fname(args.outfname)
    args.selection_metrics_to_calculate = utils.get_arg_list(args.selection_metrics_to_calculate, choices=treeutils.selection_metrics)
    args.selection_metric_plot_cfg = utils.get_arg_list(args.selection_metric_plot_cfg, choices=treeutils.all_plot_cfg)
    args.partition_plot_cfg = utils.get_arg_list(args.partition_plot_cfg, choices=utils.all_ptn_plot_cfg)
    if args.invert_affinity and args.affinity_key is None:
        raise Exception('--affinity-key must be set if setting --invert-affinity')
    args.extra_daffy_metrics = utils.get_arg_list(args.extra_daffy_metrics)
    if args.extra_daffy_metrics is not None:
        print('  --extra-daffy-metrics: adding %d metrics to treeutils.daffy_metrics (%s)' % (len(args.extra_daffy_metrics), ':'.join(args.extra_daffy_metrics)))
        treeutils.daffy_metrics += args.extra_daffy_metrics
    args.tree_inference_outdir = None  # NOTE not an actual command line arg
    if args.outfname is None and args.paired_outdir is None:
        args.tree_inference_outdir = None
    else:
        args.tree_inference_outdir = utils.fpath(utils.getprefix(utils.non_none([args.outfname, args.paired_outdir])))
    if args.tree_inference_subdir is not None:
        args.tree_inference_outdir = '%s/%s' % (utils.non_none([args.tree_inference_outdir, os.getcwd()]), args.tree_inference_subdir)
    # if args.tree_inference_outdir is not None:
    #     print('  note: set tree inference outdir to %s' % args.tree_inference_outdir)
    if args.is_simu and args.tree_inference_method is not None:
        raise Exception('can\'t set both --is-simu and --tree-inference-method, since if --is-simu is set we use the simulation tree without inferring anything')
    args.mutation_label_cfg = utils.get_arg_list(args.mutation_label_cfg, choices=['all', 'leaf', 'mut-strs', 'short'])

    if args.plot_annotation_performance:
        if args.plotdir is None and args.print_n_worst_annotations is None:
            raise Exception('doesn\'t make sense to set --plot-annotation-performance but not either of --plotdir or --print-n-worst-annotations (we\'ll spend all the cycles counting things up but then they\'ll just disappear from memory without being recorded).')
        if not args.is_simu:
            raise Exception('can\'t plot performance unless --is-simu is set (and this is simulation)')
    if args.print_n_worst_annotations is not None and not args.plot_annotation_performance:
        raise Exception('--plot-annotation-performance must be set if you\'re setting --print-worst-annotations')
    if not args.paired_loci and (args.action=='plot-partitions' or args.action=='annotate' and args.plot_partitions) and args.plotdir is None:
        raise Exception('--plotdir must be specified if plotting partitions')
    if args.input_partition_fname is not None:
        if args.action not in ['annotate', 'partition']:
            print('  %s --input-partition-fname only makes sense/has an effect for actions \'annotate\' and \'partition\' (at least at the moment)' % utils.wrnstr())
    if args.action == 'annotate' and args.plot_partitions and args.input_partition_fname is None:  # could set this up to use e.g. --simultaneous-true-clonal-seqs as well, but it can't atm
        print('  %s running annotate with --plot-partitions, but --input-partition-fname is not set, which likely means the partitions will be trivial/singleton partitions' % utils.color('yellow', 'warning'))

    if args.make_per_gene_per_base_plots and not args.make_per_gene_plots:  # the former doesn't do anything unless the latter is turned on
        args.make_per_gene_plots = True

    if not args.paired_loci and args.count_correlations and args.action != 'cache-parameters' and not args.count_parameters:
        raise Exception('--count-correlations has no effect for action \'%s\' unless --paired-loci is set, or you also turn on --count-parameters' % args.action)

    if args.action == 'simulate':
        if args.n_trees is None and not args.paired_loci:
            args.n_trees = max(1, int(float(args.n_sim_events) / args.n_procs))
        if args.n_procs > args.n_sim_events:
            print('  note: reducing --n-procs to %d (was %d) so it isn\'t bigger than --n-sim-events' % (args.n_sim_events, args.n_procs))
            args.n_procs = args.n_sim_events
        if args.n_max_queries != -1:
            print('  note: --n-max-queries is not used when simulating (use --n-sim-events to set the simulated number of rearrangemt events)')

        if args.outfname is None and args.paired_outdir is None:
            print('  note: no %s specified, so nothing will be written to disk' % ('--paired-outdir' if args.paired_loci else '--outfname'))
            args.outfname = get_dummy_outfname(args.workdir)  # hackey, but otherwise I have to rewrite the whole run_simulation() in bin/partis to handle None type outfname

        if args.airr_output:
            raise Exception('--airr-output isn\'t implemented for \'simulate\', but you can instead convert to airr tsv afterwards with for instance \'parse-output simu.yaml simu.tsv --airr-output\' or \'parse-output --paired <paired-simu-dir> <paired-simu-dir> --airr-output\'')

        if args.simulate_from_scratch:
            args.rearrange_from_scratch = True
            args.mutate_from_scratch = True
        if args.rearrange_from_scratch and not args.force_dont_generate_germline_set:  # i would probably just default to always generating germline sets when rearranging from scratch, but bin/test-germline-inference.py (and any other case where you want to dramatically restrict the germline set) really argue for a way to force just using the genes in the germline dir
            args.generate_germline_set = True
        if args.flat_mute_freq or args.same_mute_freq_for_all_seqs:
            assert args.mutate_from_scratch
        if args.mutate_from_scratch and not args.no_per_base_mutation:
            print('  note: setting --no-per-base-mutation since --mutate-from-scratch was set')
            args.no_per_base_mutation = True

        if args.dry_run and any(a is not None for a in [args.bulk_data_fraction, args.mean_cells_per_droplet, args.fraction_of_reads_to_remove]):
            raise Exception('can\'t use --dry-run to rerun locus combination parts of paired simulation when any of --bulk-data-fraction, --mean-cells-per-droplet, or --fraction-of-reads-to-remove are set (since the single chain files that we\'d read with --dry-run get rewritten after modification from any of these args)')
        # end result of this block: shm/reco parameter dirs are set (unless we're doing their bit from scratch), --parameter-dir is set to None (and if --parameter-dir was set but shm/reco were _not_ set, we've just used --parameter-dir for either/both as needed)
        if args.parameter_dir is not None:
            if args.rearrange_from_scratch or args.mutate_from_scratch:
                raise Exception('can\'t set --parameter-dir if rearranging or mutating from scratch (use --reco-parameter-dir and/or --shm-parameter-dir)')
            if args.reco_parameter_dir is not None or args.shm_parameter_dir is not None:
                raise Exception('can\'t set --parameter-dir if either --reco-parameter-dir or --shm-parameter-dir are also set')
            args.reco_parameter_dir = args.parameter_dir
            args.shm_parameter_dir = args.parameter_dir
            args.parameter_dir = None
        if args.rearrange_from_scratch and args.reco_parameter_dir is not None:
            raise Exception('doesn\'t make sense to set both --rearrange-from-scratch and --reco-parameter-dir')
        if args.mutate_from_scratch and args.shm_parameter_dir is not None:
            raise Exception('doesn\'t make sense to set both --mutate-from-scratch and --shm-parameter-dir')
        if args.reco_parameter_dir is None and not args.rearrange_from_scratch:
            raise Exception('have to either set --rearrange-from-scratch or --reco-parameter-dir (or --simulate-from-scratch)')
        if args.shm_parameter_dir is None and not args.mutate_from_scratch:
            raise Exception('have to either set --mutate-from-scratch or --shm-parameter-dir (or --simulate-from-scratch)')

        if args.generate_germline_set and not args.rearrange_from_scratch:
            raise Exception('can only --generate-germline-set if also rearranging from scratch (set --rearrange-from-scratch)')

        if args.constant_number_of_leaves and args.n_leaf_distribution is not None:
            raise Exception('--n-leaf-distribution has no effect if --constant-number-of-leaves is set (but both were set)')

        if args.mean_cells_per_droplet is not None and args.mean_cells_per_droplet == 'None':  # it's nice to be able to unset this from the command line (even though the defualt is off) e.g. when calling from cf-paired-loci.py
            args.mean_cells_per_droplet = None
        if args.mean_cells_per_droplet is not None:  # ick
            args.mean_cells_per_droplet = float(args.mean_cells_per_droplet)

        if args.indel_frequency > 0.:
            if args.indel_frequency < 0. or args.indel_frequency > 1.:  # < 0. is there in case previous if statement or default value change
                raise Exception('--indel-frequency must be in [0., 1.] (got %f)' % args.indel_frequency)
        args.n_indels_per_indeld_seq = utils.get_arg_list(args.n_indels_per_indeld_seq, intify=True)
        if args.indel_location not in [None, 'v', 'cdr3']:  # if it's not default (None) and also not set to a region ('v', 'cdr3')
            if int(args.indel_location) in range(500):
                args.indel_location = int(args.indel_location)
                if any(n > 1 for n in args.n_indels_per_indeld_seq):
                    print('  note: removing entries from --n-indels-per-indeld-seq (%s), since --indel-location was set to a single position.' % [n for n in args.n_indels_per_indeld_seq if n > 1])
                    args.n_indels_per_indeld_seq = [n for n in args.n_indels_per_indeld_seq if n <= 1]
            else:
                raise Exception('--indel-location \'%s\' neither one of None, \'v\' or \'cdr3\', nor an integer less than 500' % args.indel_location)

        if args.simulation_germline_dir is not None:
            raise Exception('--simulation-germline-dir has no effect on simulation (maybe you meant --initial-germline-dir?)')
        if args.generate_germline_set:
            args.snp_positions = None  # if you want to control the exact positions, you have to use bin/test-germline-inference.py
            args.indel_positions = None
            process_gls_gen_args(args)

        if args.generate_trees:
            assert args.n_procs == 1  # not set up to handle output, and also no need

        if args.treefname is not None:
            raise Exception('--treefname was set for simulation action (probably meant to use --input-simulation-treefname)')
        if args.fraction_of_reads_to_remove is not None:
            assert args.fraction_of_reads_to_remove >= 0. and args.fraction_of_reads_to_remove < 1.

        if args.remove_nonfunctional_seqs and args.paired_loci and not args.mutate_stop_codons:
            print('  note: removing nonfunctional seqs with paired loci may remove a very large fraction of sequences unless you also set --mutate-stop-codons')

        # ----------------------------------------------------------------------------------------
        def process_corr_values(cvals, estr=''):
            if cvals is None:
                return cvals
            if not args.rearrange_from_scratch:
                raise Exception('--%scorrelation-values has no effect unless --rearrange-from-scratch (or --simulate-from-scratch) is set' % estr)
            cvals = utils.get_arg_list(cvals, key_val_pairs=True, floatify=True)
            for kpstr in list(cvals.keys()):
                ppair = tuple(kpstr.split('.'))
                avail_corrs = utils.available_simu_correlations if estr=='' else utils.paired_available_simu_correlations
                if ppair not in avail_corrs:
                    raise Exception('parameter pair %s in --%scorrelation-values not among allowed ones (%s)' % (ppair, estr, ', '.join(str(p) for p in avail_corrs)))
                if cvals[kpstr] < 0. or cvals[kpstr] > 1.:
                    raise Exception('correlation value %f not in [0, 1]' % cvals[kpstr])
                cvals[ppair] = cvals[kpstr]
                del cvals[kpstr]
            return cvals
        # ----------------------------------------------------------------------------------------
        args.correlation_values = process_corr_values(args.correlation_values)
        args.paired_correlation_values = process_corr_values(args.paired_correlation_values, estr='paired-')
        if args.correlation_values is not None and args.paired_correlation_values is not None:
            raise Exception('can\'t yet handle single-chain --correlation-values and --paired-correlation-values at the same time')
        if args.locus is not None and not utils.has_d_gene(args.locus) and args.paired_correlation_values is not None and args.heavy_chain_event_fname is None:
            raise Exception('if --paired-correlation-values is set for light chain, you have to also set --heavy-chain-event-fname')

        if args.allowed_cdr3_lengths is not None and not args.paired_loci and args.rearrange_from_scratch and any(l < 30 for l in args.allowed_cdr3_lengths):  # ick
            factor = 3
            print('  note: asked for cdr3 lengths <30 with --rearrange-from-scratch set, so increasing scratch mean deletion lengths and decreasing scratch mean insertion lengths by factor of %.1f:' % factor)
            for tdict, tfact in zip([utils.scratch_mean_erosion_lengths[args.locus], utils.scratch_mean_insertion_lengths[args.locus]], [factor, 1./factor]):
                for tstr in tdict:
                    print('    %8s  %5.1f --> %5.1f' % (tstr, tdict[tstr], tfact * tdict[tstr]))
                    tdict[tstr] *= tfact

    if args.parameter_dir is not None and not args.paired_loci:  # if we're splitting loci, this isn't the normal parameter dir, it's a parent of that
        args.parameter_dir = args.parameter_dir.rstrip('/')
        if os.path.exists(args.parameter_dir):
            pdirs = [d for d in os.listdir(args.parameter_dir) if os.path.isdir(d)]
            if len(pdirs) > 0 and len(set(pdirs) & set(utils.parameter_type_choices)) == 0:
                raise Exception('couldn\'t find any expected parameter types (i.e. subdirs) in --parameter-dir \'%s\'. Allowed types: %s, found: %s. Maybe you added the parameter type to the parameter dir path?' % (args.parameter_dir, ' '.join(utils.parameter_type_choices), ' '.join(os.listdir(args.parameter_dir))))

    if os.path.exists(args.default_initial_germline_dir + '/' + args.species):  # ick that is hackey
        args.default_initial_germline_dir += '/' + args.species

    if args.species != 'human' and not args.allele_cluster and args.action == 'cache-parameters':
        print('  non-human species \'%s\', turning on allele clustering' % args.species)
        args.allele_cluster = True

    if args.n_max_snps is not None and args.n_max_mutations_per_segment is not None:
        if args.n_max_snps > args.n_max_mutations_per_segment - 10:
            raise Exception('--n-max-snps should be at least ten less than --n-max-mutations-per-segment, but I got %d and %d' % (args.n_max_snps, args.n_max_mutations_per_segment))

    if args.leave_default_germline:
        args.dont_remove_unlikely_alleles = True
        args.allele_cluster = False
        args.dont_find_new_alleles = True

    if args.action not in actions_not_requiring_input and [args.infname, args.paired_indir].count(None) == 2:
        if args.paired_loci:
            raise Exception('--infname or --paired-indir is required for action \'%s\' with --paired-loci' % args.action)
        else:
            raise Exception('--infname is required for action \'%s\'' % args.action)

    if args.action == 'get-linearham-info':
        if args.linearham_info_fname is None:  # for some reason setting required=True isn't working
            raise Exception('have to specify --linearham-info-fname')
        if args.sw_cachefname is None and args.parameter_dir is None:
            raise Exception('have to specify --sw-cachefname or --parameter-dir, since we need sw info to calculate linearham inputs')
        if args.extra_annotation_columns is None or 'linearham-info' not in args.extra_annotation_columns:
            args.extra_annotation_columns = utils.add_lists(args.extra_annotation_columns, ['linearham-info'])
