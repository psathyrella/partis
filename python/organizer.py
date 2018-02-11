import os
import numpy
import time
import random
import copy
import sys
import multiprocessing

import utils
import glutils
from partitiondriver import PartitionDriver

# ----------------------------------------------------------------------------------------
def default_parameter_dir(args):
    if args.infname is not None:
        label = args.infname[ : args.infname.rfind('.')]  # probably ok not to use the full path here
        label = label.replace('/', '_')
    else:
        label = 'xxx-dummy-xxx'  # gah, this is a terrible way to do this. But I can't set it to None, and I need some way to communicate that it isn't a valid pdir, so...
    pdir = '_output/' + label
    return pdir

# ----------------------------------------------------------------------------------------
def set_default_parameter_dir(args):
    if args.parameter_dir is None and not args.rearrange_from_scratch and not args.mutate_from_scratch:
        args.parameter_dir = default_parameter_dir(args)

# ----------------------------------------------------------------------------------------
def run_simulation(args):
    from recombinator import Recombinator

    if args.batch_system is not None and args.n_procs > 1 and not args.subsimproc:
        print '  %s setting subsimproc' % utils.color('red', 'warning')
        args.subsimproc = True
    if args.n_trees is None:
        args.n_trees = max(1, int(float(args.n_sim_events) / args.n_procs))
    if args.outfname is None:
        raise Exception('have to specify --outfname for simulation')
    if args.n_max_queries != -1:
        print '  note: --n-max-queries is not used when simulating (use --n-sim-events to set the simulated number of rearrangemt events)'

    if args.parameter_dir is None:
        if not args.rearrange_from_scratch and not args.mutate_from_scratch:
            raise Exception('either --parameter-dir must be specified, or one of the scratch simulation options')
        if args.mutate_from_scratch and not args.rearrange_from_scratch:
            raise Exception('haven\'t yet implemented mutating from scratch without rearranging from scratch')  # and maybe not ever
    else:
        if args.rearrange_from_scratch or args.mutate_from_scratch or args.generate_germline_set:
            raise Exception('you can\'t specify --parameter-dir if you also set either of the scratch options or --generate-germline-set')
        if args.initial_germline_dir is not None:
            raise Exception('you can\'t specify both --parameter-dir and --inital-germline-dir (germline info has to be read from the parameter dir)')

    utils.prep_dir(args.workdir)

    default_prevalence_fname = args.workdir + '/allele-prevalence.csv'
    if args.generate_germline_set and args.allele_prevalence_fname is None:
        args.allele_prevalence_fname = default_prevalence_fname

    if args.initial_germline_dir is None:  # if the gl info dir wasn't explicitly set on the command line, we have to decide what to use
        if args.parameter_dir is not None:  # if --parameter-dir was explicitly set, assume there's gl info in this --parameter-dir
            input_gldir = args.parameter_dir + '/' + args.parameter_type + '/' + glutils.glfo_dir
        else:  # otherwise use the default
            input_gldir = args.default_initial_germline_dir
    else:
        input_gldir = args.initial_germline_dir
    assert len(args.loci) == 1  # needs to be implemented
    glfo = glutils.read_glfo(input_gldir, args.locus, only_genes=args.only_genes)
    if not args.im_a_subproc and args.rearrange_from_scratch and args.generate_germline_set:
        glutils.generate_germline_set(glfo, args.n_genes_per_region, args.n_sim_alleles_per_gene, args.min_sim_allele_prevalence_freq, args.allele_prevalence_fname)  # NOTE removes unwanted genes from <glfo>
        print '  writing generated germline set to %s/' % args.outfname.replace('.csv', '-glfo')
        glutils.write_glfo(args.outfname.replace('.csv', '-glfo'), glfo)
    if args.subsimproc:
        working_gldir = args.workdir + '/' + glutils.glfo_dir
        glutils.write_glfo(working_gldir, glfo)

    # ----------------------------------------------------------------------------------------
    def get_subproc_cmd_str(n_events, iproc, workdir, outfname):
        clist = copy.deepcopy(sys.argv)
        utils.remove_from_arglist(clist, '--n-procs', has_arg=True)
        utils.remove_from_arglist(clist, '--subsimproc')
        clist.append('--im-a-subproc')
        utils.replace_in_arglist(clist, '--seed', str(args.seed + iproc))
        utils.replace_in_arglist(clist, '--workdir', workdir)
        utils.replace_in_arglist(clist, '--outfname', outfname)
        utils.replace_in_arglist(clist, '--n-sim-events', str(n_events))
        utils.replace_in_arglist(clist, '--initial-germline-dir', working_gldir)
        if args.allele_prevalence_fname is not None:
            utils.replace_in_arglist(clist, '--allele-prevalence-fname', args.allele_prevalence_fname)
        return ' '.join(clist)

    # ----------------------------------------------------------------------------------------
    def make_events(n_events, iproc, workdir, outfname, random_ints):
        reco = Recombinator(args, glfo, seed=args.seed+iproc, workdir=workdir, outfname=outfname)
        for ievt in range(n_events):
            reco.combine(random_ints[ievt])
        # reco.print_validation_values()

    def get_workdir(iproc):
        return args.workdir + '/sub-' + str(iproc)
    def get_outfname(iproc):
        return get_workdir(iproc) + '/' + os.path.basename(args.outfname)

    if not args.im_a_subproc:
        print 'simulating'

    set_default_parameter_dir(args)

    n_per_proc = int(float(args.n_sim_events) / args.n_procs)

    # generate all the random seeds NOTE these aren't used if <args.subsimproc> is set, i.e. results will be different with and without <args.subsimproc> set.
    all_random_ints = []
    for iproc in range(args.n_procs):  # have to do it all at once, 'cause each of the subprocesses is going to reset its seed and god knows what happens to our seed at that point
        all_random_ints.append([random.randint(0, numpy.iinfo(numpy.int32).max) for i in range(n_per_proc)])

    if args.n_procs == 1:  # multiprocessing is kind of messy
        make_events(n_per_proc, 0, args.workdir, args.outfname, all_random_ints[0])
    else:  # start the processes and wait for 'em to finish
        cmdfos = [{'cmd_str' : get_subproc_cmd_str(n_per_proc, iproc, get_workdir(iproc), get_outfname(iproc)) if args.subsimproc else None,
                   'workdir' : get_workdir(iproc),
                   'logdir' : args.workdir + '/log-' + str(iproc),  # have to be different than <workdirs> since ./bin/partis barfs if its workdir already exists (as it should)
                   'outfname' : get_outfname(iproc)}
                  for iproc in range(args.n_procs)]
        if args.subsimproc:
            utils.run_cmds(cmdfos, batch_system=args.batch_system, batch_options=args.batch_options, batch_config_fname=args.batch_config_fname, debug='print')
        else:
            for iproc in range(args.n_procs):
                proc = multiprocessing.Process(target=make_events, args=(n_per_proc, iproc, cmdfos[iproc]['workdir'], cmdfos[iproc]['outfname'], all_random_ints[iproc]))
                proc.start()
            while len(multiprocessing.active_children()) > 0:
                sys.stdout.flush()
                time.sleep(1)

        # check and merge output
        n_total_events = 0
        for iproc in range(args.n_procs):
            fname = cmdfos[iproc]['outfname']
            if not os.path.exists(fname):
                raise Exception('output simulation file %s d.n.e.' % fname)
            n_events = int(check_output('grep -v unique_id ' + fname + ' | awk -F, \'{print $2}\' | uniq | wc -l', shell=True).split()[0])  # if *adjacent* events have the same rearrangement parameters (and hence reco id), they'll show up as the same event. It has happened!
            n_total_events += n_events
            if n_events < n_per_proc - 3:  # give ourselves a little breathing room for adjacent events that have the same rearrangement parameters (it's only happened once, ever, so maybe 3 is enough?)
                raise Exception('only found %d events (expected %d) in output file %s' % (n_events, n_per_proc, fname))
        utils.merge_csvs(args.outfname, [cmdfos[i]['outfname'] for i in range(args.n_procs)])
        print '   read %d events from %d files' % (n_total_events, args.n_procs)

    if not args.im_a_subproc and args.allele_prevalence_fname is not None:  # check final prevalence freqs
        glutils.check_allele_prevalence_freqs(args.outfname, glfo, args.allele_prevalence_fname, only_region='v')
        if args.allele_prevalence_fname == default_prevalence_fname:
            os.remove(default_prevalence_fname)

    if args.subsimproc:
        glutils.remove_glfo_files(working_gldir, args.locus)
        for iproc in range(args.n_procs):
            os.rmdir(cmdfos[iproc]['logdir'])

    if not args.im_a_subproc:
        try:
            os.rmdir(args.workdir)
        except OSError:
            raise Exception('workdir (%s) not empty: %s' % (args.workdir, ' '.join(os.listdir(args.workdir))))  # hm... you get weird recursive exceptions if you get here. Oh, well, it still works

# ----------------------------------------------------------------------------------------
def run_partitiondriver(args):
    if len(args.n_max_per_region) != 3:
        raise Exception('n-max-per-region should be of the form \'x:y:z\', but I got ' + str(args.n_max_per_region))
    if len(args.initial_match_mismatch) != 2:
        raise Exception('--initial-match-mismatch should be of the form \'match:mismatch\', but I got ' + str(args.n_max_per_region))
    if args.infname is None:
        raise Exception('--infname is required for the \'%s\' action' % args.action)

    set_default_parameter_dir(args)

    if args.action == 'cache-parameters':
        parter = PartitionDriver(args, 'cache-parameters', args.default_initial_germline_dir if args.initial_germline_dir is None else args.initial_germline_dir)  # if we're caching parameters, read initial gl info from <args.initial_germline_dir> (then the final glfo gets written to the parameter dir)
        parter.cache_parameters()
        parter.clean()
        return

    if not os.path.exists(args.parameter_dir):
        print '  parameter dir \'%s\' does not exist, so caching a new set of parameters before running action \'%s\'' % (args.parameter_dir, args.action)
        newargv = copy.deepcopy(sys.argv)
        for argconf in [ac for ac in subargs[args.action] if ac['name'] in newargv]:  # remove args that only make sense for <args.action>
            index = newargv.index(argconf['name'])
            newargv.remove(argconf['name'])
            if 'action' not in argconf['kwargs'] or argconf['kwargs']['action'] != 'store_true':  # also remove the argument's argument
                newargv.pop(index)
        newargv = newargv[:1] + ['cache-parameters', ] + newargv[2:]
        check_call(newargv)

    if args.action == 'annotate':
        parter = PartitionDriver(args, args.action, args.parameter_dir + '/' + args.parameter_type + '/' + glutils.glfo_dir)
        parter.annotate()
    elif args.action == 'partition':
        if args.n_partition_subsets is None:
            parter = PartitionDriver(args, args.action, args.parameter_dir + '/' + args.parameter_type + '/' + glutils.glfo_dir)
            parter.partition()
        else:
            if args.abbreviate:
                raise Exception('not supported (different runs have the same uid for different sequences, which breaks everything)')
            old_outfname = args.outfname
            if args.n_max_queries == -1:
                raise Exception('--n-max-queries must be set in order to use --n-partition-subsets')
            n_per_subset = args.n_max_queries / args.n_partition_subsets
            args.n_max_queries = -1
            outfnames = [args.workdir + '/outfiles/partition-' + str(isub) + '.csv' for isub in range(args.n_partition_subsets)]
            istart = 0
            for isub in range(args.n_partition_subsets):
                istop = istart + n_per_subset
                args.outfname = outfnames[isub]
                args.istartstop = [istart, istop]  # and here I break my rule about not modifying <args>, sigh.
                parter = PartitionDriver(args, args.action, args.parameter_dir + '/' + args.parameter_type + '/' + glutils.glfo_dir)
                parter.partition()
                istart = istop
            partplotter = PartitionPlotter()
            if args.plotdir is not None:
                partplotter.plot(args.plotdir + '/partitions', infiles=outfnames, only_csv=args.only_csv_plots)
            for fn in outfnames:
                os.remove(fn)
            os.rmdir(os.path.dirname(outfnames[0]))
    else:
        raise Exception('bad action ' + args.action)
    parter.clean()

# ----------------------------------------------------------------------------------------
def view_existing_output(args):
    if args.outfname is None:
        if args.infname is not None:  # can specify either outfname or infname
            args.outfname = args.infname
            args.infname = None  # otherwise seqfileopener tries to read it (you can specify *both* of them on the command line, though, if you want this step to know about the reco info)
        else:
            raise Exception('--outfname is required for action \'%s\'' % args.action)
    if not os.path.exists(args.outfname):
        raise Exception('--outfname \'%s\' does not exist' % args.outfname)

    set_default_parameter_dir(args)  # all the calls to this function could probably be combined

    if args.initial_germline_dir is None:
        if 'xxx-dummy-xxx' in args.parameter_dir:
            args.initial_germline_dir = args.default_initial_germline_dir
        else:
            args.initial_germline_dir = args.parameter_dir + '/' + args.parameter_type + '/' + glutils.glfo_dir

    parter = PartitionDriver(args, args.action, args.initial_germline_dir)
    if args.action == 'view-annotations':
        parter.read_existing_annotations(debug=True)
    elif args.action == 'view-partitions':
        parter.read_existing_partitions(debug=True)
    elif args.action == 'plot-partitions':
        parter.plot_existing_partitions()
    elif args.action == 'view-alternative-naive-seqs':
        parter.view_alternative_naive_seqs()
    else:
        assert False
    parter.clean()

