import time
import copy
import sys
import math
import re
import os
import itertools
import operator
import pysam
import contextlib
from collections import OrderedDict
import csv

import utils
import glutils
from opener import opener
from parametercounter import ParameterCounter
from performanceplotter import PerformancePlotter
from allelefinder import AlleleFinder
from alleleremover import AlleleRemover
from clusterpath import ClusterPath

# ----------------------------------------------------------------------------------------
class Waterer(object):
    """ Run smith-waterman on the query sequences in <infname> """
    def __init__(self, args, input_info, reco_info, glfo, count_parameters=False, parameter_out_dir=None, remove_less_likely_alleles=False, find_new_alleles=False, plot_performance=False, simglfo=None, itry=None):
        self.args = args
        self.input_info = input_info
        self.reco_info = reco_info
        self.glfo = glfo
        self.simglfo = simglfo
        self.parameter_out_dir = parameter_out_dir
        self.debug = self.args.debug if self.args.sw_debug is None else self.args.sw_debug

        self.max_insertion_length = 35  # if an insertion is longer than this, we skip the proposed annotation (i.e. rerun it)
        self.absolute_max_insertion_length = 120  # but if it's longer than this, we always skip the annotation

        self.remaining_queries = set([q for q in self.input_info.keys()])  # we remove queries from this set when we're satisfied with the current output (in general we may have to rerun some queries with different match/mismatch scores)
        self.new_indels = 0  # number of new indels that were kicked up this time through

        self.match_mismatch = copy.deepcopy(self.args.initial_match_mismatch)  # don't want to modify it!
        self.gap_open_penalty = self.args.gap_open_penalty  # not modifying it now, but just to make sure we don't in the future

        self.info = {}
        self.info['queries'] = []  # list of queries that *passed* sw, i.e. for which we have information
        self.info['all_best_matches'] = set()  # every gene that was a best match for at least one query
        self.info['all_matches'] = {r : set() for r in utils.regions}  # every gene that was *any* match (up to <self.args.n_max_per_region[ireg]>) for at least one query NOTE there is also an 'all_matches' in each query's info
        self.info['indels'] = {}  # NOTE if we find shm indels in a sequence, we store the indel info in here, and rerun sw with the reversed sequence (i.e. <self.info> contains the sw inference on the reversed sequence -- if you want the original sequence, get that from <self.input_info>)
        self.info['mute-freqs'] = None  # kind of hackey... but it's to allow us to keep this info around when we don't want to keep the whole waterer object around, and we don't want to write all the parameters to disk

        self.nth_try = 1
        self.skipped_unproductive_queries, self.kept_unproductive_queries = set(), set()

        self.my_gldir = self.args.workdir + '/' + glutils.glfo_dir

        self.alremover, self.alfinder, self.pcounter, self.true_pcounter, self.perfplotter = None, None, None, None, None
        if remove_less_likely_alleles:
            self.alremover = AlleleRemover(self.glfo, self.args, AlleleFinder(self.glfo, self.args, itry=0))
        if find_new_alleles:  # NOTE *not* the same as <self.args.find_new_alleles>
            self.alfinder = AlleleFinder(self.glfo, self.args, itry)
        if count_parameters:  # NOTE *not* the same as <self.args.cache_parameters>
            self.pcounter = ParameterCounter(self.glfo, self.args)
            if not self.args.is_data:
                self.true_pcounter = ParameterCounter(self.simglfo, self.args)
        if plot_performance:  # NOTE *not* the same as <self.args.plot_performance>
            self.perfplotter = PerformancePlotter('sw')

        if not os.path.exists(self.args.ig_sw_binary):
            raise Exception('ig-sw binary d.n.e: %s' % self.args.ig_sw_binary)

    # ----------------------------------------------------------------------------------------
    def run(self, cachefname=None):
        start = time.time()
        base_infname = 'query-seqs.fa'
        base_outfname = 'query-seqs.sam'
        sys.stdout.flush()

        n_procs = self.args.n_fewer_procs
        initial_queries_per_proc = float(len(self.remaining_queries)) / n_procs
        while len(self.remaining_queries) > 0:  # we remove queries from <self.remaining_queries> as we're satisfied with their output
            if self.nth_try > 1 and float(len(self.remaining_queries)) / n_procs < initial_queries_per_proc:
                n_procs = int(max(1., float(len(self.remaining_queries)) / initial_queries_per_proc))
            self.write_vdjalign_input(base_infname, n_procs)
            self.execute_commands(base_infname, base_outfname, n_procs)
            self.read_output(base_outfname, n_procs)
            if self.nth_try > 3:
                break
            self.nth_try += 1  # it's set to 1 before we begin the first try, and increases to 2 just before we start the second try

        self.finalize(cachefname)
        print '        water time: %.1f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def read_cachefile(self, cachefname):
        start = time.time()
        print '        reading sw results from %s' % cachefname

        if os.path.exists(cachefname.replace('.csv', '-glfo')):
            self.glfo = glutils.read_glfo(cachefname.replace('.csv', '-glfo'), self.args.chain)
        else:
            print '    %s didn\'t find a germline info dir along with sw cache file, but trying to read it anyway' % utils.color('red', 'warning')

        with open(cachefname) as cachefile:
            reader = csv.DictReader(cachefile)
            for line in reader:
                utils.process_input_line(line)
                assert len(line['unique_ids']) == 1
                for region in utils.regions:  # uh... should do this more cleanly at some point
                    del line[region + '_per_gene_support']
                utils.add_implicit_info(self.glfo, line)
                if line['indelfos'][0]['reversed_seq'] != '':
                    self.info['indels'][line['unique_ids'][0]] = line['indelfos'][0]
                self.add_to_info(line)

        self.finalize(cachefname=None, just_read_cachefile=True)
        print '        water time: %.1f' % (time.time()-start)

    # ----------------------------------------------------------------------------------------
    def finalize(self, cachefname=None, just_read_cachefile=False):
        print '      info for %d / %d = %.3f' % (len(self.info['queries']), len(self.input_info), float(len(self.info['queries'])) / len(self.input_info))
        if len(self.kept_unproductive_queries) > 0:
            print '      kept %d (%.3f) unproductive' % (len(self.kept_unproductive_queries), float(len(self.kept_unproductive_queries)) / len(self.input_info))

        if just_read_cachefile:  # it's past tense!
            pass #print ''
        else:
            if len(self.skipped_unproductive_queries) > 0:
                print '         skipped %d unproductive' % len(self.skipped_unproductive_queries)
            if len(self.remaining_queries) > 0:
                # printstr = '   %s %d missing %s' % (utils.color('red', 'warning'), len(self.remaining_queries), utils.plural_str('annotation', len(self.remaining_queries)))
                if len(self.remaining_queries) < 15:
                    print '            missing annotations: ' + ' '.join(self.remaining_queries)
            if self.debug and len(self.info['indels']) > 0:
                print '      indels: %s' % ':'.join(self.info['indels'].keys())
            assert len(self.info['queries']) + len(self.skipped_unproductive_queries) + len(self.remaining_queries) == len(self.input_info)
            if self.debug and not self.args.is_data and len(self.remaining_queries) > 0:
                print 'true annotations for remaining events:'
                for qry in self.remaining_queries:
                    utils.print_reco_event(self.glfo['seqs'], self.reco_info[qry], extra_str='      ', label='true:')

        self.info['remaining_queries'] = self.remaining_queries

        found_germline_changes = False  # set to true if either alremover or alfinder found changes to the germline info
        if self.alremover is not None:
            self.alremover.finalize(self.pcounter, self.info, debug=self.args.debug_allele_finding)
            self.info['genes-to-remove'] = self.alremover.genes_to_remove
            if len(self.info['genes-to-remove']) > 0:
                found_germline_changes = True
        if self.alfinder is not None:
            cpath = None
            if not self.args.dont_collapse_clones:
                naive_seq_list = [(q, self.info[q]['naive_seq']) for q in self.info['queries']]
                cpath = ClusterPath(partition=utils.get_partition_from_collapsed_naive_seqs(naive_seq_list, self.info))
            self.alfinder.prepare_to_finalize(self.info, cpath=cpath, debug=self.args.debug_allele_finding)
            for query in self.info['queries']:
                self.alfinder.increment(self.info[query])  # it needs to know the distribution of 3p deletions before it can increment, so it has to be here
            self.alfinder.finalize(debug=self.args.debug_allele_finding)
            self.info['new-alleles'] = self.alfinder.new_allele_info
            self.info['alleles-with-evidence'] = self.alfinder.alleles_with_evidence
            if self.args.plotdir is not None:
                self.alfinder.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)
            if len(self.info['new-alleles']) > 0:
                found_germline_changes = True

        if self.args.seed_unique_id is not None:  # TODO I should really get the seed cdr3 length before running anything, and then not add seqs with different cdr3 length to start with, so those other sequences' gene matches don't get mixed in
            if self.args.seed_unique_id in self.info['queries']:
                seed_cdr3_length = self.info[self.args.seed_unique_id]['cdr3_length']
            else:  # if it failed
                print '    seed unique id \'%s\' not in final s-w queries, so removing all queries' % self.args.seed_unique_id
                seed_cdr3_length = -1
            for query in copy.deepcopy(self.info['queries']):
                if self.info[query]['cdr3_length'] != seed_cdr3_length:
                    self.remove_query(query)
            initial_n_queries = len(self.info['queries'])
            n_removed = initial_n_queries - len(self.info['queries'])
            if n_removed > 0:
                print '      removed %d / %d = %.2f sequences with cdr3 length different from seed sequence (leaving %d)' % (n_removed, initial_n_queries, float(n_removed) / initial_n_queries, len(self.info['queries']))

        # want to do this *before* we pad sequences, so that when we read the cache file we're reading unpadded sequences and can pad them below
        if cachefname is not None and not found_germline_changes:
            if self.args.write_trimmed_and_padded_seqs_to_sw_cachefname:  # hackey workaround, shouldn't be used in general
                self.pad_seqs_to_same_length()

            print '        writing sw results to %s' % cachefname
            glutils.write_glfo(cachefname.replace('.csv', '-glfo'), self.glfo)
            with open(cachefname, 'w') as outfile:
                writer = csv.DictWriter(outfile, utils.annotation_headers + utils.sw_cache_headers)
                writer.writeheader()
                missing_input_keys = set(self.input_info.keys())  # all the keys we originially read from the file
                for query in self.info['queries']:
                    missing_input_keys.remove(query)
                    outline = utils.get_line_for_output(self.info[query])  # convert lists to colon-separated strings and whatnot (doens't modify input dictionary)
                    outline = {k : v for k, v in outline.items() if k in utils.annotation_headers + utils.sw_cache_headers}  # remove the columns we don't want to output
                    writer.writerow(outline)

        self.pad_seqs_to_same_length()  # NOTE this uses all the gene matches (not just the best ones), so it has to come before we call pcounter.write(), since that fcn rewrites the germlines removing genes that weren't best matches. But NOTE also that I'm not sure what but that the padding actually *needs* all matches (rather than just all *best* matches)

        if self.perfplotter is not None:
            self.perfplotter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)

        if self.pcounter is not None:
            self.info['mute-freqs'] = {rstr : self.pcounter.mfreqer.mean_rates[rstr].get_mean() for rstr in ['all', ] + utils.regions}
            if self.parameter_out_dir is not None and not found_germline_changes:
                if self.args.plotdir is not None:
                    self.pcounter.plot(self.args.plotdir + '/sw', only_csv=self.args.only_csv_plots)
                    if self.true_pcounter is not None:
                        self.true_pcounter.plot(self.args.plotdir + '/sw-true', only_csv=self.args.only_csv_plots)
                self.pcounter.write(self.parameter_out_dir)
                if self.true_pcounter is not None:
                    self.true_pcounter.write(self.parameter_out_dir + '-true')

        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def subworkdir(self, iproc, n_procs):
        if n_procs == 1:
            return self.args.workdir
        else:
            return self.args.workdir + '/sw-' + str(iproc)

    # ----------------------------------------------------------------------------------------
    def execute_commands(self, base_infname, base_outfname, n_procs):

        def get_cmd_str(iproc):
            return self.get_ig_sw_cmd_str(self.subworkdir(iproc, n_procs), base_infname, base_outfname, n_procs)
            # return self.get_vdjalign_cmd_str(self.subworkdir(iproc, n_procs), base_infname, base_outfname, n_procs)

        cmdfos = [{'cmd_str' : get_cmd_str(iproc),
                   'workdir' : self.subworkdir(iproc, n_procs),
                   'outfname' : self.subworkdir(iproc, n_procs) + '/' + base_outfname}
                  for iproc in range(n_procs)]
        utils.run_cmds(cmdfos, batch_system=self.args.batch_system, batch_options=self.args.batch_options)

        for iproc in range(n_procs):
            os.remove(self.subworkdir(iproc, n_procs) + '/' + base_infname)
        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def write_vdjalign_input(self, base_infname, n_procs):
        queries_per_proc = float(len(self.remaining_queries)) / n_procs
        n_queries_per_proc = int(math.ceil(queries_per_proc))
        written_queries = set()  # make sure we actually write each query NOTE I should be able to remove this when I work out where they're disappearing to. But they don't seem to be disappearing any more, *sigh*
        if n_procs == 1:  # double check for rounding problems or whatnot
            assert n_queries_per_proc == len(self.remaining_queries)
        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            if n_procs > 1:
                utils.prep_dir(workdir)
            with opener('w')(workdir + '/' + base_infname) as sub_infile:
                iquery = 0
                for query_name in self.remaining_queries:  # NOTE this is wasteful to loop of all the remaining queries for each process... but maybe not that wasteful
                    if iquery >= len(self.remaining_queries):
                        break
                    if iquery < iproc*n_queries_per_proc or iquery >= (iproc + 1)*n_queries_per_proc:  # not for this process
                        iquery += 1
                        continue
                    sub_infile.write('>' + query_name + ' NUKES\n')

                    assert len(self.input_info[query_name]['seqs']) == 1  # sw can't handle multiple simultaneous sequences, but it's nice to have the same headers/keys everywhere, so we use the plural versions (with lists) even here
                    seq = self.input_info[query_name]['seqs'][0]
                    if query_name in self.info['indels']:
                        seq = self.info['indels'][query_name]['reversed_seq']  # use the query sequence with shm insertions and deletions reversed
                    sub_infile.write(seq + '\n')
                    written_queries.add(query_name)
                    iquery += 1
        not_written = self.remaining_queries - written_queries
        if len(not_written) > 0:
            raise Exception('didn\'t write %s to %s' % (':'.join(not_written), self.args.workdir))

    # ----------------------------------------------------------------------------------------
    def get_vdjalign_cmd_str(self, workdir, base_infname, base_outfname, n_procs=None):
        """
        Run smith-waterman alignment (from Connor's ighutils package) on the seqs in <base_infname>, and toss all the top matches into <base_outfname>.
        """
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        cmd_str = os.getenv('HOME') + '/.local/bin/vdjalign align-fastq -q'
        cmd_str += ' --locus ' + 'IG' + self.args.chain.upper()
        cmd_str += ' --max-drop 50'
        match, mismatch = self.match_mismatch
        cmd_str += ' --match ' + str(match) + ' --mismatch ' + str(mismatch)
        cmd_str += ' --gap-open ' + str(self.gap_open_penalty)
        cmd_str += ' --vdj-dir ' + self.my_gldir + '/' + self.args.chain
        cmd_str += ' --samtools-dir ' + self.args.partis_dir + '/packages/samtools'
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname
        return cmd_str

    # ----------------------------------------------------------------------------------------
    def get_ig_sw_cmd_str(self, workdir, base_infname, base_outfname, n_procs=None):
        """
        Run smith-waterman alignment (from Connor's ighutils package) on the seqs in <base_infname>, and toss all the top matches into <base_outfname>.
        """
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        cmd_str = self.args.ig_sw_binary
        cmd_str += ' -l ' + 'IG' + self.args.chain.upper()  # locus
        cmd_str += ' -d 50'  # max drop
        match, mismatch = self.match_mismatch
        cmd_str += ' -m ' + str(match) + ' -u ' + str(mismatch)
        cmd_str += ' -o ' + str(self.gap_open_penalty)
        cmd_str += ' -p ' + self.my_gldir + '/' + self.args.chain + '/'  # NOTE needs the trailing slash
        cmd_str += ' ' + workdir + '/' + base_infname + ' ' + workdir + '/' + base_outfname
        return cmd_str

    # ----------------------------------------------------------------------------------------
    def read_output(self, base_outfname, n_procs=1):
        queries_to_rerun = OrderedDict()  # This is to keep track of every query that we don't add to self.info (i.e. it does *not* include unproductive queries that we ignore/skip entirely because we were told to by a command line argument)
                                          # ...whereas <self.skipped_unproductive_queries> is to keep track of the queries that were definitively unproductive (i.e. we removed them from self.remaining_queries) when we were told to skip unproductives by a command line argument
        for reason in ['unproductive', 'no-match', 'weird-annot.', 'nonsense-bounds', 'invalid-codon', 'indel-fails', 'super-high-mutation']:
            queries_to_rerun[reason] = set()

        self.new_indels = 0
        queries_read_from_file = set()  # should be able to remove this, eventually
        for iproc in range(n_procs):
            outfname = self.subworkdir(iproc, n_procs) + '/' + base_outfname
            with contextlib.closing(pysam.Samfile(outfname)) as sam:  # changed bam to sam because ig-sw outputs sam files
                grouped = itertools.groupby(iter(sam), operator.attrgetter('qname'))
                for _, reads in grouped:  # loop over query sequences
                    try:
                        readlist = list(reads)
                    except:
                        print 'failed!'
                        # print 'len', len(readlist)
                        for thing in reads:
                            print thing
                        assert False
                    qinfo = self.read_query(sam.references, readlist)
                    self.summarize_query(qinfo, queries_to_rerun)  # returns before adding to <self.info> if it thinks we should rerun the query
                    queries_read_from_file.add(qinfo['name'])

        not_read = self.remaining_queries - queries_read_from_file
        if len(not_read) > 0:
            raise Exception('didn\'t read %s from %s' % (':'.join(not_read), self.args.workdir))

        if self.nth_try == 1:
            print '        processed       remaining      new-indels          rerun: ' + '      '.join([reason for reason in queries_to_rerun])
        print '      %8d' % len(queries_read_from_file),
        if len(self.remaining_queries) > 0:
            printstr = '       %8d' % len(self.remaining_queries)
            printstr += '       %8d' % self.new_indels
            printstr += '            '
            n_to_rerun = 0
            for reason in queries_to_rerun:
                printstr += '        %8d' % len(queries_to_rerun[reason])
                n_to_rerun += len(queries_to_rerun[reason])
            print printstr,
            if n_to_rerun + self.new_indels != len(self.remaining_queries):
                print ''
                raise Exception('numbers don\'t add up in sw output reader (n_to_rerun + new_indels != remaining_queries): %d + %d != %d   (look in %s)' % (n_to_rerun, self.new_indels, len(self.remaining_queries), self.args.workdir))
            if self.nth_try < 2 or self.new_indels == 0:  # increase the mismatch score if it's the first try, or if there's no new indels
                print '            increasing mismatch score (%d --> %d) and rerunning them' % (self.match_mismatch[1], self.match_mismatch[1] + 1)
                self.match_mismatch[1] += 1
            elif self.new_indels > 0:  # if there were some indels, rerun with the same parameters (but when the input is written the indel will be "reversed' in the sequences that's passed to ighutil)
                print '            rerunning for indels'
                self.new_indels = 0
            else:  # shouldn't get here
                assert False
        else:
            print '        all done'

        for iproc in range(n_procs):
            workdir = self.subworkdir(iproc, n_procs)
            os.remove(workdir + '/' + base_outfname)
            if n_procs > 1:  # still need the top-level workdir
                os.rmdir(workdir)

        sys.stdout.flush()

    # ----------------------------------------------------------------------------------------
    def get_indel_info(self, query_name, cigarstr, qrseq, glseq, gene):
        cigars = re.findall('[0-9][0-9]*[A-Z]', cigarstr)  # split cigar string into its parts
        cigars = [(cstr[-1], int(cstr[:-1])) for cstr in cigars]  # split each part into the code and the length

        codestr = ''
        qpos = 0  # position within query sequence
        indelfo = utils.get_empty_indel()  # replacement_seq: query seq with insertions removed and germline bases inserted at the position of deletions
        tmp_indices = []
        for code, length in cigars:
            codestr += length * code
            if code == 'I':  # advance qr seq but not gl seq
                indelfo['indels'].append({'type' : 'insertion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # insertion begins at <pos>
                tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
            elif code == 'D':  # advance qr seq but not gl seq
                indelfo['indels'].append({'type' : 'deletion', 'pos' : qpos, 'len' : length, 'seqstr' : ''})  # first deleted base is <pos> (well, first base which is in the position of the first deleted base)
                tmp_indices += [len(indelfo['indels']) - 1  for _ in range(length)]  # indel index corresponding to this position in the alignment
            else:
                tmp_indices += [None  for _ in range(length)]  # indel index corresponding to this position in the alignment
            qpos += length

        qrprintstr, glprintstr = '', ''
        iqr, igl = 0, 0
        for icode in range(len(codestr)):
            code = codestr[icode]
            if code == 'M':
                qrbase = qrseq[iqr]
                if qrbase != glseq[igl]:
                    qrbase = utils.color('red', qrbase)
                qrprintstr += qrbase
                glprintstr += glseq[igl]
                indelfo['reversed_seq'] += qrseq[iqr]  # add the base to the overall sequence with all indels reversed
            elif code == 'S':
                continue
            elif code == 'I':
                qrprintstr += utils.color('light_blue', qrseq[iqr])
                glprintstr += utils.color('light_blue', '*')
                indelfo['indels'][tmp_indices[icode]]['seqstr'] += qrseq[iqr]  # and to the sequence of just this indel
                igl -= 1
            elif code == 'D':
                qrprintstr += utils.color('light_blue', '*')
                glprintstr += utils.color('light_blue', glseq[igl])
                indelfo['reversed_seq'] += glseq[igl]  # add the base to the overall sequence with all indels reversed
                indelfo['indels'][tmp_indices[icode]]['seqstr'] += glseq[igl]  # and to the sequence of just this indel
                iqr -= 1
            else:
                raise Exception('unhandled code %s' % code)

            iqr += 1
            igl += 1

        dbg_str_list = ['      indels in %s' % query_name,
                        '          %20s %s' % (gene, glprintstr),
                        '          %20s %s' % ('query', qrprintstr)]
        for idl in indelfo['indels']:
            dbg_str_list.append('          %10s: %d bases at %d (%s)' % (idl['type'], idl['len'], idl['pos'], idl['seqstr']))
        indelfo['dbg_str'] = '\n'.join(dbg_str_list)

        return indelfo

    # ----------------------------------------------------------------------------------------
    def remove_query(self, query):
        # NOTE you're iterating over a deep copy of <self.info['queries']>, right? you better be!
        del self.info[query]  # this still leaves this query's gene matches in <self.info> (there may also be other traces of it)
        self.info['queries'].remove(query)
        if query in self.info['indels']:
            del self.info['indels'][query]

    # ----------------------------------------------------------------------------------------
    def add_dummy_d_match(self, qinfo, first_v_qr_end):
        dummy_d = glutils.dummy_d_genes[self.args.chain]
        qinfo['matches']['d'].append((1, dummy_d))
        qinfo['qrbounds'][dummy_d] = (first_v_qr_end, first_v_qr_end)
        qinfo['glbounds'][dummy_d] = (1, 1)

    # ----------------------------------------------------------------------------------------
    def read_query(self, references, reads):
        """ convert bam crap to python dict """
        primary = next((r for r in reads if not r.is_secondary), None)
        qinfo = {
            'name' : primary.qname,
            'seq' : primary.seq,
            'matches' : {r : [] for r in utils.regions},
            'qrbounds' : {},
            'glbounds' : {},
            'new_indels' : {}
        }

        last_scores = {r : None for r in utils.regions}
        for read in reads:  # loop over the matches found for each query sequence
            read.seq = qinfo['seq']  # only the first one has read.seq set by default, so we need to set the rest by hand
            gene = references[read.tid]
            region = utils.get_region(gene)
            qrbounds = (read.qstart, read.qend)
            glbounds = (read.pos, read.aend)
            score = read.tags[0][1]
            if last_scores[region] is not None and score > last_scores[region]:
                raise Exception('[sb]am file from smith-waterman not ordered by match score')
            last_scores[region] = score

            # NOTE it is very important not to ever skip the best match for a region, since we need its qrbounds later to know how much was trimmed

            if len(qinfo['matches'][region]) >= self.args.n_max_per_region[utils.regions.index(region)]:
                assert len(qinfo['matches'][region]) == self.args.n_max_per_region[utils.regions.index(region)]  # there better not be a way to get more than we asked for
                continue

            if 'I' in read.cigarstring or 'D' in read.cigarstring:  # shm indels!
                if len(qinfo['matches'][region]) > 0:  # skip any gene matches with indels after the first one for each region (if we want to handle [i.e. reverse] an indel, we will have stored the indel info for the first match, and we'll be rerunning)
                    continue
                assert region not in qinfo['new_indels']  # only to double-check the continue just above
                qinfo['new_indels'][region] = self.get_indel_info(qinfo['name'], read.cigarstring, qinfo['seq'][qrbounds[0] : qrbounds[1]], self.glfo['seqs'][region][gene][glbounds[0] : glbounds[1]], gene)
                qinfo['new_indels'][region]['reversed_seq'] = qinfo['seq'][ : qrbounds[0]] + qinfo['new_indels'][region]['reversed_seq'] + qinfo['seq'][qrbounds[1] : ]

            # and finally add this match's information
            qinfo['matches'][region].append((score, gene))  # NOTE it is important that this is ordered such that the best match is first
            qinfo['qrbounds'][gene] = qrbounds
            qinfo['glbounds'][gene] = glbounds

        if self.args.chain != 'h':
            _, first_v_match = qinfo['matches']['v'][0]
            self.add_dummy_d_match(qinfo, first_v_qr_end=qinfo['qrbounds'][first_v_match][1])

        return qinfo

    # ----------------------------------------------------------------------------------------
    def print_match(self, region, gene, score, qseq, glbounds, qrbounds, skipping=False):
        out_str_list = []
        buff_str = (19 - len(gene)) * ' '
        out_str_list.append('%8s%s%s%9s%3s %6.0f        ' % (' ', utils.color_gene(gene), '', '', buff_str, score))
        out_str_list.append('%4d%4d   %s\n' % (glbounds[0], glbounds[1], self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]]))
        out_str_list.append('%46s  %4d%4d' % ('', qrbounds[0], qrbounds[1]))
        out_str_list.append('   %s ' % (utils.color_mutants(self.glfo['seqs'][region][gene][glbounds[0]:glbounds[1]], qseq[qrbounds[0]:qrbounds[1]])))
        if skipping:
            out_str_list.append('skipping!')

        print ''.join(out_str_list)

    # ----------------------------------------------------------------------------------------
    def super_high_mutation(self, qinfo, best):
        gl_seqs, qr_seqs = [], []
        for region in ['v', 'j']:  # insertion mutation rates are meaningless, and we expect D to frequently be really high
            gl_bounds = qinfo['glbounds'][best[region]]
            gl_seqs.append(self.glfo['seqs'][region][best[region]][gl_bounds[0] : gl_bounds[1]])
            qr_bounds = qinfo['qrbounds'][best[region]]
            qr_seqs.append(qinfo['seq'][qr_bounds[0] : qr_bounds[1]])

        vj_mute_freq = utils.hamming_fraction(''.join(gl_seqs), ''.join(qr_seqs))
        if vj_mute_freq > self.args.max_vj_mut_freq:
            if self.debug:
                print '      super high vj mutation (%.3f > %.3f)' % (vj_mute_freq, self.args.max_vj_mut_freq)
            return True
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def get_overlap_and_available_space(self, rpair, best, qrbounds):
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]
        overlap = qrbounds[l_gene][1] - qrbounds[r_gene][0]
        available_space = qrbounds[r_gene][1] - qrbounds[l_gene][0]
        return overlap, available_space

    # ----------------------------------------------------------------------------------------
    def check_boundaries(self, rpair, qinfo, best, recursed=False, debug=False):
        # NOTE this duplicates code in shift_overlapping_boundaries(), which makes me cranky, but this setup avoids other things I dislike more
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qinfo['qrbounds'])

        if debug:
            print '  %s %s    overlap %d    available space %d' % (l_reg, r_reg, overlap, available_space)

        status = 'ok'
        if overlap > 0:  # positive overlap means they actually overlap
            status = 'overlap'
        if overlap > available_space or overlap == 1 and available_space == 1:  # call it nonsense if the boundaries are really whack (i.e. there isn't enough space to resolve the overlap) -- we'll presumably either toss the query or rerun with different match/mismatch
            status = 'nonsense'

        if debug:
            print '  overlap status: %s' % status

        if not recursed and status == 'nonsense' and l_reg == 'd' and self.nth_try > 2:  # on rare occasions with very high mutation, vdjalign refuses to give us a j match that's at all to the right of the d match
            assert l_reg == 'd' and r_reg == 'j'
            if debug:
                print '  %s: synthesizing d match' % qinfo['name']
            leftmost_position = min(qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][r_gene][0])
            qinfo['qrbounds'][l_gene] = (leftmost_position, leftmost_position + 1)  # swap whatever crummy nonsense d match we have now for a one-base match at the left end of things (things in practice should be left end of j match)
            qinfo['glbounds'][l_gene] = (0, 1)
            status = self.check_boundaries(rpair, qinfo, best, recursed=True, debug=debug)
            if status == 'overlap':
                if debug:
                    print '  \'overlap\' status after synthesizing d match. Setting to \'nonsense\', I can\'t deal with this bullshit'
                status = 'nonsense'

        return status

    # ----------------------------------------------------------------------------------------
    def shift_overlapping_boundaries(self, rpair, qinfo, best, debug=False):
        """
        s-w allows d and j matches (and v and d matches) to overlap... which makes no sense, so apportion the disputed territory between the two regions.
        Note that this still works if, say, v is the entire sequence, i.e. one match is entirely subsumed by another.
        """
        l_reg = rpair['left']
        r_reg = rpair['right']
        l_gene = best[l_reg]
        r_gene = best[r_reg]

        overlap, available_space = self.get_overlap_and_available_space(rpair, best, qinfo['qrbounds'])

        if overlap <= 0:  # nothing to do, they're already consistent
            print 'shouldn\'t get here any more if there\'s no overlap'
            return

        if overlap > available_space:
            raise Exception('overlap %d bigger than available space %d between %s and %s for %s' % (overlap, available_space, l_reg, r_reg, qinfo['name']))

        if debug:
            print '%s%s:  %d-%d overlaps with %d-%d by %d' % (l_reg, r_reg, qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][l_gene][1], qinfo['qrbounds'][r_gene][0], qinfo['qrbounds'][r_gene][1], overlap)

        l_length = qinfo['qrbounds'][l_gene][1] - qinfo['qrbounds'][l_gene][0]  # initial length of lefthand gene match
        r_length = qinfo['qrbounds'][r_gene][1] - qinfo['qrbounds'][r_gene][0]  # and same for the righthand one
        l_portion, r_portion = 0, 0  # portion of the initial overlap that we give to each side
        if debug:
            print '    lengths        portions     '
        while l_portion + r_portion < overlap:
            if debug:
                print '  %4d %4d      %4d %4d' % (l_length, r_length, l_portion, r_portion)
            if l_length <= 1 and r_length <= 1:  # don't want to erode match (in practice it'll be the d match) all the way to zero
                raise Exception('both lengths went to one without resolving overlap for %s: %s %s' % (qinfo['name'], qinfo['qrbounds'][l_gene], qinfo['qrbounds'][r_gene]))
            elif l_length > 1 and r_length > 1:  # if both have length left, alternate back and forth
                if (l_portion + r_portion) % 2 == 0:
                    l_portion += 1  # give one base to the left
                    l_length -= 1
                else:
                    r_portion += 1  # and one to the right
                    r_length -= 1
            elif l_length > 1:
                l_portion += 1
                l_length -= 1
            elif r_length > 1:
                r_portion += 1
                r_length -= 1

        if debug:
            print '  %4d %4d    %4d %4d      %s %s' % (l_length, r_length, l_portion, r_portion, '', '')
            print '      %s apportioning %d bases between %s (%d) match and %s (%d) match' % (qinfo['name'], overlap, l_reg, l_portion, r_reg, r_portion)
        assert l_portion + r_portion == overlap
        qinfo['qrbounds'][l_gene] = (qinfo['qrbounds'][l_gene][0], qinfo['qrbounds'][l_gene][1] - l_portion)
        qinfo['glbounds'][l_gene] = (qinfo['glbounds'][l_gene][0], qinfo['glbounds'][l_gene][1] - l_portion)
        qinfo['qrbounds'][r_gene] = (qinfo['qrbounds'][r_gene][0] + r_portion, qinfo['qrbounds'][r_gene][1])
        qinfo['glbounds'][r_gene] = (qinfo['glbounds'][r_gene][0] + r_portion, qinfo['glbounds'][r_gene][1])

    # ----------------------------------------------------------------------------------------
    def remove_probably_spurious_deletions(self, qinfo, best, debug=False):  # remove probably-spurious v_5p and j_3p deletions
        for erosion in utils.effective_erosions:
            region = erosion[0]
            if region == 'v':
                d_len = qinfo['glbounds'][best['v']][0]  # NOTE corresponds to *best* match (but it's what we subtract from everybody's bounds)
                insertion = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
            elif region == 'j':
                d_len = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]  # NOTE corresponds to *best* match (but it's what we subtract from everybody's bounds)
                insertion = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]
            else:
                assert False

            if d_len == 0 or len(insertion) == 0:
                continue
            if len(insertion) >= d_len:
                spurious_bases = insertion[len(insertion) - d_len:]  # i.e. the last <v_5p_del> bases of the fv insertion (or first <j_3p_del> bases of the jf insertion)
            else:
                spurious_bases = insertion  # the whole damn thing
            if spurious_bases.count(utils.ambiguous_bases[0]) == len(spurious_bases):  # don't do it if it's all Ns
                continue
            if debug:
                print 'EXPANDING %s %s d_len: %d   insertion: %s (len %d)   spurious: %s (len %d)' % (qinfo['name'], region, d_len, insertion, len(insertion), spurious_bases, len(spurious_bases))
            for gene in [g for g in qinfo['glbounds'] if utils.get_region(g) == region]:  # it's ok to assume glbounds and qrbounds have the same keys
                if region == 'v':
                    qinfo['glbounds'][gene] = (qinfo['glbounds'][gene][0] - len(spurious_bases), qinfo['glbounds'][gene][1])
                    qinfo['qrbounds'][gene] = (qinfo['qrbounds'][gene][0] - len(spurious_bases), qinfo['qrbounds'][gene][1])
                elif region == 'j':
                    qinfo['glbounds'][gene] = (qinfo['glbounds'][gene][0], qinfo['glbounds'][gene][1] + len(spurious_bases))
                    qinfo['qrbounds'][gene] = (qinfo['qrbounds'][gene][0], qinfo['qrbounds'][gene][1] + len(spurious_bases))
                else:
                    assert False

    # ----------------------------------------------------------------------------------------
    def convert_qinfo(self, qinfo, best, codon_positions):
        """ convert <qinfo> (which is from reading sam files) to format for <self.info> (this is so add_to_info() can be used by the cache file reader, as well) """
        qname = qinfo['name']
        assert qname not in self.info

        infoline = {}
        infoline['unique_ids'] = [qname, ]  # redundant, but used somewhere down the line
        infoline['seqs'] = [qinfo['seq'], ]  # NOTE this is the seq output by vdjalign, i.e. if we reversed any indels it is the reversed sequence, also NOTE many, many things depend on this list being of length one

        # erosion, insertion, mutation info for best match
        infoline['v_5p_del'] = qinfo['glbounds'][best['v']][0]
        infoline['v_3p_del'] = len(self.glfo['seqs']['v'][best['v']]) - qinfo['glbounds'][best['v']][1]  # len(germline v) - gl_match_end
        infoline['d_5p_del'] = qinfo['glbounds'][best['d']][0]
        infoline['d_3p_del'] = len(self.glfo['seqs']['d'][best['d']]) - qinfo['glbounds'][best['d']][1]
        infoline['j_5p_del'] = qinfo['glbounds'][best['j']][0]
        infoline['j_3p_del'] = len(self.glfo['seqs']['j'][best['j']]) - qinfo['glbounds'][best['j']][1]

        infoline['fv_insertion'] = qinfo['seq'][ : qinfo['qrbounds'][best['v']][0]]
        infoline['vd_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['v']][1] : qinfo['qrbounds'][best['d']][0]]
        infoline['dj_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['d']][1] : qinfo['qrbounds'][best['j']][0]]
        infoline['jf_insertion'] = qinfo['seq'][qinfo['qrbounds'][best['j']][1] : ]

        infoline['indelfos'] = [self.info['indels'].get(qname, utils.get_empty_indel()), ]

        infoline['all_matches'] = {r : [g for _, g in qinfo['matches'][r]] for r in utils.regions}  # get lists with no scores, just the names (still ordered by match quality, though)
        for region in utils.regions:
            infoline[region + '_gene'] = best[region]

        infoline['cdr3_length'] = codon_positions['j'] - codon_positions['v'] + 3
        infoline['codon_positions'] = copy.deepcopy(codon_positions)

        return infoline

    # ----------------------------------------------------------------------------------------
    def check_simulation_kbounds(self, line, true_line):
        true_kbounds = {'v' : {}, 'd' : {}}
        true_kbounds['v']['best'] = true_line['regional_bounds']['v'][1]
        true_kbounds['d']['best'] = true_line['regional_bounds']['d'][1] - true_line['regional_bounds']['v'][1]

        def print_kbound_warning():
            print '  %s true kset (%s) not within kbounds (%s) for %s' % (utils.color('red', 'warning'), utils.kbound_str(true_kbounds), utils.kbound_str({r : line['k_' + r] for r in true_kbounds}), ':'.join(line['unique_ids']))

        for region in true_kbounds:
            if true_kbounds[region]['best'] < line['k_' + region]['min'] or true_kbounds[region]['best'] >= line['k_' + region]['max']:
                print_kbound_warning()

    # ----------------------------------------------------------------------------------------
    def add_to_info(self, infoline):
        assert len(infoline['unique_ids'])
        qname = infoline['unique_ids'][0]

        self.info['queries'].append(qname)
        self.info[qname] = infoline

        # add this query's matches into the overall gene match sets
        for region in utils.regions:
            self.info['all_best_matches'].add(self.info[qname][region + '_gene'])
            self.info['all_matches'][region] |= set(self.info[qname]['all_matches'][region])  # NOTE there's an 'all_matches' in this query's info, and also in <self.info>

        # everything this is flagging seems to be cases where the insertion (or the overlapping germline region) is the same as the germline base, i.e. the kbound getter is working
        # if not self.args.is_data:
        #     self.check_simulation_kbounds(self.info[qname], self.reco_info[qname])

        if self.debug:
            inf_label = ''
            if not self.args.is_data:
                inf_label = 'inferred:'
                utils.print_reco_event(self.glfo['seqs'], self.reco_info[qname], extra_str='      ', label='true:')
            utils.print_reco_event(self.glfo['seqs'], self.info[qname], extra_str='      ', label=inf_label)

        if self.pcounter is not None:
            self.pcounter.increment(self.info[qname])
            if self.true_pcounter is not None:
                self.true_pcounter.increment(self.reco_info[qname])
        if self.perfplotter is not None:
            if qname in self.info['indels']:
                print '    skipping performance evaluation of %s because of indels' % qname  # I just have no idea how to handle naive hamming fraction when there's indels
            else:
                self.perfplotter.evaluate(self.reco_info[qname], self.info[qname])

        if not utils.is_functional(self.info[qname]):
            self.kept_unproductive_queries.add(qname)
        self.remaining_queries.remove(qname)

    # ----------------------------------------------------------------------------------------
    def summarize_query(self, qinfo, queries_to_rerun):
        """
        Fiddle with a few things, but mostly decide whether we're satisfied with the current matches.
        If not, return without calling add_to_info().
        """
        qname = qinfo['name']
        qseq = qinfo['seq']
        assert qname not in self.info

        # do we have a match for each region?
        for region in utils.getregions(self.args.chain):
            if len(qinfo['matches'][region]) == 0:
                if self.debug:
                    print '      no', region, 'match found for', qname  # TODO if no d match found, we should really just assume entire d was eroded
                queries_to_rerun['no-match'].add(qname)
                return

        best = {r : qinfo['matches'][r][0][1] for r in utils.regions}  # already made sure there's at least one match for each region

        if len(qinfo['new_indels']) > 0:  # if any of the best matches had new indels this time through (in practice: only v or j best matches)
            if self.nth_try < 2:
                if self.debug:
                    print 's-w called some indels, but we don\'t allow them first try (rerun with different match/mismatch)'
                queries_to_rerun['indel-fails'].add(qname)
                return

            if qname in self.info['indels']:
                if self.debug:
                    print '      s-w called indels after already calling some before -- but we don\'t allow multiple cycles'
                queries_to_rerun['indel-fails'].add(qname)
                return

            # NOTE (probably) important to look for v indels first (since a.t.m. we only take the first one)
            for region in [r for r in ['v', 'j', 'd'] if r in qinfo['new_indels']]:  # TODO this doesn't allow indels in more than one region
                self.info['indels'][qinfo['name']] = qinfo['new_indels'][region]  # the next time through, when we're writing ig-sw input, we look to see if each query is in <self.info['indels']>, and if it is we pass ig-sw the reversed sequence
                self.info['indels'][qinfo['name']]['reversed_seq'] = qinfo['new_indels'][region]['reversed_seq']
                self.new_indels += 1
                if self.debug:
                    print self.info['indels'][qinfo['name']]['dbg_str']
                return

            print '%s fell through indel block for %s' % (utils.color('red', 'warning'), qname)
            return  # otherwise something's weird/wrong, so we want to re-run (I think we actually can't fall through to here)

        if self.debug >= 2:
            print qname
            for region in utils.regions:
                for score, gene in qinfo['matches'][region]:  # sorted by decreasing match quality
                    self.print_match(region, gene, score, qseq, qinfo['glbounds'][gene], qinfo['qrbounds'][gene], skipping=False)

        if self.super_high_mutation(qinfo, best):
            queries_to_rerun['super-high-mutation'].add(qname)
            return

        # s-w allows d and j matches to overlap, so we need to apportion the disputed bases
        for rpair in utils.region_pairs():
            overlap_status = self.check_boundaries(rpair, qinfo, best)
            if overlap_status == 'overlap':
                self.shift_overlapping_boundaries(rpair, qinfo, best)
            elif overlap_status == 'nonsense':
                queries_to_rerun['nonsense-bounds'].add(qname)
                return
            else:
                assert overlap_status == 'ok'

        # force v 5p and j 3p matches to (in most cases) go to the end of the sequence
        self.remove_probably_spurious_deletions(qinfo, best)

        # check for suspiciously bad annotations
        for rp in utils.region_pairs():
            insertion_length = qinfo['qrbounds'][best[rp['right']]][0] - qinfo['qrbounds'][best[rp['left']]][1]  # start of right match minus end of left one
            if insertion_length > self.absolute_max_insertion_length or (self.nth_try < 2 and insertion_length > self.max_insertion_length):
                if self.debug:
                    print '      suspiciously long insertion in %s, rerunning' % qname
                queries_to_rerun['weird-annot.'].add(qname)
                return

        if self.debug == 1:
            print qname

        # set and check conserved codon positions
        codon_positions = {}
        for region, codon in utils.conserved_codons[self.args.chain].items():
            # position within original germline gene, minus the position in that germline gene at which the match starts, plus the position in the query sequence at which the match starts
            pos = self.glfo[codon + '-positions'][best[region]] - qinfo['glbounds'][best[region]][0] + qinfo['qrbounds'][best[region]][0]
            if pos < 0 or pos >= len(qseq):
                if self.debug:
                    print '      invalid %s codon position (%d in seq of length %d), rerunning' % (codon, pos, len(qseq))
                queries_to_rerun['invalid-codon'].add(qname)
                return
            codon_positions[region] = pos

        # check for unproductive rearrangements
        codons_ok = utils.both_codons_ok(self.args.chain, qseq, codon_positions)
        cdr3_length = codon_positions['j'] - codon_positions['v'] + 3

        if cdr3_length < 6:  # NOTE six is also hardcoded in utils
            if self.debug:
                print '      negative cdr3 length %d' % (cdr3_length)
            queries_to_rerun['invalid-codon'].add(qname)
            return

        in_frame_cdr3 = (cdr3_length % 3 == 0)
        stop_codon = utils.is_there_a_stop_codon(qseq, fv_insertion=qseq[ : qinfo['qrbounds'][best['v']][0]], jf_insertion=qseq[qinfo['qrbounds'][best['j']][1] : ], cyst_position=codon_positions['v'])
        if not codons_ok or not in_frame_cdr3 or stop_codon:  # TODO should use the new utils.is_functional() here
            if self.debug:
                print '       unproductive rearrangement:',
                if not codons_ok:
                    print '  bad codons',
                if not in_frame_cdr3:
                    print '  out of frame cdr3',
                if stop_codon:
                    print '  stop codon'
                print ''

            if self.nth_try < 2 and (not codons_ok or not in_frame_cdr3):  # rerun with higher mismatch score (sometimes unproductiveness is the result of a really screwed up annotation rather than an actual unproductive sequence). Note that stop codons aren't really indicative of screwed up annotations, so they don't count.
                if self.debug:
                    print '            ...rerunning'
                queries_to_rerun['unproductive'].add(qname)
                return
            elif self.args.skip_unproductive:
                if self.debug:
                    print '            ...skipping'
                self.skipped_unproductive_queries.add(qname)
                self.remaining_queries.remove(qname)
                return
            else:
                pass  # this is here so you don't forget that if neither of the above is true, we fall through and add the query to self.info

        infoline = self.convert_qinfo(qinfo, best, codon_positions)
        try:
            utils.add_implicit_info(self.glfo, infoline)
        except:
            if self.debug:
                print '      implicit info adding failed for %s, rerunning' % qname
            queries_to_rerun['weird-annot.'].add(qname)
            return

        kbounds = self.get_kbounds(infoline, qinfo, best, codon_positions)  # gets the boundaries of the non-best matches from <qinfo>
        if kbounds is None:
            if self.debug:
                print '      nonsense kbounds for %s, rerunning' % qname
            queries_to_rerun['weird-annot.'].add(qname)
            return
        infoline['k_v'] = kbounds['v']
        infoline['k_d'] = kbounds['d']

        self.add_to_info(infoline)

    # ----------------------------------------------------------------------------------------
    def get_kbounds(self, line, qinfo, best, codon_positions, debug=False):
        # NOTE
        #  - k_v is index of first d/dj insert base (i.e. length of v match)
        #  - k_v + k_d is index of first j/dj insert base (i.e. length of v + vd insert + d match)

        best_k_v = line['regional_bounds']['v'][1]  # end of v match
        best_k_d = line['regional_bounds']['d'][1] - line['regional_bounds']['v'][1]  # end of d minus end of v

        k_v_min, k_v_max = best_k_v, best_k_v
        k_d_min, k_d_max = best_k_d, best_k_d

        n_matched_to_break = 3  # once we see this many consecutive unmutated bases, assume we're actually within a germline v/d/j match

        # first make sure the hmm will check for cases in which sw over-expanded v
        if debug:
            print 'k_v_min --- %d' % k_v_min
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq = line['fv_insertion'] + line['v_qr_seqs'][0]
        glseq = line['fv_insertion'] + line['v_gl_seq']
        if k_v_min > len(qrseq):
            if debug:
                print 'k_v_min too big %d %d' % (k_v_min, len(qrseq))
            return None
        icheck = k_v_min
        while icheck > codon_positions['v'] + 3:  # i.e. stop when the last v base is the last base of the cysteine
            icheck -= 1
            if debug:
                print '    check %d' % icheck
            if qrseq[icheck] == glseq[icheck]:
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_v_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_v_min to %d' % icheck
                k_v_min = icheck
                n_matched = 0


        # then make sure the hmm will check for cases in which sw over-expanded j
        if debug:
            print 'k_d_max --- %d' % k_d_max
        k_v_min_CHK = k_v_min  # make sure we don't accidentally change it
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq, glseq = line['j_qr_seqs'][0], line['j_gl_seq']
        j_start = line['regional_bounds']['j'][0]
        icheck = k_d_max  # <icheck> is what we're considering changing k_d_max to
        while k_v_min + icheck + 1 < codon_positions['j']:  # i.e. stop when the first j/dj base is the first base of the tryp (even for the smallest k_v)
            icheck += 1
            if debug:
                print '    check %d' % icheck
            if k_v_min + icheck < j_start:  # make sure we're at least as far as the start of the j (i.e. let the hmm arbitrate between d match and dj insertion)
                if debug:
                    print '      not yet to start of j (%d + %d < %d)' % (k_v_min, icheck, j_start)
                k_d_max = icheck + 1
            elif qrseq[k_v_min + icheck - j_start] == glseq[k_v_min + icheck - j_start]:
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_d_max> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_d_max to %d' % (icheck + 1)
                k_d_max = icheck + 1  # i.e. if <icheck> doesn't match, we want the first d/dj base to be the *next* one (whereas for k_v, if <icheck> didn't match, we wanted <icheck> to be the first d/vd match)
                n_matched = 0

        assert k_v_min_CHK == k_v_min

        if debug:
            print 'k_d_min --- %d' % k_d_min
        k_v_max_CHK = k_v_max  # make sure we don't accidentally change it
        n_matched = 0  # if <n_matched_to_break> bases match, assume we're definitely within the germline
        qrseq, glseq = line['d_qr_seqs'][0], line['d_gl_seq']
        d_start = line['regional_bounds']['d'][0]
        icheck = k_d_min  # <icheck> is what we're considering changing k_d_min to
        while k_v_max + icheck > d_start:
            icheck -= 1
            if debug:
                print '    check %d' % icheck
            if qrseq[k_v_max + icheck - d_start] == glseq[k_v_max + icheck - d_start]:  # note that it's <k_v_max> here, since even with the longest k_v we want to make sure to check as short of a d as we mean to
                n_matched += 1
                if debug:
                    print '      match number %d' % n_matched
                if n_matched >= n_matched_to_break:
                    if debug:
                        print '      break on %d matches' % n_matched_to_break
                    break
            else:  # set <k_d_min> to <icheck>, i.e. move first possible d/dj insert base to index <icheck>
                if debug:
                    print '      set k_d_min to %d' % icheck
                k_d_min = icheck
                n_matched = 0

        assert k_v_max_CHK == k_v_max

        # NOTE I think there isn't any reason to increase k_v_max or decrease k_d_min -- sw already pretty much expands the v and j matches as far as it can (i.e. these fuzzing algorithms are mostly trying to decide if the last region to be matched by sw (d) should've been given more bases)
        # UPDATE not so sure, maybe I should, but it's not so important

        # also take the OR of k-space for the best <self.args.n_max_per_region[ireg]> matches (note that this happens after the expanding above because that expanding can break if it gets bounds from other matches)

        # k_v
        tmpreg = 'v'
        n_v_genes = min(self.args.n_max_per_region[utils.regions.index(tmpreg)], len(qinfo['matches'][tmpreg]))  # NOTE don't need the n_max_per_region thing any more
        for igene in range(n_v_genes):
            _, gene = qinfo['matches'][tmpreg][igene]
            this_k_v = qinfo['qrbounds'][gene][1]  # NOTE even if the v match doesn't start at the left hand edge of the query sequence, we still measure k_v from there. In other words, sw doesn't tell the hmm about it
            k_v_min = min(this_k_v, k_v_min)
            k_v_max = max(this_k_v, k_v_max)
            if debug:
                print '    %s %d %d' % (utils.color_gene(gene, width=15), k_v_min, k_v_max)

        # k_d
        tmpreg = 'd'
        n_d_genes = min(self.args.n_max_per_region[utils.regions.index(tmpreg)], len(qinfo['matches'][tmpreg]))  # NOTE don't need the n_max_per_region thing any more
        for igene in range(n_d_genes):
            _, gene = qinfo['matches'][tmpreg][igene]
            this_k_d = qinfo['qrbounds'][gene][1] - qinfo['qrbounds'][best['v']][1]  # end of d minus end of first/best v
            k_d_min = min(this_k_d, k_d_min)
            k_d_max = max(this_k_d, k_d_max)
            if debug:
                print '    %s %d %d' % (utils.color_gene(gene, width=15), k_d_min, k_d_max)


        # ----------------------------------------------------------------------------------------
        # switch to usual indexing conventions, i.e. that start with min and do not include max NOTE would be clearer to do this more coherently
        # NOTE i.e. k_[vd]_max means different things before here and after here
        # ----------------------------------------------------------------------------------------
        k_v_max += 1
        k_d_max += 1

        if self.args.chain != 'h':
            best_k_d = 1
            k_d_min = 1
            k_d_max = 2

        if best_k_v < k_v_min or best_k_v > k_v_max or best_k_d < k_d_min or best_k_d > k_d_max:
            print '  %s inconsistent best kset for %s (v: %d (%d %d)  d: %d (%d %d)' % (utils.color('red', 'error'), qinfo['name'], best_k_v, k_v_min, k_v_max, best_k_d, k_d_min, k_d_max)
            return None
        if k_v_min <= 0 or k_d_min <= 0 or k_v_min >= k_v_max or k_d_min >= k_d_max:
            print '  %s nonsense k bounds for %s (v: %d %d  d: %d %d)' % (utils.color('red', 'error'), qinfo['name'], k_v_min, k_v_max, k_d_min, k_d_max)
            return None

        kbounds = {'v' : {'best' : best_k_v, 'min' : k_v_min, 'max' : k_v_max},
                   'd' : {'best' : best_k_d, 'min' : k_d_min, 'max' : k_d_max}}
        if self.debug:
            print '    %s' % utils.kbound_str(kbounds)
        return kbounds

    # ----------------------------------------------------------------------------------------
    def remove_framework_insertions(self, debug=False):
        for query in self.info['queries']:
            swfo = self.info[query]
            assert len(swfo['seqs']) == 1
            if query in self.info['indels']:
                if self.info['indels'][query] != swfo['indelfos'][0]:
                    print '%s indelfos not equal for %s' % (utils.color('red', 'warning'), query)
                if swfo['indelfos'][0]['reversed_seq'] != swfo['seqs'][0]:
                    print '%s reversed seq not same as seq:\n%s\n%s' % (utils.color('red', 'warning'), swfo['indelfos'][0]['reversed_seq'], swfo['seqs'][0])

            # utils.print_reco_event(self.glfo['seqs'], swfo)
            utils.remove_all_implicit_info(swfo)
            fv_len = len(swfo['fv_insertion'])
            jf_len = len(swfo['jf_insertion'])

            swfo['seqs'][0] = swfo['seqs'][0][fv_len : len(swfo['seqs'][0]) - jf_len]
            if swfo['indelfos'][0]['reversed_seq'] != '':
                swfo['indelfos'][0]['reversed_seq'] = swfo['seqs'][0]
            for indel in reversed(swfo['indelfos'][0]['indels']):
                indel['pos'] -= fv_len
            for key in swfo['k_v']:
                swfo['k_v'][key] -= fv_len
            swfo['fv_insertion'] = ''
            swfo['jf_insertion'] = ''

            utils.add_implicit_info(self.glfo, swfo)
            # utils.print_reco_event(self.glfo['seqs'], swfo)

            # *sigh* not super happy about it, but I think the best way to handle this is to also remove these bases from the simulation info
            if self.reco_info is not None:
                raise Exception('needs fixing (and maybe actually shouldn\'t be fixed)')
                simfo = self.reco_info[query]
                utils.remove_all_implicit_info(simfo)
                simfo['seqs'][0] = simfo['seqs'][0][fv_len : len(simfo['seqs'][0]) - jf_len]
                if simfo['indelfos'][0]['reversed_seq'] != '':
                    simfo['indelfos'][0]['reversed_seq'] = simfo['seqs'][0]
                for indel in reversed(simfo['indelfos'][0]['indels']):
                    indel['pos'] -= fv_len
                simfo['fv_insertion'] = ''
                simfo['jf_insertion'] = ''
                utils.add_implicit_info(self.glfo, simfo)

    # ----------------------------------------------------------------------------------------
    def remove_duplicate_sequences(self, debug=False):
        uids_to_pre_keep = set()  # add these uids/seqs before looping through all the queries
        if self.args.seed_unique_id is not None:
            uids_to_pre_keep.add(self.args.seed_unique_id)
        if self.args.queries is not None:
            uids_to_pre_keep |= set(self.args.queries)

        seqs_to_keep = set()
        for utpk in uids_to_pre_keep:
            if utpk not in self.info:
                print 'requested uid %s not in sw info (probably failed above)' % utpk
                continue
            assert len(self.info[utpk]['seqs']) == 1
            seqs_to_keep.add(self.info[utpk]['seqs'][0])

        n_kept, n_removed = len(uids_to_pre_keep), 0
        for query in copy.deepcopy(self.info['queries']):
            if query in uids_to_pre_keep:  # already added it
                continue
            seq = self.info[query]['seqs'][0]
            if seq in seqs_to_keep:
                self.remove_query(query)
                n_removed += 1
            else:
                seqs_to_keep.add(seq)
                n_kept += 1

        if n_removed > 0:
            print '      removed %d / %d = %.2f duplicate sequences after trimming framework insertions (leaving %d)' % (n_removed, n_removed + n_kept, n_removed / float(n_removed + n_kept), len(self.info['queries']))

    # ----------------------------------------------------------------------------------------
    def get_padding_parameters(self, debug=False):
        padnames = ['gl_cpos', 'gl_cpos_to_j_end']

        def get_empty_maxima():
            return {pn : None for pn in padnames}

        def check_set_maxima(name, val, cdr3):
            if maxima[name] is None or val > maxima[name]:
                maxima[name] = val
            if cdr3 not in per_cdr3_maxima:
                per_cdr3_maxima[cdr3] = get_empty_maxima()
            if per_cdr3_maxima[cdr3][name] is None or val > per_cdr3_maxima[cdr3][name]:
                per_cdr3_maxima[cdr3][name] = val

        maxima = get_empty_maxima()
        per_cdr3_maxima = {}
        for query in self.info['queries']:
            swfo = self.info[query]

            # find biggest cyst position among all gl matches
            fvstuff = max(0, len(swfo['fv_insertion']) - swfo['v_5p_del'])  # we always want to pad out to the entire germline sequence, so don't let this go negative
            # loop over all matches for all sequences (up to n_max_per_region), because we want bcrham to be able to compare any sequence to any other (although, could probably use all *best* matches rather than all *all* UPDATE no, I kinda think not)
            for v_match in self.info['all_matches']['v']:
                gl_cpos = self.glfo['cyst-positions'][v_match] + fvstuff
                check_set_maxima('gl_cpos', gl_cpos, swfo['cdr3_length'])

            # Since we only store j_3p_del for the best match, we can't loop over all of 'em. But j stuff doesn't vary too much, so it works ok.
            cpos = swfo['codon_positions']['v']  # cyst position in query sequence (as opposed to gl_cpos, which is in germline allele)
            jfstuff = max(0, len(swfo['jf_insertion']) - swfo['j_3p_del'])
            gl_cpos_to_j_end = len(swfo['seqs'][0]) - cpos + swfo['j_3p_del'] + jfstuff
            check_set_maxima('gl_cpos_to_j_end', gl_cpos_to_j_end, swfo['cdr3_length'])

        if debug:
            print '    maxima:',
            for k in padnames:
                print '%s %d    ' % (k, maxima[k]),
            print ''

            print '    per-cdr3 maxima:'
            print '         %s  %s  %s' % ('cdr3', padnames[0], padnames[1])
            for cdr3 in per_cdr3_maxima:
                print '         %3d' % cdr3,
                for k in padnames:
                    print '   %d    ' % (per_cdr3_maxima[cdr3][k]),
                print ''

        return maxima, per_cdr3_maxima

    # ----------------------------------------------------------------------------------------
    def pad_seqs_to_same_length(self, debug=False):
        """
        Pad all sequences in <seqinfo> to the same length to the left and right of their conserved cysteine positions.
        Next, pads all sequences further out (if necessary) such as to eliminate all v_5p and j_3p deletions.
        """

        cluster_different_cdr3_lengths = False  # if you want glomerator.cc to try to cluster different cdr3 lengths, you need to pass it *everybody* with the same N padding... but then you're padding way more than you need to on almost every sequence, which is really wasteful and sometimes confuses bcrham

        # NOTE that an additional reason not to do this in simulation is that it will screw up the purity/completeness calculation
        if not self.args.dont_remove_framework_insertions and self.args.is_data:  # don't want to do this on simulation -- it's too much trouble to keep things consistent with the simulation info
            self.remove_framework_insertions(debug=debug)
            self.remove_duplicate_sequences(debug=debug)

        maxima, per_cdr3_maxima = self.get_padding_parameters(debug=debug)

        for query in self.info['queries']:
            swfo = self.info[query]
            assert len(swfo['seqs']) == 1

            cpos = swfo['codon_positions']['v']
            if cluster_different_cdr3_lengths:
                padleft = maxima['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
                padright = maxima['gl_cpos_to_j_end'] - (len(swfo['seqs'][0]) - cpos)
            else:
                padleft = per_cdr3_maxima[swfo['cdr3_length']]['gl_cpos'] - cpos  # left padding: biggest germline cpos minus cpos in this sequence
                padright = per_cdr3_maxima[swfo['cdr3_length']]['gl_cpos_to_j_end'] - (len(swfo['seqs'][0]) - cpos)
            if padleft < 0 or padright < 0:
                raise Exception('bad padding %d %d for %s' % (padleft, padright, query))

            leftstr = padleft * utils.ambiguous_bases[0]
            rightstr = padright * utils.ambiguous_bases[0]
            swfo['fv_insertion'] = leftstr + swfo['fv_insertion']
            swfo['jf_insertion'] = swfo['jf_insertion'] + rightstr
            swfo['seqs'][0] = leftstr + swfo['seqs'][0] + rightstr
            swfo['naive_seq'] = leftstr + swfo['naive_seq'] + rightstr  # NOTE I should eventually rewrite this to remove all implicit info, then change things, then re-add implicit info (like in remove_framework_insertions)
            if query in self.info['indels']:  # also pad the reversed sequence
                self.info['indels'][query]['reversed_seq'] = leftstr + self.info['indels'][query]['reversed_seq'] + rightstr
            for key in swfo['k_v']:
                swfo['k_v'][key] += padleft
            swfo['codon_positions']['v'] += padleft
            swfo['codon_positions']['j'] += padleft
            for region in utils.regions:
                swfo['regional_bounds'][region] = tuple([rb + padleft for rb in swfo['regional_bounds'][region]])  # I kind of want to just use a list now, but a.t.m. don't much feel like changing it everywhere else
            swfo['padlefts'] = [padleft, ]
            swfo['padrights'] = [padright, ]
            utils.add_implicit_info(self.glfo, swfo)  # check to make sure we modified everything in a consistent manner

            if debug:
                print '      pad %d %d   %s' % (padleft, padright, query)

        if debug:
            for query in self.info['queries']:
                print '%20s %3d %s' % (query, self.info[query]['cdr3_length'], self.info[query]['seqs'][0])
