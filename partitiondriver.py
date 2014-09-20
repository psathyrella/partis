import time
import sys
import json
import os
import glob
import csv
import re
from subprocess import Popen, check_call, PIPE

import_start = time.time()
from utils import utils
from utils.opener import opener
from hmmwriter import HmmWriter
from clusterer import Clusterer
from waterer import Waterer
from parametercounter import ParameterCounter

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, args):
        self.args = args
        self.germline_seqs = utils.read_germlines(self.args.datadir)

        self.safety_buffer_that_I_should_not_need = 35  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this

        with opener('r')(self.args.datadir + '/v-meta.json') as json_file:  # get location of <begin> cysteine in each v region
            self.cyst_positions = json.load(json_file)

        self.precluster_info = {}
        self.input_info = {}
        self.reco_info = {}  # generator truth information

    # ----------------------------------------------------------------------------------------
    def from_same_event(self, query_name, second_query_name):
        if self.args.pair:
            return self.reco_info[query_name]['reco_id'] == self.reco_info[second_query_name]['reco_id']
        else:
            return False

    # ----------------------------------------------------------------------------------------
    def read_input_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        with opener('r')(self.args.seqfile) as seqfile:
            delimiter = ','
            name_column = 'unique_id'
            seq_column = 'seq'
            # n_nukes_skip = 0  # skip all the TTTTTTTTs in data. TODO er that seems hackey doesn't it?
            if '.tsv' in self.args.seqfile:
                delimiter = '\t'
                name_column = 'name'
                seq_column = 'nucleotide'
                # n_nukes_skip = 9  # TODO ehhhhh I don't like this
            reader = csv.DictReader(seqfile, delimiter=delimiter)
            n_queries = 0
            for line in reader:
                # if command line specified query or reco ids, skip other ones
                if self.args.queries != None and line[name_column] not in self.args.queries:
                    continue
                if self.args.reco_ids != None and line['reco_id'] not in self.args.reco_ids:
                    continue

                self.input_info[line[name_column]] = {'unique_id':line[name_column], 'seq':line[seq_column]}  #[n_nukes_skip:]}
                if not self.args.is_data:
                    self.reco_info[line['unique_id']] = line
                n_queries += 1
                if self.args.n_total_queries > 0 and n_queries >= self.args.n_total_queries:
                    break
        assert len(self.input_info) > 0
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        sys.exit()
        self.read_input_file()  # read simulation info and write sw input file

        # run smith-waterman
        from_scratch, write_parameters = True, True
        sw_parameter_dir = ''
        if self.args.cache_parameters:
            from_scratch = True
            write_parameters = True
            sw_parameter_dir = self.args.parameter_dir
        elif self.args.read_cached_parameters:
            from_scratch = False
            write_parameters = False
            sw_parameter_dir = self.args.parameter_dir
        else:
            from_scratch = True
            write_parameters = True
            sw_parameter_dir = self.args.workdir + '/sw_parameters'
        sw_plotdir = os.getenv('www') + '/partis/sw_parameters'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs,
                          from_scratch=from_scratch, parameter_dir=sw_parameter_dir, write_parameters=write_parameters, plotdir=sw_plotdir)
        waterer.run()

        # cdr3 length partitioning
        cdr3_cluster = False  # don't precluster on cdr3 length for the moment -- I cannot accurately infer cdr3 length in some sequences, so I need a way to pass query seqs to the clusterer with several possible cdr3 lengths. TODO fix that
        cdr3_length_clusters = None
        if cdr3_cluster:
            cdr3_length_clusters = self.cdr3_length_precluster(waterer)

        parameter_dir = sw_parameter_dir  # NOTE this notation kinda sucks, but I'm just making the point that if we've read cached parameters they're not necessarily from smith-waterman. TODO unhackify it
        if not self.args.read_cached_parameters:
            print 'writing hmm files'
            if self.args.only_genes != None:
                self.write_specified_hmms(sw_parameter_dir, self.args.only_genes)
            else:
                self.write_all_hmms(sw_parameter_dir)

        if self.args.cache_parameters:
            print 'cached the parameters, now exiting'
            sys.exit()

        if self.args.pair and self.args.algorithm == 'forward':
            assert False  # oh man this is getting complicated need to fix
            # hamming_clusters = self.hamming_precluster(cdr3_length_clusters)
            # stripped_clusters = self.run_hmm(waterer, parameter_dir=parameter_dir, preclusters=hamming_clusters, stripped=True)
            # final_clusters = self.run_hmm(waterer, parameter_dir=parameter_dir, preclusters=stripped_clusters, stripped=False)
        else:
            self.run_hmm(waterer, parameter_in_dir=parameter_dir, count_parameters=True, plotdir=os.getenv('www') + '/partis/hmm_parameters')

    # ----------------------------------------------------------------------------------------
    def clean(self, waterer):
        if self.args.no_clean:
            return
        if not self.args.read_cached_parameters:
            waterer.clean()  # remove all the parameter files that the waterer's ParameterCounter made
            for fname in glob.glob(waterer.parameter_dir + '/hmms/*.hmm'):  # remove all the hmm files from the same directory
                os.remove(fname)
            os.rmdir(waterer.parameter_dir + '/hmms')
            os.rmdir(waterer.parameter_dir)

        # check to see if we missed anything
        leftovers = [ fname for fname in os.listdir(self.args.workdir) ]
        for over in leftovers:
            print self.args.workdir + '/' + over
        os.rmdir(self.args.workdir)

    # ----------------------------------------------------------------------------------------
    def cache_parameters(self):
        self.read_input_file()  # read simulation info and write sw input file

        write_parameters = True
        sw_parameter_dir = self.args.parameter_dir + '/sw_parameters'
        sw_plotdir = os.getenv('www') + '/partis/sw_parameters'
        waterer = Waterer(self.args, self.input_info, self.reco_info, self.germline_seqs,
                          from_scratch=True, parameter_dir=sw_parameter_dir, write_parameters=True, plotdir=sw_plotdir)
        waterer.run()

        print 'writing hmm files'
        self.write_hmms(sw_parameter_dir, self.args.only_genes)

        hmm_parameter_dir = self.args.parameter_dir + '/hmm_parameters'
        assert not self.args.pair and self.args.algorithm == 'viterbi'  # er, could probably do it for forward pair, but I'm not really sure what it would mean
        self.run_hmm(waterer.info, parameter_in_dir=sw_parameter_dir, parameter_out_dir=hmm_parameter_dir, count_parameters=True, plotdir=os.getenv('www') + '/partis/hmm_parameters')

    # ----------------------------------------------------------------------------------------
    def run_hmm(self, sw_info, parameter_in_dir, parameter_out_dir='', preclusters=None, stripped=False, prefix='', count_parameters=False, plotdir=''):
        pcounter = None
        if count_parameters:
            pcounter = ParameterCounter(self.germline_seqs, parameter_out_dir, plotdir=plotdir)
        if plotdir != '':
            utils.prep_dir(plotdir + '/plots', '*.svg')

        if prefix == '' and stripped:
            prefix = 'stripped'
        print prefix,'hmm'
        csv_infname = self.args.workdir + '/' + prefix + '_hmm_input.csv'
        csv_outfname = self.args.workdir + '/' + prefix + '_hmm_output.csv'
        pairscorefname = self.args.workdir + '/' + prefix + '_hmm_pairscores.csv'
        self.write_hmm_input(csv_infname, sw_info, preclusters=preclusters, stripped=stripped)
        self.run_stochhmm(csv_infname, csv_outfname, parameter_dir=parameter_in_dir)
        self.read_hmm_output(csv_outfname, pairscorefname, pcounter)

        if count_parameters:
            pcounter.write_counts()

        clusters = None
        if self.args.pair and self.args.algorithm == 'forward':
            clusters = Clusterer(0, greater_than=True)
            clusters.cluster(pairscorefname, debug=False, reco_info=self.reco_info)
            if preclusters != None:
                for query_name in sw_info:  # check for singletons that got split out in the preclustering step
                    if query_name not in clusters.query_clusters:
                        print '    singleton ', query_name

        if not self.args.no_clean:
            os.remove(csv_infname)
            os.remove(csv_outfname)
            os.remove(pairscorefname)

        return clusters

    # ----------------------------------------------------------------------------------------
    def cdr3_length_precluster(self, waterer, preclusters=None):
        cdr3lengthfname = self.args.workdir + '/cdr3lengths.csv'
        with opener('w')(cdr3lengthfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'cdr3_length', 'second_cdr3_length', 'score'))
            writer.writeheader()
            for query_name, second_query_name in self.get_pairs(preclusters):
                cdr3_length = waterer.info[query_name]['cdr3_length']
                second_cdr3_length = waterer.info[second_query_name]['cdr3_length']
                same_length = cdr3_length == second_cdr3_length
                if not self.args.is_data:
                    assert cdr3_length == int(self.reco_info[query_name]['cdr3_length'])
                    if second_cdr3_length != int(self.reco_info[second_query_name]['cdr3_length']):
                        print 'WARNING did not infer correct cdr3 length'
                        assert False
                writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'cdr3_length':cdr3_length, 'second_cdr3_length':second_cdr3_length, 'score':int(same_length)})

        clust = Clusterer(0.5, greater_than=True)  # i.e. cluster together if same_length == True
        clust.cluster(cdr3lengthfname, debug=False)
        os.remove(cdr3lengthfname)
        return clust

    # ----------------------------------------------------------------------------------------
    def get_pairs(self, preclusters=None):
        """ Get all unique the pairs of sequences in input_info, skipping where preclustered out """
        pair_scores = {}
        pairs = []
        n_lines, n_preclustered, n_previously_preclustered, n_removable, n_singletons = 0, 0, 0, 0, 0
        for query_name in self.input_info:
            for second_query_name in self.input_info:
                if second_query_name == query_name:
                    continue
                if utils.get_key(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                    continue
                pair_scores[utils.get_key(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                if preclusters != None:  # if we've already run preclustering, skip the pairs that we know aren't matches
                    if query_name not in preclusters.query_clusters or second_query_name not in preclusters.query_clusters:  # singletons (i.e. they were already preclustered into their own group)
                        n_singletons += 1
                        continue
                    if utils.get_key(query_name, second_query_name) not in preclusters.pairscores:  # preclustered out in a previous preclustering step
                        n_previously_preclustered += 1
                        continue
                    if preclusters.query_clusters[query_name] != preclusters.query_clusters[second_query_name]:  # not in same cluster
                        n_preclustered += 1
                        continue
                    if preclusters.is_removable(preclusters.pairscores[utils.get_key(query_name, second_query_name)]):  # in same cluster, but score (link) is long. i.e. *this* pair is far apart, but other seqs to which they are linked are close to each other
                        n_removable += 1
                        continue
                pairs.append((query_name, second_query_name))
                n_lines += 1

        print '    %d lines (%d preclustered out, %d removable links, %d singletons, %d previously preclustered)' % (n_lines, n_preclustered, n_removable, n_singletons, n_previously_preclustered)
        return pairs

    # ----------------------------------------------------------------------------------------
    def hamming_precluster(self, preclusters=None):
        start = time.time()
        print 'hamming clustering'
        hammingfname = self.args.workdir + '/fractional-hamming-scores.csv'
        with opener('w')(hammingfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'score'))
            writer.writeheader()
            for query_name, second_query_name in self.get_pairs(preclusters):
                query_seq = self.input_info[query_name]['seq']
                second_query_seq = self.input_info[second_query_name]['seq']
                mutation_frac = utils.hamming(query_seq, second_query_seq) / float(len(query_seq))
                writer.writerow({'unique_id':query_name, 'second_unique_id':second_query_name, 'score':mutation_frac})
                if self.args.debug:
                    print '    %20s %20s %8.2f' % (query_name, second_query_name, mutation_frac)

        clust = Clusterer(0.5, greater_than=False)  # TODO this 0.5 number isn't gonna be the same if the query sequences change length
        clust.cluster(hammingfname, debug=False)
        os.remove(hammingfname)
        print '    hamming time: %.3f' % (time.time()-start)
        return clust
    
    # ----------------------------------------------------------------------------------------
    def write_hmms(self, parameter_dir, gene_list=None):
        hmm_dir = parameter_dir + '/hmms'
        utils.prep_dir(hmm_dir, '*.hmm')

        if gene_list == None:  # if specific genes weren't specified, do 'em all
            for region in utils.regions:
                gene_list.append(self.germline_seqs[region])

        for gene in gene_list:
            if self.args.debug:
                print '   ',gene
            if len(re.findall('J[123]P', gene)) > 0:  # pretty sure these versions are crap
                print '  poof'
                continue
            if 'OR16' in gene or 'OR21' in gene or 'V3/OR15' in gene:
                print '  poof'
                continue
            writer = HmmWriter(parameter_dir, hmm_dir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=self.args.v_right_length)
            writer.write()

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, csv_fname, sw_info, preclusters=None, stripped=False):  # TODO use different input files for the two hmm steps
        with opener('w')(csv_fname) as csvfile:
            # write header
            header = ['name', 'second_name', 'k_v_min', 'k_v_max', 'k_d_min', 'k_d_max', 'only_genes', 'seq', 'second_seq']  # I wish I had a good c++ csv reader 
            csvfile.write(' '.join(header) + '\n')

            # then write a line for each query sequence (or pair of them)
            if self.args.pair:
                for query_name, second_query_name in self.get_pairs(preclusters):  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
                    info = sw_info[query_name]
                    second_info = sw_info[second_query_name]
    
                    # make sure we don't need to rewrite the hmm model files
                    if abs(info['v_right_length'] - self.args.v_right_length) >= self.safety_buffer_that_I_should_not_need:
                        print 'WARNING VRIGHT ', info['v_right_length'], self.args.v_right_length, (info['v_right_length']-self.args.v_right_length)
                    # assert abs(info['v_right_length'] - self.default_v_right_length) < self.safety_buffer_that_I_should_not_need
    
                    # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                    assert len(re.findall('J[123]P', info['best']['j'])) == 0
    
                    # TODO the current method of expanding the k space bounds to encompass all the gene matches, and both sequences, is *correct*, but I think it's unnecessarily slow
                    k_v, k_d = {}, {}
                    k_v['min'] = min(info['k_v']['min'], second_info['k_v']['min'])
                    k_v['max'] = max(info['k_v']['max'], second_info['k_v']['max'])
                    k_d['min'] = min(info['k_d']['min'], second_info['k_d']['min'])
                    k_d['max'] = max(info['k_d']['max'], second_info['k_d']['max'])

                    only_genes = info['all']
                    second_only_genes = second_info['all']
                    if stripped:  # strip down the hmm -- only use the single best gene for each sequence, and don't fuzz at all
                        only_genes = info['best']['v'] + ':' + info['best']['d'] + ':' + info['best']['j']
                        second_only_genes = second_info['best']['v'] + ':' + second_info['best']['d'] + ':' + second_info['best']['j']
                        k_v['min'] = info['k_v']['best']
                        k_v['max'] = k_v['min'] + 1
                        k_d['min'] = info['k_d']['best']
                        k_d['max'] = k_d['min'] + 1

                    only_genes = ':'.join(set(only_genes.split(':')) | set(second_only_genes.split(':')))  # NOTE using both sets of genes (from both query seqs) like this *really* helps, 
                    csvfile.write('%s %s %d %d %d %d %s %s %s\n' %
                                  (query_name, second_query_name, k_v['min'], k_v['max'], k_d['min'], k_d['max'], only_genes,
                                   self.input_info[query_name]['seq'], self.input_info[second_query_name]['seq']))
            else:
                for query_name in self.input_info:
                    info = sw_info[query_name]
    
                    # make sure we don't need to rewrite the hmm model files
                    if abs(info['v_right_length'] - self.args.v_right_length) >= self.safety_buffer_that_I_should_not_need:
                        print 'WARNING VRIGHT ', info['v_right_length'], self.args.v_right_length, (info['v_right_length']-self.args.v_right_length)
                    # assert abs(info['v_right_length'] - self.default_v_right_length) < self.safety_buffer_that_I_should_not_need
    
                    # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                    assert len(re.findall('J[123]P', info['j_gene'])) == 0
    
                    only_genes = info['all']
                    # assert self.args.algorithm == 'viterbi'  # TODO hm, was that really all I had to do to allow non-pair forward?
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, info['k_v']['min'], info['k_v']['max'], info['k_d']['min'], info['k_d']['max'], only_genes, self.input_info[query_name]['seq']))

    # ----------------------------------------------------------------------------------------
    def run_stochhmm(self, csv_infname, csv_outfname, parameter_dir):
        start = time.time()

        # build the command line
        cmd_str = self.args.stochhmm_dir + '/stochhmm'
        cmd_str += ' --algorithm ' + self.args.algorithm
        if self.args.pair:
            cmd_str += ' --pair 1'
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + parameter_dir + '/hmms'
        cmd_str += ' --datadir ' + self.args.datadir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        # cmd_str += ' >/dev/null'

        if self.args.debug == 2:  # not sure why, but popen thing hangs with debug 2
            check_call(cmd_str, shell=True)  # um, not sure which I want here, but it doesn't really matter. TODO kinda
            sys.exit()

        # run hmm
        hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        hmm_proc.wait()
        hmm_out, hmm_err = hmm_proc.communicate()
        print hmm_out,
        print hmm_err,
        if hmm_proc.returncode != 0:
            print 'aaarrrrrrrgggggh\n', hmm_out, hmm_err
            sys.exit()

        print '    hmm run time: %.3f' % (time.time() - start)

    # ----------------------------------------------------------------------------------------
    def read_hmm_output(self, hmm_csv_outfname, pairscorefname, pcounter):
        # TODO the input and output files for this function are almost identical at this point

        # write header for pairwise score file
        with opener('w')(pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id,second_unique_id,score\n')
        with opener('r')(hmm_csv_outfname) as hmm_csv_outfile:
            reader = csv.DictReader(hmm_csv_outfile)
            last_id = None
            n_boundary_errors = 0
            for line in reader:

                # check for errors
                boundary_error = False
                if 'boundary' in line['errors'].split(':'):
                    n_boundary_errors += 1
                    boundary_error = True
                else:
                    assert len(line['errors']) == 0

                # print the results
                if self.args.algorithm == 'viterbi':  # for viterbi print the inferred reco events
                    if last_id != utils.get_key(line['unique_id'], line['second_unique_id']):  # if this is the first line for this query (or query pair), print the true event
                        # increment counters
                        if pcounter != None:
                            utils.get_match_seqs(self.germline_seqs, line)
                            pcounter.increment(line)

                        # then print stuff
                        print '%20s %20s   %d' % (line['unique_id'], line['second_unique_id'], self.from_same_event(line['unique_id'], line['second_unique_id']))
                        print '    true:'
                        cpos = self.cyst_positions[self.reco_info[line['unique_id']]['v_gene']]['cysteine-position']
                        # TODO fix this, i.e. figure out the number to add to it so it's correct: tpos = int(self.tryp_positions[self.reco_info[line['unique_id']]['j_gene']]) + 
                        utils.print_reco_event(self.germline_seqs, self.reco_info[line['unique_id']], cpos, 0, extra_str='    ')
                        if self.args.pair:
                            utils.print_reco_event(self.germline_seqs, self.reco_info[line['second_unique_id']], 0, 0, self.from_same_event(line['unique_id'], line['second_unique_id']), '    ')
                        print '    inferred:'
                    utils.print_reco_event(self.germline_seqs, line, 0, 0, extra_str='    ')
                    if self.args.pair:
                        tmpseq = line['seq']  # temporarily set 'seq' to the second query's seq. TODO oh, man, that's a cludge
                        line['seq'] = line['second_seq']
                        utils.print_reco_event(self.germline_seqs, line, 0, 0, True, extra_str='    ')
                        line['seq'] = tmpseq
                else:  # for forward, write the pair scores to file to be read by the clusterer
                    with opener('a')(pairscorefname) as pairscorefile:
                        pairscorefile.write('%s,%s,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))

                last_id = utils.get_key(line['unique_id'], line['second_unique_id'])
        if n_boundary_errors > 0:
            print '    %d boundary errors' % n_boundary_errors
