import time
import_start = time.time()
import sys
import os
import csv
import re
import itertools
import StringIO
import operator
import pysam
import contextlib
from subprocess import Popen, check_call, PIPE
import utils
from opener import opener
from hmmwriter import HmmWriter

print 'import time: %.3f' % (time.time() - import_start)

# ----------------------------------------------------------------------------------------
class PartitionDriver(object):
    def __init__(self, datadir, args, default_v_right_length=90):
        self.datadir = datadir
        self.args = args
        self.germline_seqs = utils.read_germlines()
        self.workdir = '/tmp/' + os.getenv('USER') + '/hmms/' + str(os.getpid())  # use a tmp dir specific to this process for the hmm input file
        if not os.path.exists(self.workdir):
            os.makedirs(self.workdir)
        self.reco_info = self.read_sim_file()  # read simulation info and write sw input file
        self.default_v_fuzz = 5  # TODO play around with these default fuzzes
        self.default_d_fuzz = 10  # note that the k_d guess can be almost worthless, since the s-w step tends to expand the j match over a lot of the d
        self.safety_buffer_that_I_should_not_need = 15  # the default one is plugged into the cached hmm files, so if this v_right_length is really different, we'll have to rewrite the hmms TODO fix this
        self.pairscorefname = 'pairwise-scores.csv'
        self.default_hmm_dir = 'bcell/hmms/' + self.args.human + '/' + self.args.naivety
        self.default_v_right_length = default_v_right_length
        self.sw_info = {}
        self.precluster_info = {}
        self.match_names_for_all_queries = set()
            
    # ----------------------------------------------------------------------------------------
    def get_score_index(self, query_name, second_query_name):
        """
        Return a hashable combination of the two query names that's the same if we reverse their order.
        At the moment, just add 'em up.
        """
        assert query_name != ''
        if second_query_name == '':
            second_query_name = '0'
        return int(query_name) + int(second_query_name)

    # ----------------------------------------------------------------------------------------
    def read_sim_file(self):
        """ Read simulator info and write it to input file for sw step. Returns dict of simulation info. """
        reco_info = {}  # generator truth information
        with opener('r')(self.args.simfile) as simfile:
            last_reco_id = -1  # only run on the first seq in each reco event. they're all pretty similar
            reader = csv.DictReader(simfile)
            n_queries = 0
            for line in reader:
                reco_info[line['unique_id']] = line
                last_reco_id = line['reco_id']
                n_queries += 1
                if self.args.n_total_queries > 0 and n_queries >= self.args.n_total_queries:
                    break
        return reco_info
    
    # ----------------------------------------------------------------------------------------
    def run(self):
        # write header for pairwise score file
        with opener('w')(self.pairscorefname) as pairscorefile:
            pairscorefile.write('unique_id,second_unique_id,score\n')

        # run smith-waterman
        swinfname = self.workdir + '/seq.fa'  # file for input to s-w step
        swoutfname = swinfname.replace('.fa', '.bam')
        if not self.args.skip_sw:
            sys.stdout.flush()
            self.run_smith_waterman(swinfname, swoutfname)
            os.remove(swinfname)
        else:
            swoutfname = swoutfname.replace(os.path.dirname(swoutfname), '.')
        self.read_smith_waterman(swoutfname)  # read sw bam output, collate it, and write to csv for hmm input
        os.remove(swoutfname)

        # run hmm
        single_clusters = self.single_gene_score_precluster()  # TODO what should n_max_per_region be for this step?

        # TODO then do another precluster step with a super stripped down pair forward hmm (zero fuzz, etc)

        hmm_csv_infname = self.workdir + '/hmm_input.csv'
        self.write_hmm_input(hmm_csv_infname, single_gene_precluster=False, preclusters=single_clusters)  # TODO man those variable names suck. wtf dude?
        hmm_csv_outfname = self.workdir + '/hmm_output.csv'
        self.run_stochhmm(hmm_csv_infname, hmm_csv_outfname)
        self.read_stochhmm_output(hmm_csv_outfname);

        from clusterer import Clusterer
        clust = Clusterer(0, greater_than=True)
        clust.cluster(self.pairscorefname, debug=True)
        for query_name in self.sw_info:
            if query_name not in clust.query_clusters:
                print 'singleton ',query_name
        
        os.remove(hmm_csv_infname)
        os.remove(hmm_csv_outfname)
        os.rmdir(self.workdir)

    # ----------------------------------------------------------------------------------------
    def run_smith_waterman(self, infname, outfname):
        """
        Run smith-waterman alignment on the seqs in <infname>, and toss all the top matches into <outfname>.
        Then run through <outfname> to get the top hits and their locations to pass to the hmm.
        Then run the hmm on each gene set.
        """
        with opener('w')(infname) as swinfile:  # write *all* the input seqs to file, i.e. run s-w on all of 'em at once
            for query_name,line in self.reco_info.iteritems():
                swinfile.write('>' + query_name + ' NUKES\n')
                swinfile.write(line['seq'] + '\n')
        start = time.time()
        # large gap-opening penalty: we want *no* gaps in the middle of the alignments
        # match score larger than (negative) mismatch score: we want to *encourage* some level of shm. If they're equal, we tend to end up with short unmutated alignments, which screws everything up
        check_call('/home/dralph/.local/bin/vdjalign align-fastq --j-subset adaptive --max-drop 50 --match 3 --mismatch 1 --gap-open 100 ' + infname + ' ' + outfname + ' 2>/dev/null', shell=True)
        print 's-w time: %.3f' % (time.time()-start)
    
    # ----------------------------------------------------------------------------------------
    def read_smith_waterman(self, infname):
        """
        Read bamfile output by s-w step, and write the info (after a bit of collation) to a csv.
        Note that the only reason to bother writing the csv is so you can avoid rerunning the s-w step every time you run the hmm.
        """
        start = time.time()
        gene_choice_probs = utils.read_overall_gene_prob(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety)
        with contextlib.closing(pysam.Samfile(infname)) as bam:
            grouped = itertools.groupby(iter(bam), operator.attrgetter('qname'))
            for _, reads in grouped:  # loop over query sequences
                reads = list(reads)
                primary = next((r for r in reads if not r.is_secondary), None)
                query_seq = primary.seq
                query_name = primary.qname
                raw_best = {}
                all_match_names = {}
                for region in utils.regions:
                    all_match_names[region] = []
                all_query_bounds, all_germline_bounds = {}, {}
                for read in reads:  # loop over the matches found for each query sequence
                    read.seq = query_seq  # only the first one has read.seq set by default, so we need to set the rest by hand
                    gene = bam.references[read.tid]
                    region = utils.get_region(gene)
    
                    if region not in raw_best:  # best v, d, and j before multiplying by gene choice probs. needed 'cause *these* are the v and j that get excised
                        raw_best[region] = gene
    
                    if 'J1P' in gene or 'J3P' in gene:
                        continue
    
                    raw_score = read.tags[0][1]  # raw because they don't include the gene choice probs. TODO oh wait shit this isn't right. raw_score isn't a prob. what the hell is it, anyway?
                    choice_prob = 0.0  # set to zero if we didn't see it in data. kinda hacky, I suppose
                    if gene in gene_choice_probs[region]:
                        choice_prob = gene_choice_probs[region][gene]
                    score = choice_prob * raw_score  # multiply by the probability to choose this gene
                    all_match_names[region].append((score,gene))
                    all_query_bounds[gene] = (read.qstart, read.qend)
                    all_germline_bounds[gene] = (read.pos, read.aend)
    
                # append the best matches to match_names and work out how much v_fuzz we need
                best = {}
                match_names = {}
                n_matches = {'v':0, 'd':0, 'j':0}
                n_used = {'v':0, 'd':0, 'j':0}
                v_fuzz = self.default_v_fuzz  # TODO I have the defaults adjusted in two places at the moment
                d_fuzz = self.default_d_fuzz
                for region in utils.regions:
                    all_match_names[region] = sorted(all_match_names[region], reverse=True)
                    match_names[region] = []
                for region in utils.regions:
                    for score,gene in all_match_names[region]:
                        n_matches[region] += 1
                        if n_matches[region] > self.args.n_max_per_region:  # only take the top few from each region
                            # TODO should use *lots* of d matches, but fewer vs and js
                            # assert False  # TODO also should loop over *way* more k_d than k_v
                            continue
                        n_used[region] += 1
                        match_names[region].append(gene)
    
                        # check consistency with best match (since the best match is excised in s-w code, and because stochhmm is run with *one* k_v k_d set)
                        if region not in best:
                            best[region] = gene
                        else:
                            if region == 'v':
                                necessary_fuzz = abs(all_query_bounds[gene][1] - all_query_bounds[raw_best[region]][1])
                                if necessary_fuzz > v_fuzz-1:  # if this v match ends at a wildly different position than did the best v match, expand v_fuzz accordingly
                                    v_fuzz = necessary_fuzz + 1
    
                        if self.args.debug:
                            buff_str = (17 - len(gene)) * ' '
                            print '%8s%s%s%6.1f' % (' ', utils.color_gene(gene), buff_str, score),
                            glbounds = all_germline_bounds[gene]
                            qrbounds = all_query_bounds[gene]
                            print ' %4d%4d   %s' % (glbounds[0], glbounds[1], self.germline_seqs[region][gene][glbounds[0]:glbounds[1]])
                            print '%31s  %4d%4d   %s' % ('', qrbounds[0], qrbounds[1], utils.color_mutants(self.germline_seqs[region][gene][glbounds[0]:glbounds[1]], query_seq[qrbounds[0]:qrbounds[1]]))
                                
                # print how many of the available matches we used
                if self.args.debug:
                    print '  used',
                    for region in utils.regions:
                        if region != 'v':
                            print '      ',
                        print ' %d / %d in %s' % (n_used[region], n_matches[region], region)
    
                # write sw info to file for hmm
                k_v = all_query_bounds[best['v']][1]  # end of v match
                k_d = all_query_bounds[best['d']][1] - all_query_bounds[best['v']][1]  # end of d minus end of v
                if k_d < 5:  # since the s-w step matches to the longest possible j and then excises it, this sometimes gobbles up the d, resulting in a very short d alignment. 
                    k_d = max(5, k_d)
                    d_fuzz = max(d_fuzz, 10)
                    
                v_right_length = len(self.germline_seqs['v'][best['v']]) - all_germline_bounds[best['v']][0]  # germline v length minus (germline) start of v match
                # if 'IGHJ4*' in best_genes['j'] and self.germline_seqs['d'][best_genes['d']][-5:] == 'ACTAC':  # the end of some d versions is the same as the start of some j versions, so the s-w frequently kicks out the 'wrong' alignment
                #     d_fuzz = 10

                self.match_names_for_all_queries |= set(match_names['v'] + match_names['d'] + match_names['j'])

                assert query_name not in self.sw_info
                self.sw_info[query_name] = {}
                self.sw_info[query_name]['k_v'] = k_v
                self.sw_info[query_name]['k_d'] = k_d
                self.sw_info[query_name]['v_fuzz'] = v_fuzz
                self.sw_info[query_name]['d_fuzz'] = d_fuzz
                self.sw_info[query_name]['v_right_length'] = v_right_length
                self.sw_info[query_name]['best'] = best
                self.sw_info[query_name]['all'] = ':'.join(match_names['v'] + match_names['d'] + match_names['j'])

                assert k_v > 0 and k_d > 0 and v_fuzz > 0 and d_fuzz > 0 and v_right_length > 0
        print 'sw read time: %.3f' % (time.time() - start)
    
    # ----------------------------------------------------------------------------------------
    def single_gene_score_precluster(self):
        hmm_csv_infname = self.workdir + '/single-gene-hmm-input.csv'  # list of single query seqs and genes for stochhmm to check
        single_gene_prob_fname = self.workdir + '/single-gene-probs.csv'  # list of scores for each gene output by stochhmm
        precluster_outfname = self.workdir + '/scorespace-distances.csv'  # single gene scorespace distances for input to the clusterer

        self.write_hmm_input(hmm_csv_infname, single_gene_precluster=True)
        self.run_stochhmm(hmm_csv_infname, single_gene_prob_fname, single_gene_precluster=True)
        self.read_stochhmm_output(single_gene_prob_fname, precluster=True);
        self.calculate_scorespace_distances(precluster_outfname)
        from clusterer import Clusterer
        clust = Clusterer(0.4, greater_than=False)  # TODO will need a way to get the cut value automatically. *sigh*
        clust.cluster(precluster_outfname, debug=True)
        # for cluster_id in clust.cluster_ids:
        #     for query in iter(clust.query_clusters):
        #         if clust.query_clusters[query] == cluster_id:
        #             print '%s,%d' % (query, cluster_id)
        os.remove(single_gene_prob_fname)
        os.remove(hmm_csv_infname)
        os.remove(precluster_outfname)
        return clust

    # ----------------------------------------------------------------------------------------
    def write_specified_hmms(self, hmmdir, gene_list, v_right_length):
        for gene in gene_list:
            if len(re.findall('J[123]P', gene)) > 0:  # pretty sure these versions are crap
                print '  poof'
                continue
            if 'OR16' in gene or 'OR21' in gene or 'V3/OR15' in gene:
                print '  poof'
                continue
            writer = HmmWriter(self.datadir + '/human-beings/' + self.args.human + '/' + self.args.naivety,
                               hmmdir, gene, self.args.naivety, self.germline_seqs[utils.get_region(gene)][gene], v_right_length=v_right_length)
            writer.write()

    # ----------------------------------------------------------------------------------------
    def write_all_hmms(self, v_right_length):
        if not os.path.exists(self.default_hmm_dir):
            os.makedirs(self.default_hmm_dir)
        for region in utils.regions:
            self.write_specified_hmms(self.default_hmm_dir, self.germline_seqs[region], v_right_length)

    # ----------------------------------------------------------------------------------------
    def write_hmm_input(self, csv_fname, single_gene_precluster=False, preclusters=None):  # TODO use different input files for the two hmm steps
        pair_scores = {}
        with opener('w')(csv_fname) as csvfile:
            # write header
            header = ['name', 'second_name', 'k_v_guess', 'k_d_guess', 'v_fuzz', 'd_fuzz', 'only_genes', 'seq', 'second_seq']  # I wish I had a good c++ csv reader 
            csvfile.write(' '.join(header) + '\n')

            # write a line for each query sequence (or pair of them)
            for query_name in self.reco_info:
                info = self.sw_info[query_name]

                # make sure we don't need to rewrite the hmm model files
                if abs(info['v_right_length'] - self.default_v_right_length) >= self.safety_buffer_that_I_should_not_need:
                    print 'WARNING VRIGHT ', info['v_right_length'], self.default_v_right_length, (info['v_right_length']-self.default_v_right_length)
                # assert abs(info['v_right_length'] - self.default_v_right_length) < self.safety_buffer_that_I_should_not_need
                # I think we can remove these versions (we never see them), but I'm putting a check in here just in case
                assert len(re.findall('J[123]P', info['best']['j'])) == 0

                only_genes = info['all']
                if single_gene_precluster:
                    only_genes = ':'.join(self.match_names_for_all_queries)

                if self.args.pair and not single_gene_precluster:
                    for second_query_name in self.reco_info:  # TODO I can probably get away with skipping a lot of these pairs -- if A clusters with B and B with C, don't run A against C
                        # second_info = self.sw_info[second_query_name]  # TODO dammit I'm still only using info from the first query
                        if second_query_name == query_name:
                            continue
                        if self.get_score_index(query_name, second_query_name) in pair_scores:  # already wrote this pair to the file
                            continue

                        if preclusters != None:  # if we've already run preclustering, skip the pairs that we know aren't matches
                            if preclusters.query_clusters[query_name] != preclusters.query_clusters[second_query_name]:
                                continue

                        pair_scores[self.get_score_index(query_name, second_query_name)] = 0  # set the value to zero so we know we alrady added this pair to the csv file
                        csvfile.write('%s %s %d %d %d %d %s %s %s\n' %
                                             (query_name, second_query_name, info['k_v'], info['k_d'], info['v_fuzz'], info['d_fuzz'], only_genes,
                                              self.reco_info[query_name]['seq'], self.reco_info[second_query_name]['seq']))
                else:
                    # assert self.args.algorithm == 'viterbi'  # TODO allow non-pair forward
                    csvfile.write('%s x %d %d %d %d %s %s x\n' % (query_name, info['k_v'], info['k_d'], info['v_fuzz'], info['d_fuzz'], only_genes, self.reco_info[query_name]['seq']))

    # ----------------------------------------------------------------------------------------
    def run_stochhmm(self, csv_infname, csv_outfname, single_gene_precluster=False):
        start = time.time()

        # build the command line
        cmd_str = './stochhmm'
        cmd_str += ' --algorithm ' + self.args.algorithm
        if self.args.pair and not single_gene_precluster:
            cmd_str += ' --pair 1'
        cmd_str += ' --n_best_events ' + str(self.args.n_best_events)
        cmd_str += ' --debug ' + str(self.args.debug)
        cmd_str += ' --hmmdir ' + self.default_hmm_dir
        cmd_str += ' --infile ' + csv_infname
        cmd_str += ' --outfile ' + csv_outfname
        if single_gene_precluster:  # write single gene probabilities for preclustering
            cmd_str += ' --single_gene_probs 1'

        if self.args.debug == 2:  # not sure why, but popen thing hangs with debug 2
            check_call(cmd_str, shell=True)  # um, not sure which I want here, but it doesn't really matter. TODO kinda
            sys.exit()

        # run hmm
        hmm_proc = Popen(cmd_str, shell=True, stdout=PIPE, stderr=PIPE)
        hmm_proc.wait()
        hmm_out, hmm_err = hmm_proc.communicate()
        print 'OUT\n',hmm_out
        if hmm_proc.returncode != 0:
            print 'aaarrrrrrrgggggh\n',hmm_err,hmm_out
            sys.exit()

        print 'hmm run time: %.3f' % (time.time() - start)

    # ----------------------------------------------------------------------------------------
    def read_stochhmm_output(self, csv_outfname, precluster=False):
        """ Read hmm output file """
        if precluster:
            with opener('r')(csv_outfname) as csv_outfile:
                reader = csv.DictReader(csv_outfile)
                for line in reader:
                    self.precluster_info[line['unique_id']] = {}
                    score_list = line['scores'].split(';')
                    for score_str in score_list:
                        gene = score_str.split(':')[0]
                        score = float(score_str.split(':')[1])
                        self.precluster_info[line['unique_id']][gene] = score
        else:
            with opener('r')(csv_outfname) as csv_outfile:
                reader = csv.DictReader(csv_outfile)
                last_id = None
                for line in reader:
                    if last_id != self.get_score_index(line['unique_id'], line['second_unique_id']):
                        from_same_event = False
                        if self.args.pair:
                            from_same_event = self.reco_info[line['unique_id']]['reco_id'] == self.reco_info[line['second_unique_id']]['reco_id']
                        print '%20s %20s   %d' % (line['unique_id'], line['second_unique_id'], from_same_event),
                    if self.args.algorithm == 'viterbi':
                        print ''
                        # if this is the first line for this query (or query pair), print the true event
                        if last_id != self.get_score_index(line['unique_id'], line['second_unique_id']):
                            print '    true:'
                            utils.print_reco_event(self.germline_seqs, self.reco_info[line['unique_id']], 0, 0, extra_str='    ')
                            if self.args.pair:
                                print '      and maybe'
                                utils.print_reco_event(self.germline_seqs, self.reco_info[line['second_unique_id']], 0, 0, True, '    ')
                            print '    inferred:'
                        utils.print_reco_event(self.germline_seqs, line, 0, 0, extra_str='    ')
                        if self.args.pair:
                            tmpseq = line['seq']  # TODO oh, man, that's a cludge
                            line['seq'] = line['second_seq']
                            utils.print_reco_event(self.germline_seqs, line, 0, 0, True, extra_str='    ')
                            line['seq'] = tmpseq
                    else:
                        print '  ',line['score']
                        with opener('a')(self.pairscorefname) as pairscorefile:
                            pairscorefile.write('%s,%s,%f\n' % (line['unique_id'], line['second_unique_id'], float(line['score'])))
                    last_id = self.get_score_index(line['unique_id'], line['second_unique_id'])
    
    # ----------------------------------------------------------------------------------------
    def calculate_scorespace_distances(self, outfname):
        start = time.time()
        precluster_pairscores = {}
        mean_per_gene = {}
        dbgfname = 'output/' + self.args.human + '-scorespace-preclustering.txt'
        if os.path.exists(dbgfname):
            os.remove(dbgfname)
        with opener('w')(outfname) as outfile:
            writer = csv.DictWriter(outfile, ('unique_id', 'second_unique_id', 'score'))
            writer.writeheader()
            for query,info in self.precluster_info.iteritems():
                total = 0.0
                entries = 0
                for gene,score in info.iteritems():
                    if score == float('-inf'):
                        continue
                    total += score
                    entries += 1
                mean_per_gene[query] = total / entries
            for query_1,info_1 in self.precluster_info.iteritems():
                for query_2,info_2 in self.precluster_info.iteritems():
                    if query_1 == query_2:
                        continue
                    if self.get_score_index(query_1, query_2) in precluster_pairscores:
                        continue
    
                    gs_1 = set(info_1.keys())
                    gs_2 = set(info_2.keys())
                    assert gs_1.issubset(gs_2) and gs_2.issubset(gs_1)
                    total = 0.0
                    for gene in gs_1:
                        score_1 = info_1[gene]
                        score_2 = info_2[gene]
                        if score_1 == float('-inf') and score_2 != float('-inf'):  # if one, but not the other, is -inf, assume they're not in the same cluster
                            continue
                        if score_2 == float('-inf') and score_1 != float('-inf'):
                            continue
                        if score_1 != float('-inf') and score_2 != float('-inf'):
                            total += (score_1 / mean_per_gene[query_1] - score_2 / mean_per_gene[query_2])**2  # NOTE rescaling doesn't really seem to help... oh, well
                    from_same_event = self.reco_info[query_1]['reco_id'] == self.reco_info[query_2]['reco_id']
                    with opener('a')(dbgfname) as dbgfile:
                        dbgfile.write('%20s %20s %2d %-.4f\n' % (query_1,query_2,from_same_event,total))
                    precluster_pairscores[self.get_score_index(query_1, query_2)] = total
                    writer.writerow({'unique_id':query_1, 'second_unique_id':query_2, 'score':total})
        print 'scorespace cluster time: %.3f' % (time.time() - start)
