#!/usr/bin/env python
import argparse
import csv
import sys
from subprocess import check_call
import re
from bs4 import BeautifulSoup
import os

sys.path.insert(1, './python')
from opener import opener
import utils
from joinparser import resolve_overlapping_matches  # why the hell is this import so slow?
# from ihhhmmmparser import FileKeeper
from performanceplotter import PerformancePlotter
from imgtparser import equivalent_genes, just_always_friggin_skip, get_genes_to_skip, genes_to_skip, genes_to_use

# ----------------------------------------------------------------------------------------
def find_qr_bounds(global_qr_start, global_qr_end, gl_match_seq):
    """ Return the start and end of this match in the query coordinate system """
    # find first matching character
    # print 'find_qr_bounds', global_qr_start, global_qr_end, gl_match_seq
    istart, iend = 0, len(gl_match_seq)
    for ic in range(len(gl_match_seq)):
        ch = gl_match_seq[ic]
        if ch == '.' or ch in utils.nukes:
            istart = ic
            break
    # and first non-matching character after end of match
    for ic in range(istart, len(gl_match_seq)):
        ch = gl_match_seq[ic]
        if ch != '.' and ch not in utils.nukes:
            iend = ic
            break
    assert istart >= 0 and iend >= 0
    assert istart <= len(gl_match_seq) and iend <= len(gl_match_seq)
    # print '    %d + %d = %d, %d - (%d - %d) = %d' % (global_qr_start, istart, global_qr_start+istart, global_qr_end, len(gl_match_seq), iend, global_qr_end - (len(gl_match_seq)-iend))
    return (global_qr_start + istart, global_qr_end - (len(gl_match_seq) - iend))

# ----------------------------------------------------------------------------------------
def clean_alignment_crap(query_seq, match_seq):
    if len(query_seq) != len(match_seq):
        print 'ERROR in clean_alignment_crap(): not the same length'
        print '    %s' % query_seq
        print '    %s' % match_seq
        sys.exit()

    final_match_seq = []
    for inuke in range(len(query_seq)):
        mnuke = match_seq[inuke]
        qnuke = query_seq[inuke]
        if mnuke == '-':
            continue
        elif mnuke == '.':
            final_match_seq.append(qnuke)
        else:
            # if mnuke == 'N':
            #     mnuke = 'A'
            #     print 'WARNING replacing N with A in germlines'
            if mnuke not in utils.nukes:
                print 'ERROR unexpected character %s' % mnuke
                sys.exit()
            final_match_seq.append(mnuke)
            
    return ''.join(final_match_seq)

# ----------------------------------------------------------------------------------------
genes_actually_skipped = {}
class IgblastParser(object):
    def __init__(self, args):
        self.args = args

        self.germline_seqs = utils.read_glfo(self.args.datadir, remove_N_nukes=True)['seqs']

        self.perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir, 'igblast')
        self.n_total, self.n_partially_failed, self.n_skipped = 0, 0, 0

        # get sequence info that was passed to igblast
        self.seqinfo = {}
        with opener('r')(self.args.simfname) as simfile:
            reader = csv.DictReader(simfile)
            iline = 0
            for line in reader:
                if self.args.n_queries > 0 and iline >= self.args.n_queries:
                    break
                iline += 1
                if self.args.queries != None and int(line['unique_id']) not in self.args.queries:
                    continue
                if len(re.findall('_[FP]', line['j_gene'])) > 0:
                    line['j_gene'] = line['j_gene'].replace(re.findall('_[FP]', line['j_gene'])[0], '')
                self.seqinfo[int(line['unique_id'])] = line

        print 'reading', self.args.infname

        get_genes_to_skip(self.args.infname, self.germline_seqs, method='igblast', debug=False)

        paragraphs = None
        info = {}
        with opener('r')(self.args.infname) as infile:
            line = infile.readline()
            # first find the start of the next query's section
            while line.find('<b>Query=') != 0:
                line = infile.readline()
            # then keep going till eof
            iquery = 0
            while line != '':
                if self.args.n_queries > 0 and iquery >= self.args.n_queries:
                    break
                # first find the query name
                query_name = int(line.split()[1])
                # and collect the lines for this query
                query_lines = []
                line = infile.readline()
                while line.find('<b>Query=') != 0:
                    query_lines.append(line.strip())
                    line = infile.readline()
                    if line == '':
                        break
                iquery += 1
                # then see if we want this query
                if self.args.queries != None and query_name not in self.args.queries:
                    continue
                if query_name not in self.seqinfo:
                    print 'ERROR %d not in reco info' % query_name
                    sys.exit()
                if self.args.debug:
                    print query_name
                # and finally add the query to <info[query_name]>
                info[query_name] = {'unique_id':query_name}
                self.n_total += 1
                self.process_query(info[query_name], query_name, query_lines)

        self.perfplotter.plot()
        print 'partially failed: %d / %d = %f' % (self.n_partially_failed, self.n_total, float(self.n_partially_failed) / self.n_total)
        print 'skipped: %d / %d = %f' % (self.n_skipped, self.n_total, float(self.n_skipped) / self.n_total)
        for g, n in genes_actually_skipped.items():
            print '  %d %s' % (n, utils.color_gene(g))

    # ----------------------------------------------------------------------------------------
    def process_query(self, qr_info, query_name, query_lines):
        # split query_lines up into blocks
        blocks = []
        for line in query_lines:
            if line.find('Query_') == 0:
                blocks.append([])
            if len(line) == 0:
                continue
            if len(re.findall('<a name=#_[0-9][0-9]*_IGH', line)) == 0 and line.find('Query_') != 0:
                continue
            if len(blocks) == 0:
                print 'wtf? %s' % query_name  # it's probably kicking a reverse match
                self.perfplotter.add_partial_fail(self.seqinfo[query_name], qr_info)  # NOTE that's really a total failure
                self.n_partially_failed += 1
                return
            blocks[-1].append(line)

        # then process each block
        for block in blocks:
            self.process_single_block(block, query_name, qr_info)
            if 'skip_gene' in qr_info:
                self.n_skipped += 1
                return
            if 'fail' in qr_info:
                self.perfplotter.add_partial_fail(self.seqinfo[query_name], qr_info)
                self.n_partially_failed += 1
                return

        for region in utils.regions:
            if region + '_gene' not in qr_info:
                print '    %d: no %s match' % (query_name, region)
                self.perfplotter.add_partial_fail(self.seqinfo[query_name], qr_info)
                self.n_partially_failed += 1
                return

        # expand v match to left end and j match to right end
        qr_info['v_5p_del'] = 0
        qr_info['fv_insertion'] = ''
        if qr_info['match_start'] > 0:
            if self.args.debug:
                print '    add to v left:', self.seqinfo[query_name]['seq'][ : qr_info['match_start']]
            qr_info['seq'] = self.seqinfo[query_name]['seq'][ : qr_info['match_start']] + qr_info['seq']

        qr_info['j_3p_del'] = 0
        qr_info['jf_insertion'] = ''
        if len(self.seqinfo[query_name]['seq']) > qr_info['match_end']:
            if self.args.debug:
                print '    add to j right:', self.seqinfo[query_name]['seq'][ qr_info['match_end'] - len(self.seqinfo[query_name]['seq']) : ]
            qr_info['seq'] = qr_info['seq'] + self.seqinfo[query_name]['seq'][ qr_info['match_end'] - len(self.seqinfo[query_name]['seq']) : ]

        for boundary in utils.boundaries:
            start = qr_info[boundary[0] + '_qr_bounds'][1]
            end = qr_info[boundary[1] + '_qr_bounds'][0]
            qr_info[boundary + '_insertion'] = qr_info['seq'][start : end]

        for region in utils.regions:
            start = qr_info[region + '_qr_bounds'][0]
            end = qr_info[region + '_qr_bounds'][1]
            qr_info[region + '_qr_seq'] = qr_info['seq'][start : end]

        try:
            resolve_overlapping_matches(qr_info, self.args.debug, self.germline_seqs)
        except AssertionError:
            print '    %s: apportionment failed' % query_name
            self.perfplotter.add_partial_fail(self.seqinfo[query_name], qr_info)
            self.n_partially_failed += 1
            return

        if self.args.debug:
            print '  query seq:', qr_info['seq']
            for region in utils.regions:
                true_gene = self.seqinfo[query_name][region + '_gene']
                infer_gene = qr_info[region + '_gene']
                if utils.are_alleles(infer_gene, true_gene):
                    regionstr = utils.color('bold', utils.color('blue', region))
                    truestr = ''  #'(originally %s)' % match_name
                else:
                    regionstr = utils.color('bold', utils.color('red', region))
                    truestr = '(true: %s)' % utils.color_gene(true_gene).replace(region, '')
                # print '  %s %s %s' % (regionstr, utils.color_gene(infer_gene).replace(region, ''), truestr)

                print '    %s %3d %3d %s %s %s' % (regionstr, qr_info[region + '_qr_bounds'][0], qr_info[region + '_qr_bounds'][1], utils.color_gene(infer_gene).replace(region, ''), truestr, qr_info[region + '_gl_seq'])
        for boundary in utils.boundaries:
            start = qr_info[boundary[0] + '_qr_bounds'][1]
            end = qr_info[boundary[1] + '_qr_bounds'][0]
            qr_info[boundary + '_insertion'] = qr_info['seq'][start : end]
            if self.args.debug:
                print '   ', boundary, qr_info[boundary + '_insertion']

        self.perfplotter.evaluate(self.seqinfo[query_name], qr_info)
        # for key, val in qr_info.items():
        #     print key, val
        if self.args.debug:
            utils.print_reco_event(self.germline_seqs, self.seqinfo[query_name], label='true:', extra_str='  ')
            utils.print_reco_event(self.germline_seqs, qr_info, extra_str=' ')
            
    # ----------------------------------------------------------------------------------------
    def process_single_block(self, block, query_name, qr_info):
        assert block[0].find('Query_') == 0
        vals = block[0].split()
        qr_start = int(vals[1]) - 1  # converting from one-indexed to zero-indexed
        qr_seq = vals[2]
        qr_end = int(vals[3])  # ...and from inclusive of both bounds to normal programming conventions
        if qr_seq not in self.seqinfo[query_name]['seq']:
            if '-' in qr_seq:
                print '    %s: insertion inside query seq, treating as partial failure' % query_name
                qr_info['fail'] = True
                return
            else:
                print '  ERROR query seq from igblast info not found in original query seq for %d' % query_name
                print '    %s' % qr_seq
                print '    %s' % self.seqinfo[query_name]['seq']
                sys.exit()

        if 'seq' in qr_info:
            qr_info['seq'] += qr_seq
        else:
            qr_info['seq'] = qr_seq


        # keep track of the absolute first and absolute last bases matched so we can later work out the fv and jf insertions
        if 'match_start' not in qr_info or qr_start < qr_info['match_start']:
            qr_info['match_start'] = qr_start
        if 'match_end' not in qr_info or qr_end > qr_info['match_end']:
            qr_info['match_end'] = qr_end

        # ----------------------------------------------------------------------------------------
        # skipping bullshit
        def skip_gene(gene):
            if self.args.debug:
                print '    %s in list of genes to skip' % utils.color_gene(gene)
            if gene not in genes_actually_skipped:
                genes_actually_skipped[gene] = 0
            genes_actually_skipped[gene] += 1
            qr_info['skip_gene'] = True

        if self.args.debug:
            print '      query: %3d %3d %s' % (qr_start, qr_end, qr_seq)
        for line in block[1:]:
            gene = line[line.rfind('IGH') : line.rfind('</a>')]
            region = utils.get_region(gene)
            true_gene = self.seqinfo[query_name][region + '_gene']

            for gset in equivalent_genes:
                if gene in gset and true_gene in gset and gene != true_gene:  # if the true gene and the inferred gene are in the same equivalence set, treat it as correct, i.e. just pretend it inferred the right name
                    if self.args.debug:
                        print '   %s: replacing name %s with true name %s' % (query_name, gene, true_gene)
                    gene = true_gene

            if gene in just_always_friggin_skip:
                continue  # go on to the next match

            if not self.args.dont_skip_or15_genes and '/OR1' in true_gene:
                skip_gene(true_gene)
                return

            if self.args.skip_missing_genes:
                if gene in genes_to_skip:
                    continue  # go on to the next match
                    # skip_gene(gene)
                    # return
                if true_gene in genes_to_skip:
                    skip_gene(true_gene)
                    return

            if gene not in self.germline_seqs[region]:
                print '    %s: %s not in germlines (skipping)' % (query_name, gene)
                skip_gene(gene)
                return
                
            vals = line.split()
            gl_start = int(vals[-3]) - 1  # converting from one-indexed to zero-indexed
            gl_seq = vals[-2]
            gl_end = int(vals[-1])  # ...and from inclusive of both bounds to normal programming conventions

            if region + '_gene' in qr_info:
                if qr_info[region + '_gene'] == gene:
                    if self.args.debug:
                        print '        %s match: %s' % (region, clean_alignment_crap(qr_seq, gl_seq))
                    qr_info[region + '_gl_seq'] = qr_info[region + '_gl_seq'] + clean_alignment_crap(qr_seq, gl_seq)
                    if gl_end > len(self.germline_seqs[region][gene]):  # not really sure what's wrong... but it seems to be rare
                        qr_info['fail'] = True
                        return
                    qr_info[region + '_3p_del'] = len(self.germline_seqs[region][gene]) - gl_end
                    qr_info[region + '_qr_bounds'] = (qr_info[region + '_qr_bounds'][0], find_qr_bounds(qr_start, qr_end, gl_seq)[1])
                else:
                    continue
            else:
                qr_info[region + '_gene'] = gene
                qr_info[region + '_gl_seq'] = clean_alignment_crap(qr_seq, gl_seq)
                # deletions
                qr_info[region + '_5p_del'] = gl_start
                assert gl_end <= len(self.germline_seqs[region][gene])
                qr_info[region + '_3p_del'] = len(self.germline_seqs[region][gene]) - gl_end
                # bounds
                qr_info[region + '_qr_bounds'] = find_qr_bounds(qr_start, qr_end, gl_seq)
                if self.args.debug:
                    print '        %s match: %s' % (region, clean_alignment_crap(qr_seq, gl_seq))

# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n-queries', type=int, default=-1)
    parser.add_argument('--queries')
    parser.add_argument('--plotdir', required=True)
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
    parser.add_argument('--datadir', default='data/imgt')
    parser.add_argument('--infname', required=True)  #, default='data/performance/igblast/igblast.html')
    parser.add_argument('--simfname', required=True)
    parser.add_argument('--skip-missing-genes', action='store_true')
    parser.add_argument('-dont-skip-or15-genes', action='store_true', help='by default skip all the genes with the /OR1[56] bullshit, since they don\'t seem to be in imgt\'s output')
    args = parser.parse_args()
    args.queries = utils.get_arg_list(args.queries, intify=True)
    
    # if os.path.isdir('data/performance/igblast'):
    #     print 'skipping tar xzf \'cause output\'s already there'
    # else:
    #     print 'untgzing...'
    #     check_call(['tar', 'xzf', 'data/performance/igblast.tgz', '-C', 'data/performance/'])  # untar the igblast output
    igblastparser = IgblastParser(args)
