#!/usr/bin/env python
import argparse
import csv
import sys
import re
from bs4 import BeautifulSoup
import os

from opener import opener
import utils
# import joinparser  # why the hell is this import so slow?
# from ihhhmmmparser import FileKeeper
from performanceplotter import PerformancePlotter

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
            assert mnuke in utils.nukes
            final_match_seq.append(mnuke)
            
    return ''.join(final_match_seq)

# ----------------------------------------------------------------------------------------
class IgblastParser(object):
    def __init__(self, args):
        self.args = args

        self.germline_seqs = utils.read_germlines(self.args.datadir)

        perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir, 'igblast')

        # get sequence info that was passed to igblast
        self.seqinfo = {}
        with opener('r')(self.args.simfname) as simfile:
            reader = csv.DictReader(simfile)
            iline = 0
            for line in reader:
                if self.args.queries != None and int(line['unique_id']) not in self.args.queries:
                    continue
                if len(re.findall('_[FP]', line['j_gene'])) > 0:  # TODO remove this
                    line['j_gene'] = line['j_gene'].replace(re.findall('_[FP]', line['j_gene'])[0], '')
                self.seqinfo[int(line['unique_id'])] = line
                iline += 1
                if self.args.n_max_queries > 0 and iline >= self.args.n_max_queries:
                    break

        paragraphs = None
        print 'reading', self.args.infname
        info = {}
        with opener('r')(self.args.infname) as infile:
            line = infile.readline()
            # first find the start of the next query's section
            while line.find('<b>Query=') != 0:
                line = infile.readline()
            # then keep going till eof
            iquery = 0
            while line != '':
                query_name = int(line.split()[1])
                if query_name not in self.seqinfo:
                    print 'ERROR %d not in reco info' % query_name
                    sys.exit()
                if self.args.debug:
                    print query_name
                    utils.print_reco_event(self.germline_seqs, self.seqinfo[query_name], label='true:')
                info[query_name] = {}
                query_lines = []
                line = infile.readline()
                while line.find('<b>Query=') != 0:
                    query_lines.append(line.strip())
                    line = infile.readline()
                # then add the query to <info[query_name]>
                self.process_query(info[query_name], query_name, query_lines)
                iquery += 1
                if self.args.n_max_queries > 0 and iquery >= self.args.n_max_queries:
                    break

        # n_failed, n_total, n_not_found, n_found = 0, 0, 0, 0
        # for unique_id in self.seqinfo:
        #     if self.args.debug:
        #         print unique_id,
        #     igblastinfo = []
        #     for tag in paragraphs:  # NOTE this loops over everything an awful lot of times. Shouldn't really matter for now, though
        #         assert tag.text == 'Query='
        #         floating_info = tag.next_sibling.split()  # I wish they'd put this inside a tag of some sort
        #         query_name = int(floating_info[0])
        #         length = int(floating_info[1].split('=')[1])
        #         assert length == len(self.seqinfo[query_name]['seq'])
                

    # ----------------------------------------------------------------------------------------
    def process_query(self, qr_info, query_name, query_lines):
        # split query_lines up into blocks
        blocks = []
        for line in query_lines:
            if line.find('Query_') == 0:
                blocks.append([])
            if len(line) == 0:
                continue
            if '<a name=#_0_IGH' not in line and line.find('Query_') != 0:
                continue
            blocks[-1].append(line)

        # then process each block
        for block in blocks:
            self.process_single_block(block, query_name, qr_info)

        if self.args.debug:
            print '  query seq:', qr_info['seq']
        for region in utils.regions:
            print '    %s %3d %3d %s %s' % (region, qr_info[region + '_qr_bounds'][0], qr_info[region + '_qr_bounds'][1], utils.color_gene(qr_info[region + '_gene']), qr_info[region + '_gl_seq'])
        for boundary in utils.boundaries:
            start = qr_info[boundary[0] + '_qr_bounds'][1]
            end = qr_info[boundary[1] + '_qr_bounds'][0]
            qr_info[boundary + '_insertion'] = qr_info['seq'][start : end]
            print '   ', boundary, qr_info[boundary + '_insertion']
        qr_info['fv_insertion'] = qr_info['seq'][ : qr_info['v_5p_del']]
        if qr_info['j_3p_del'] > 0:  # I *still* wish slices behaved this way for -0
            qr_info['jf_insertion'] = qr_info['seq'][ : -qr_info['j_3p_del']]
        else:
            qr_info['jf_insertion'] = ''
        utils.print_reco_event(self.germline_seqs, qr_info)
            
    # ----------------------------------------------------------------------------------------
    def process_single_block(self, block, query_name, qr_info):
        assert block[0].find('Query_') == 0
        vals = block[0].split()
        qr_start = int(vals[1]) - 1  # converting from one-indexed to zero-indexed
        qr_seq = vals[2]
        qr_end = int(vals[3])  # ...and from inclusive of both bounds to normal programming conventions
        assert qr_seq in self.seqinfo[query_name]['seq']
        if 'seq' in qr_info:
            qr_info['seq'] += qr_seq
        else:
            qr_info['seq'] = qr_seq
        if self.args.debug:
            print '      query: %3d %3d %s' % (qr_start, qr_end, qr_seq)
        for line in block[1:]:
            gene = line[line.rfind('IGH') : line.rfind('</a>')]
            region = utils.get_region(gene)
            vals = line.split()
            gl_start = int(vals[-3]) - 1  # converting from one-indexed to zero-indexed
            gl_seq = vals[-2]
            gl_end = int(vals[-1])  # ...and from inclusive of both bounds to normal programming conventions

            if region + '_gene' in qr_info:
                if qr_info[region + '_gene'] == gene:
                    if self.args.debug:
                        print '        %s match: %s' % (region, clean_alignment_crap(qr_seq, gl_seq))
                    qr_info[region + '_gl_seq'] = qr_info[region + '_gl_seq'] + clean_alignment_crap(qr_seq, gl_seq)
                    assert gl_end <= len(self.germline_seqs[region][gene])
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
    parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
    parser.add_argument('--label', default='check-new-imgt')
    parser.add_argument('--n-max-queries', type=int, default=-1)
    parser.add_argument('--queries')
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
    parser.add_argument('--datadir', default='data/imgt')
    parser.add_argument('--infname', default='data/performance/igblast/igblast.html')
    args = parser.parse_args()
    args.queries = utils.get_arg_list(args.queries)
    
    args.simfname = 'caches/recombinator/performance/' + args.label + '/simu.csv'
    args.plotdir = os.getenv('www') + '/partis/performance/igblast/' + args.label
    igblastparser = IgblastParser(args)
