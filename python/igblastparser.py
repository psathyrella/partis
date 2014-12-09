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

# # ----------------------------------------------------------------------------------------
# def merge_chunks(all_info, line, gene):
#     assert gene in all_info
    
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
            while line != '':
                query_name = int(line.split()[1])
                if self.args.debug:
                    print query_name
                if query_name not in self.seqinfo:
                    print 'ERROR %d not in reco info' % query_name
                    sys.exit()
                info[query_name] = {}
                query_lines = []
                line = infile.readline()
                while line.find('<b>Query=') != 0:
                    query_lines.append(line.strip())
                    line = infile.readline()
                # then add the query to <info[query_name]>
                self.process_query(info[query_name], query_name, query_lines)

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

        for block in blocks:
            self.process_single_block(block, query_name, qr_info)

        print '    complete matches:'
        for region in utils.regions:
            print '    %s %s' % (utils.color_gene(qr_info[region + '_gene']), qr_info[region + '_gl_seq'])
        sys.exit()
            
        # query_lines = []
        # while 
        # query_lines.append
        
            
        # soup = BeautifulSoup(infile)
        # paragraphs = soup.find_all('b')

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
                

    def process_single_block(self, block, query_name, qr_info):
        assert block[0].find('Query_') == 0
        vals = block[0].split()
        qr_start = int(vals[1])
        qr_seq = vals[2]
        qr_end = int(vals[3])
        assert qr_seq in self.seqinfo[query_name]['seq']
        if self.args.debug:
            print '  query %3d %3d %s' % (qr_start, qr_end, qr_seq)
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
                        print '    %s match: %s' % (region, clean_alignment_crap(qr_seq, gl_seq))
                    qr_info[region + '_gl_seq'] = qr_info[region + '_gl_seq'] + clean_alignment_crap(qr_seq, gl_seq)
                else:
                    continue
            else:
                qr_info[region + '_gene'] = gene
                qr_info[region + '_gl_seq'] = clean_alignment_crap(qr_seq, gl_seq)
                if self.args.debug:
                    print '    %s match: %s' % (region, clean_alignment_crap(qr_seq, gl_seq))

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
