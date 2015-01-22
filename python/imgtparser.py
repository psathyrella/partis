#!/usr/bin/env python
import csv
import sys
import glob
import argparse
import os
from subprocess import check_call
import re
from bs4 import BeautifulSoup

sys.path.insert(1, './python')
from opener import opener
import utils
import joinparser

from performanceplotter import PerformancePlotter
# - convert recombinator csv to fasta: /home/dralph/Dropbox/bin/csv2fasta

# ----------------------------------------------------------------------------------------
class IMGTParser(object):
    def __init__(self, args):
        self.args = args

        self.germline_seqs = utils.read_germlines(self.args.datadir)

        perfplotter = PerformancePlotter(self.germline_seqs, self.args.plotdir, 'imgt')

        # get sequence info that was passed to imgt
        self.seqinfo = {}
        with opener('r')(self.args.simfname) as simfile:
            reader = csv.DictReader(simfile)
            iline = 0
            for line in reader:
                if self.args.queries != None and line['unique_id'] not in self.args.queries:
                    continue
                if len(re.findall('_[FP]', line['j_gene'])) > 0:
                    line['j_gene'] = line['j_gene'].replace(re.findall('_[FP]', line['j_gene'])[0], '')
                self.seqinfo[line['unique_id']] = line
                iline += 1
                if self.args.n_max_queries > 0 and iline >= self.args.n_max_queries:
                    break

        paragraphs, csv_info = None, None
        if self.args.infname != None and '.html' in self.args.infname:
            print 'reading', self.args.infname
            with opener('r')(self.args.infname) as infile:
                soup = BeautifulSoup(infile)
                paragraphs = soup.find_all('pre')

        n_failed, n_total, n_not_found, n_found = 0, 0, 0, 0
        for unique_id in self.seqinfo:
            if self.args.debug:
                print unique_id,
            imgtinfo = []
            # print 'true'
            # utils.print_reco_event(self.germline_seqs, self.seqinfo[unique_id])
            if self.args.infname != None and '.html' in self.args.infname:
                for pre in paragraphs:  # NOTE this loops over everything an awful lot of times. Shouldn't really matter for now, though
                    if unique_id in pre.text:
                        imgtinfo.append(pre.text)
            else:
                n_total += 1
                assert self.args.infname == None
                infnames = glob.glob(self.args.indir + '/' + unique_id + '*')
                assert len(infnames) <= 1
                if len(infnames) != 1:
                    if self.args.debug:
                        print ' couldn\'t find it'
                    n_not_found += 1
                    continue
                n_found += 1
                with opener('r')(infnames[0]) as infile:
                    full_text = infile.read()
                    if len(re.findall('[123]. Alignment for [VDJ]-GENE', full_text)) < 3:
                        failregions = re.findall('No [VDJ]-GENE has been identified', full_text)
                        if len(failregions) > 0:
                            print '    ', failregions
                        n_failed += 1
                        continue

                    # loop over the paragraphs I want
                    position = full_text.find(unique_id)  # don't need this one
                    for ir in range(4):
                        position = full_text.find(unique_id, position+1)
                        pgraph = full_text[position : full_text.find('\n\n', position+1)]
                        if 'insertion(s) and/or deletion(s) which are not dealt in this release' in pgraph:
                            ir -= 1
                            continue
                        imgtinfo.append(pgraph)  # query seq paragraph

            if len(imgtinfo) == 0:
                print '%s no info' % unique_id
                continue
            else:
                if self.args.debug:
                    print ''
            line = self.parse_query_text(unique_id, imgtinfo)
            try:
                assert 'failed' not in line
                joinparser.add_insertions(line, debug=self.args.debug)
                joinparser.resolve_overlapping_matches(line, debug=self.args.debug, germlines=self.germline_seqs)
            except (AssertionError, KeyError):
                print '    giving up'
                n_failed += 1
                perfplotter.add_partial_fail(self.seqinfo[unique_id], line)
                print '    perfplotter: not sure what to do with a fail'
                continue
            perfplotter.evaluate(self.seqinfo[unique_id], line)
            if self.args.debug:
                utils.print_reco_event(self.germline_seqs, self.seqinfo[unique_id], label='true:')
                utils.print_reco_event(self.germline_seqs, line, label='inferred:')

        perfplotter.plot()
        print 'failed: %d / %d = %f' % (n_failed, n_total, float(n_failed) / n_total)
        if n_not_found > 0:
            print '  not found: %d / %d = %f' % (n_not_found, n_not_found + n_found, n_not_found / float(n_not_found + n_found))

    # ----------------------------------------------------------------------------------------
    def parse_query_text(self, unique_id, query_info):
        if len(query_info) == 0:  # one for the query sequence, then one for v, d, and j
            print 'no info for',unique_id
            return {}
        elif len(query_info) < 4:
            regions_ok = ''
            for info in query_info:
                for region in utils.regions:
                    if 'IGH' + region.upper() in info:
                        regions_ok += region
            for region in utils.regions:
                if region not in regions_ok:
                    print '    ERROR no %s matches' % region
                    return {}
            assert False  # shouldn't get here
        elif len(query_info) != 4:
            print 'info for', unique_id, 'all messed up'
            for info in query_info:
                print info
            sys.exit()

        full_qr_seq = query_info[0].replace('>', '').replace(unique_id, '')  # strip off the unique id
        full_qr_seq = ''.join(full_qr_seq.split()).upper()  # strip off white space and uppercase it
        assert full_qr_seq == self.seqinfo[unique_id]['seq']

        line = {}
        line['unique_id'] = unique_id
        line['seq'] = full_qr_seq
        for ireg in range(len(utils.regions)):
            region = utils.regions[ireg]
            info = query_info[ireg + 1].splitlines()
            while unique_id not in info[0]:  # remove the line marking cdr3 and framework regions
                info.pop(0)
            if len(info) <= 1:
                print info
            assert len(info) > 1
            assert len(info[0].split()) == 2
            qr_seq = info[0].split()[1].upper()  # this line should be '<unique_id> .............<query_seq>'

            imatch = 1  # which match to take
            match_name = str(info[imatch].split()[2])
            # if 'IGHV3-69' in match_name:  # it's not right anyway
            #     line['failed'] = True
            #     return line
            # while unacceptable_match(match_name, self.germline_seqs):
            #         imatch += 1
            #         match_name = str(info[imatch].split()[2])
            #         print '    new match name: %s' % match_name

            gl_seq = info[imatch].split()[4].upper()
            if qr_seq.replace('.', '') not in self.seqinfo[unique_id]['seq']:
                if self.args.debug:
                    print '    qr_seq not foundin seqinfo'
                line['failed'] = True
                return line

            if self.args.debug:
                print '  ', region, match_name
                print '    gl', gl_seq
                print '      ', qr_seq

            # replace the dots (gaps) in the gl match
            new_qr_seq, new_gl_seq = [], []
            for inuke in range(min(len(qr_seq), len(gl_seq))):
                if gl_seq[inuke] == '.':
                    pass
                else:
                    new_qr_seq.append(qr_seq[inuke])  # this should only be out of range if the v match extends through the whole query sequence, i.e. friggin never
                    new_gl_seq.append(gl_seq[inuke])
            for inuke in range(len(gl_seq), len(qr_seq)):
                new_qr_seq.append(qr_seq[inuke])
            for inuke in range(len(qr_seq), len(gl_seq)):
                new_gl_seq.append(gl_seq[inuke])
            qr_seq = ''.join(new_qr_seq)
            gl_seq = ''.join(new_gl_seq)

            # work out the erosions
            qr_ldots = qr_seq.rfind('.') + 1  # first strip off any dots on the left of query seq
            qr_seq = qr_seq[qr_ldots : ]
            gl_seq = gl_seq[qr_ldots : ]
            gl_ldots = gl_seq.rfind('.') + 1  # then remove dots on the left of the germline seq
            qr_seq = qr_seq[gl_ldots : ]
            gl_seq = gl_seq[gl_ldots : ]
            del_5p = qr_ldots + gl_ldots
            jf_insertion = ''
            if region == 'j':
                jf_insertion = qr_seq[len(gl_seq) : ]
            qr_seq = qr_seq[ : len(gl_seq)]  # then strip the right-hand portion of the query sequence that isn't aligned to the germline
            del_3p = len(gl_seq) - len(qr_seq)  # then do the same for the germline overhanging on the right of the query
            gl_seq = gl_seq[ : len(qr_seq)]
            assert len(gl_seq) == len(qr_seq)

            new_gl_seq = []
            for inuke in range(len(gl_seq)):  # replace dashes (matched bases)
                assert gl_seq[inuke] != '.'  # hoping there's no gaps in here
                if gl_seq[inuke] == '-':
                    new_gl_seq.append(qr_seq[inuke])
                else:
                    new_gl_seq.append(gl_seq[inuke])
            gl_seq = ''.join(new_gl_seq)

            if match_name not in self.germline_seqs[region]:
                print '    ERROR couldn\'t find %s in germlines' % match_name
                line['failed'] = True
                return line

            if self.germline_seqs[region][match_name].find(gl_seq) != del_5p:  # why the *@*!! can't they make this consistent?
                if self.germline_seqs[region][match_name].find(gl_seq) < 0:
                    print 'whooooaa'
                    print self.germline_seqs[region][match_name]
                    print gl_seq
                    sys.exit()
                del_5p += self.germline_seqs[region][match_name].find(gl_seq)

            try:
                assert del_5p + len(gl_seq) + del_3p + len(jf_insertion) == len(self.germline_seqs[region][match_name])
            except:
                print '    ERROR lengths failed for %s' % unique_id
                # print del_5p, len(gl_seq), del_3p, del_5p + len(gl_seq) + del_3p , len(self.germline_seqs[region][match_name])
                # print gl_seq
                # print self.germline_seqs[region][match_name]
                line['failed'] = True
                return line
                # assert False

            if self.args.debug:
                utils.color_mutants(gl_seq, qr_seq, ref_label='gl ', extra_str='    ', print_result=True, post_str='    del: %d %d' % (del_5p, del_3p))

            # try:
            #     match_name = joinparser.figure_out_which_damn_gene(self.germline_seqs, match_name, gl_seq, debug=self.args.debug)
            # except:
            #     print 'ERROR couldn\'t figure out the gene for %s' % match_name
            #     return {}

            line[region + '_gene'] = match_name
            line[region + '_qr_seq'] = qr_seq
            line[region + '_gl_seq'] = gl_seq
            line[region + '_5p_del'] = del_5p
            line[region + '_3p_del'] = del_3p
            if region == 'j':
                line['jf_insertion'] = jf_insertion
            
        return line

# ----------------------------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', action='store_true')  # passed on to ROOT when plotting
    parser.add_argument('--n-max-queries', type=int, default=-1)
    parser.add_argument('--queries')
    parser.add_argument('--plotdir', required=True)
    parser.add_argument('--debug', type=int, default=0, choices=[0, 1, 2])
    parser.add_argument('--datadir', default='data/imgt')
    parser.add_argument('--infname')  # input html file, if you chose the 'html' option on the imgt website
    args = parser.parse_args()
    args.queries = utils.get_arg_list(args.queries)
    
    args.simfname = 'data/performance/simu.csv'
    if os.path.isdir('data/performance/imgt'):
        print 'skipping tar xzf \'cause output\'s already there'
    else:
        print 'untgzing...'
        check_call(['tar', 'xzf', 'data/performance/imgt.tgz', '-C', 'data/performance/'])  # untar the imgt output
    args.indir = 'data/performance/imgt/IMGT_HighV-QUEST_individual_files_folder'  # folder with imgt result files
    imgtparser = IMGTParser(args)
