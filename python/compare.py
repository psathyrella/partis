#!/usr/bin/env python

import sys
sys.path.append('python')

import plotting

plotdir = '/var/www/sharing/dralph/partis/few-genes'
plotting.compare_directories(plotdir + '/hmm-vs-sw',
                             plotdir + '/hmm/plots', 'hmm',
                             plotdir + '/sw/plots', 'sw',
                             xtitle='inferred - true')
