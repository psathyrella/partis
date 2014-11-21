#!/usr/bin/env python

import sys
sys.path.append('python')

import plotting

label = 'check-new-imgt'
plotdir = '/var/www/sharing/dralph/partis/performance/'
plotting.compare_directories(plotdir + '/' + label + '/imgt-vs-partis',
                             plotdir + '/' + label + '/sw/plots', 'sw',
                             plotdir + '/' + label + '/hmm/plots', 'partis',
                             # dir3=plotdir + 'ihhhmmm/' + label + '/plots', name3='iHMMuneAlign')
                             dir3=plotdir + 'imgt/' + label + '/plots', name3='imgt')
