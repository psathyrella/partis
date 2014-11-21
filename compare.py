#!/usr/bin/env python

import sys
sys.path.append('python')

import plotting

# label = 'check-new-imgt'
# plotdir = '/var/www/sharing/dralph/partis/performance/'
# plotting.compare_directories(plotdir + '/' + label + '/ihhhmmm-vs-partis-vs-imgt',
#                              dirs = [plotdir + '/' + label + '/hmm/plots',
#                                      plotdir + '/' + label + '/sw/plots',
#                                      plotdir + 'ihhhmmm/' + label + '/plots',
#                                      plotdir + 'imgt/' + label + '/plots'],
#                              names = ['partis', 'sw', 'iHMMuneAlign', 'imgt'])
label = 'check-new-imgt'
plotdir = '/var/www/sharing/dralph/partis/performance/'
plotting.compare_directories(plotdir + '/' + label + '/train-vs-test',
                             dirs = [plotdir + '/' + label + '/hmm/plots',
                                     plotdir + '/' + label + '-train-on-head-test-on-tail/hmm/plots',
                                     plotdir + '/' + label + '-train-on-tail-test-on-head/hmm/plots',
                                     plotdir + '/' + label + '-train-on-tail-test-on-head/sw/plots'],
                             names = ['full', 'head-tail', 'tail-head', 'sw'])
