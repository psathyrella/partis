#!/usr/bin/env python

import sys
sys.path.append('python')

import plotting

plotdir = '/var/www/sharing/dralph/partis/performance/'
plotting.compare_directories(plotdir + '/new-imgt/ihhhmmm-vs-partis',
                             plotdir + 'new-imgt-skip-simulation/hmm/plots', 'partis',
                             plotdir + 'ihhhmmm/new-imgt/plots', 'iHMMuneAlign',
                             xtitle='inferred - true')
