from __future__ import absolute_import, division, unicode_literals
from __future__ import print_function
import sys
import re

from . import utils

# ----------------------------------------------------------------------------------------
def cluster(viterbi_info):
    """ 
    Cluster together sequences with similar rearrangement parameters

    From Vollmers paper:
        Lineage Clustering. IGH sequences were clustered into IGH lineages according
        to similarity in their junctional region. Lineages were created according to the
        following steps. A lineage is formed and populated with one IGH sequence (seed). Then, all
        IGH sequences in the lineages (initially only the seed) are compared with all
        other IGH sequences of the same length using the same V and J segments. If
        their junctional regions (untemplated nucleotides and D segments) are at
        least 90% identical, the IGH sequence is added to the lineage. This process is
        repeated until the lineage does not grow.
    """

    clusters = {}
    # il = 0
    unclustered_seqs = [ d['unique_id'] for d in viterbi_info ]

    for line in viterbi_info:
        # first make a list of all the rearrangement parameters
        # if (il%5) == 0:
        #     print '---'
        # il += 1
        # hashlist = [ line[region + '_gene'].replace(re.findall('\*[0-9][0-9]', line[region + '_gene'])[0], '') for region in utils.regions ]
        hashlist = [ line[region + '_gene'] for region in utils.regions ]
        hashlist += [ str(line[d + '_del']) for d in utils.real_erosions ]
        hashlist += [ str(len(line[b + '_insertion'])) for b in utils.boundaries ]
        print(':'.join(hashlist))
        # then hash it
        clusterhash = utils.uidhashstr(':'.join(hashlist))
        if clusterhash in clusters:
            clusters[clusterhash].append(line['unique_id'])
        else:
            clusters[clusterhash] = [line['unique_id'], ]

    for cluster_id, namelist in clusters.items():
        print(cluster_id, namelist)
    print(len(clusters))
