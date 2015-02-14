import sys

import utils

# ----------------------------------------------------------------------------------------
def cluster(viterbi_info):
    """ cluster together sequences with the same rearrangement parameters """
    clusters = {}
    for line in viterbi_info:
        # first make a list of all the rearrangement parameters
        hashlist = [ line[region + '_gene'] for region in utils.regions ]
        hashlist += [ str(line[d + '_del']) for d in utils.real_erosions ]
        hashlist += [ str(len(line[b + '_insertion'])) for b in utils.boundaries ]
        # then hash it
        clusterhash = hash(':'.join(hashlist))
        if clusterhash in clusters:
            clusters[clusterhash].append(line['unique_id'])
        else:
            clusters[clusterhash] = [line['unique_id'], ]

    for cluster_id, namelist in clusters.items():
        print cluster_id, namelist
    print len(clusters)
