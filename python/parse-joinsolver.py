#!/usr/bin/env python

import utils

import xml.etree.ElementTree as ET
tree = ET.parse('/home/dralph/multijoin.xml')
root = tree.getroot()

for query in root:
    print query.tag, query.attrib
    query_seq = query.find('vmatches').find('userSeq').find('bases').text.replace(' ','')
    print query_seq
    vmatches = [ match for match in query.find('vmatches') if match.tag == 'germline' ]
    for match in vmatches:
        print utils.color_gene(match.attrib['id'])
    # for subchild in query:
    #     print '    ', subchild.tag, subchild.attrib
    #     for subsubchild in subchild:
    #         print '      ', subsubchild.tag, subsubchild.attrib, subsubchild.text
    #         for subsubsubchild in subsubchild:
    #             print '        ', subsubsubchild.tag, subsubsubchild.attrib, subsubsubchild.text

    break
