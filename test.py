#!/usr/bin/env python

from vdjalign import bwa

index = bwa.load_index('data/ighv.fasta')
print index

print index.align('CCTGGACAAGGGCTTGAGTGGATGGGACGGATCAACCCTAACAGTGGTGGCACAAACTAT')
