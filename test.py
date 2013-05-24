#!/usr/bin/env python

from vdjalign import bwa

index = bwa.open_bwa_index('data/ighv.fasta')
print index

print index.align_sequence('CCTGGACAAGGGCTTGAGTGGATGGGACGGATCAACCCTAACAGTGGTGGCACAAACTAT')
