#!/usr/bin/env python
# encoding: utf-8
"""
beast_concatter.py

Created by Brant Faircloth on 2011-03-05.
Copyright (c) 2011 __MyCompanyName__. All rights reserved.
"""

import os
import sys
import pdb
from Bio import AlignIO
from lib.helpers import ConcatenatedAlignment


def get_loci_and_models(infile):
    loci = {}
    for line in open(infile, 'rU'):
        ls = line.strip().split('\t')
        loci.setdefault(ls[1], []).append(ls[0])
    return loci

def build_super_concat(aligns, loci):
    super_concat = ConcatenatedAlignment()
    for k,v in loci.iteritems():
        align_input = os.path.join(aligns, "{0}.nex".format(v[0]))
        align = AlignIO.read(open(align_input), "nexus")
        pos = (super_concat.get_alignment_length() + 1, align.get_alignment_length(),)
        #pdb.set_trace()
        loci[k].append(pos)
        super_concat.concat(align)
    return super_concat, loci

def write_concat_alignments(super_concat, align_output):
    o = open(align_output, 'w')
    o.write(super_concat.format("nexus"))
    o.close()

def write_concat_model_file(loci, output):
    o = open(output, 'w')
    for k,v in loci.iteritems():
        o.write("{0}\t{1}\n".format(k, v[1]))
    o.close()
        
def main():
    #models = '183_loci/183_sub_models.txt'
    #aligns = '183_loci/183/'
    #output = '183_loci/183_partitioned'
    models = '917_loci/917_partitioned/917_partitioned_sub_models.txt'
    aligns = '917_loci/917_partitioned'
    output = '917_loci/917_double_partitioned/917_super_concat.nex'
    loci = get_loci_and_models(models)
    
    super_concat, loci = build_super_concat(aligns, loci)
    write_concat_alignments(super_concat, output)
    write_concat_model_file(loci, "917_loci/917_double_partitioned/917_super_concat_models_and_positions.txt")
    #write_concat_model_file(concats, output, "917_partitioned_sub_models.txt")


if __name__ == '__main__':
    main()

