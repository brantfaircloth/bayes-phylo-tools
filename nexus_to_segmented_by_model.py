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
from align.concatentate import ConcatenatedAlignment


def get_loci_and_models(infile):
    loci = {}
    for line in open(infile, 'rU'):
        ls = line.strip().split('\t')
        loci.setdefault(ls[1], []).append(ls[0])
    return loci

def build_concat_dict(loci):
    concats = {}
    for k in loci:
        concats[k] = ConcatenatedAlignment()
    return concats

def get_and_concat_alignments(loci, concats, aligns):
    for model in concats:
        # get the loci for a particular model
        model_loci = loci[model]
        for locus in model_loci:
            # generate the alignment name/location
            align_input = os.path.join(aligns, "{0}.nex".format(locus))
            # get the alignments for the locus
            align = AlignIO.read(open(align_input), "nexus")
            concats[model].concat(align)
    return concats

def write_concat_alignments(concats, output):
    for model in concats:
        align_name = "{0}_concat.nex".format(model)
        align_output = os.path.join(output, align_name)
        o = open(align_output, 'w')
        o.write(concats[model].format("nexus"))
        o.close()

def write_concat_model_file(concats, output, model_file_name):
    model_file_output = os.path.join(output, model_file_name)
    o = open(model_file_output, 'w')
    for model in concats:
        align_name = "{0}_concat".format(model)
        o.write("{0}\t{1}\n".format(align_name, model))
    o.close()
        
def main():
    #models = '183_loci/183_sub_models.txt'
    #aligns = '183_loci/183/'
    #output = '183_loci/183_partitioned'
    models = '917_loci/917_sub_models.txt'
    aligns = '/Users/bcf/Git/brant/seqcap/Manuscripts/Tetrapods/Nature/Data/917Loci_19Species/alignments_from_loci/917Loci_19Species_alignments_nexus/'
    output = '917_loci/917_partitioned'
    loci = get_loci_and_models(models)
    concats = build_concat_dict(loci)
    concats = get_and_concat_alignments(loci, concats, aligns)
    # create the output dir if not exists
    if not os.path.exists(output):
        os.mkdir(output)
    write_concat_alignments(concats, output)
    #write_concat_model_file(concats, output, "183_partitioned_sub_models.txt")
    write_concat_model_file(concats, output, "917_partitioned_sub_models.txt")


if __name__ == '__main__':
    main()

