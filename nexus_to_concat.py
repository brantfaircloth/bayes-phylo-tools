#!/usr/bin/env python
# encoding: utf-8
"""
beast_concatter.py

Created by Brant Faircloth on 2011-03-05.
Copyright (c) 2011 Brant Faircloth. All rights reserved.
"""

import os
import sys
import pdb
import optparse
import cPickle
from Bio import AlignIO
from tools.align.concatenate import ConcatenatedAlignment
from collections import OrderedDict

def interface():
    '''Command-line interface'''
    usage = "usage: %prog [options]"

    p = optparse.OptionParser(usage)

    p.add_option('--models', dest = 'models', action='store', 
type='string', default = None, help='The path to the models file.', 
metavar='FILE')

    p.add_option('--aligns', dest = 'aligns', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    p.add_option('--concat', dest = 'concat', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    p.add_option('--metadata', dest = 'metadata', action='store', 
type='string', default = None, help='The path to the nexus files.', 
metavar='FILE')

    (options,arg) = p.parse_args()
    if not options.models or not options.aligns:
        p.print_help()
        sys.exit(2)
    return options, arg

def get_loci_and_models(infile):
    loci = OrderedDict()
    for line in open(infile, 'rU'):
        ls = line.strip().split('\t')
        group = loci.setdefault(ls[1], OrderedDict())
        group[ls[0]] = None
    return loci

def concatenate(metadata, aligns):
    concat = ConcatenatedAlignment()
    start = 1
    for model in metadata:
        for locus in metadata[model]:
            align_file = os.path.join(aligns, "{0}.nex".format(locus))
            align = AlignIO.read(open(align_file), "nexus")
            #pdb.set_trace()
            end = start + align.get_alignment_length()
            metadata[model][locus] = (start, end)
            concat.add(align)
            start = end
    return concat, metadata

def save_concat_align(concat, outfile, format = "nexus"):
    o = open(outfile, 'w')
    o.write(concat.format(format))
    o.close()

def save_concat_metadata(metadata, outfile):
    o = open(outfile, 'w')
    cPickle.dump(metadata, o)
    o.close()

def main():
    options, args = interface()
    metadata = get_loci_and_models(options.models)
    start = 1
    concat, metadata = concatenate(metadata, options.aligns)
    save_concat_align(concat, options.concat)
    save_concat_metadata(metadata, options.metadata)
    #pdb.set_trace()

if __name__ == '__main__':
    main()
