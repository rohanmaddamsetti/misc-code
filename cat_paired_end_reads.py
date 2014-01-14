#!/usr/bin/python

## cat_paired_end_reads.py by Rohan Maddamsetti

## This script takes directories containing paired end reads as input.
## It catenates paired end reads from different lanes into pairs of files
## because this script preserves paired end information.
## The output of these system cat calls go into a directory specified for
## output.

import argparse
from os import listdir, makedirs, chdir, getcwd, waitpid
from os.path import exists, abspath, basename, normpath, join
from subprocess import Popen

## Stack overflow solution for removing duplicates from a list.
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def main():
    parser = argparse.ArgumentParser(description='catenate paired-end reads')
    parser.add_argument('-i', dest='input', metavar='N', type=str, nargs='+',
                        help='directories containing reads')
    parser.add_argument('-o', dest='output_dir', metavar='output', type=str,
                        nargs=1, help='output directory')
    
    args = parser.parse_args()
    args.input = [normpath(abspath(x)) for x in args.input]
    args.output_dir = normpath(abspath(args.output_dir.pop()))
    print args.input
    print args.output_dir

    if not exists(args.output_dir):
        makedirs(args.output_dir)

    cat_args = []

    fastq_files = []
    for i in args.input:
        for x in listdir(i):
            if x.endswith(".fastq"):
                fastq_files.append(join(i,x))

    prefixes = f7([basename(x[:-18]) for x in fastq_files])
    matches = []
    for x in prefixes:
        matches = [i for i in fastq_files if basename(i).startswith(x)]
        cat_pair1 = [i for i in matches if basename(i).endswith('R1_001.fastq')]
        if len(cat_pair1) > 1:
            cat_args.append(cat_pair1)
        cat_pair2 = [i for i in matches if basename(i).endswith('R2_001.fastq')]
        if len(cat_pair2) > 1:
            cat_args.append(cat_pair2)
    for x in cat_args:
        prefix = basename(x[0][:-18])
        suffix = None
        if x[0].endswith('R1_001.fastq'):
            suffix = '_combined_R1.fastq'
        elif x[0].endswith('R2_001.fastq'):
            suffix = '_combined_R2.fastq'
        the_call = ["cat"] + x + [">", join(args.output_dir, prefix+suffix)]
        print ' '.join(the_call)
        print
#        p = Popen(the_call, shell=True)
        
main()
