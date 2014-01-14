#!/usr/bin/python

## AssembleCumberbatch.py by Rohan Maddamsetti.

## Input: a directory containing paired-end fastq reads.
## Output: results of breseq and a5 go into results/reseq/ and results/denovo/.
    
## This script runs breseq as well as the a5 assembly pipeline developed by
## Jon Eisen's group.
## breseq must be in the system path.
## The location of a5 is hardcoded as a global constant.

import argparse
from os import listdir, makedirs, chdir, getcwd
from os.path import exists, abspath, basename, normpath, join
from time import sleep
from subprocess import call, Popen
import csv

## TODO: option for just breseq
## TODO: option for just A5
## TODO: get rid of these hardcoded paths.

PROCESSMAX = 8
A5_PATH = "/Users/Rohan/Bioinformatics/ngopt_a5pipeline_macOS-x64_20130326/bin/a5_pipeline.pl"

BASE_PATH = "/Users/Rohan/Desktop/IlluminaRun2013"
REF_PATH = join(BASE_PATH, "data/REL606.gbk")

## Stack overflow solution for removing duplicates from a list.
def f7(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def remove_breseq_files(base,j):
    '''change directory to the output file, then remove intermediate files'''
    curr_dir = getcwd()
    output_path = join(BASE_PATH, "results", base, "reseq", j, j)
    chdir(output_path)
    call(["rm", "-rf", "01_sequence_conversion"])
    call(["rm", "-rf", "02_reference_alignment"])
    call(["rm", "-rf", "03_candidate_junctions"])
    call(["rm", "-rf", "04_candidate_junction_alignment"])
    call(["rm", "-rf", "05_alignment_correction"])
    call(["rm", "-rf", "06_bam"])
    call(["rm", "-rf", "07_error_calibration"])
    call(["rm", "-rf", "08_mutation_identification"])
    chdir(curr_dir)
    

def update_status_hash (base, i, calltype, check, status_hash):
    if check == 0:
        if calltype == "breseq":
            status_hash[i] = status_hash[i] + 8 # set breseq success bit
            remove_breseq_files(base, i)
        elif calltype == "a5":
            status_hash[i] = status_hash[i] + 2 # set a5 success bit
    elif check == 1:
        if calltype == "breseq":
            status_hash[i] = status_hash[i] + 4 # set breseq error bit
            remove_breseq_files(base, i)
        elif calltype == "a5":
            status_hash[i] = status_hash[i] + 1 # set a5 error bit.
    return status_hash


def check_processes(base, processes, status_hash):
    ''' remove all finished processes from the list, updating the status hash
    all the while. If still full, wait for the top process to finish,
    then update the status hash.'''
    print "check_processes function is running."
    sleep(1)
    updated_processes = []
    for x in processes:
        (i, calltype, proc) = x
        check = proc.poll()
        print proc
        check = proc.poll()
        print "Checking status"
        if check is None:
            updated_processes.append(x)
        else:
            status_hash = update_status_hash(base, i, calltype, check, status_hash)
    processes = updated_processes
    if len(processes) >= PROCESSMAX:
        top = processes[0]
        (i, calltype, proc) = top
        print proc
        check = proc.poll()
        print "Checking status"
        print check
        while check is None: ## hasn't terminated yet.
            sleep(60)
            check = proc.poll()
        status_hash = update_status_hash(base, i, calltype, check, status_hash)
        processes.pop(0)
    return (processes, status_hash)

def spawn_child(base, reads, j, process_type, queue, status):
    ''' Change to the output_path, spawn the child process, add it to the queue,
        then change back to the original working directory. We have to change to the
        desired directory for A5 to work nicely.'''
    j = normpath(j)
    print "the length of the queue is:", len(queue)
    while len(queue) >= PROCESSMAX:
        (queue, status) = check_processes(base, queue, status)
 
    if process_type == "breseq":
        output_path = join(BASE_PATH, "results", base, "reseq", j)
        args = ["breseq", "-r", REF_PATH, "-o", j] + reads
    elif process_type == "a5":
        output_path = join(BASE_PATH, "results", base, "denovo", j)
        args = [A5_PATH] + reads + [j+".out"]
    if not exists(output_path):
        makedirs(output_path)
    working_dir = getcwd()
    chdir(output_path)
    print ' '.join(args)
    process = Popen(args)
    queue.append((j,process_type,process))
    chdir(working_dir)
    return(queue, status)

def get_status(statfile, genome_ids):
    '''read in status file if it exists, otherwise create an empty hash.'''
    status = dict()
    if exists(statfile):
        with open(statfile, 'r') as csvfile:
            the_reader = csv.reader(csvfile)
            for row in the_reader:
                if row == ["Genome","Status"]:
                    continue
                else:
                    print row
                    k,v = row
                    status[k] = ( int(v) & 10) # remove error bits.
    else:
        status = dict([(x,0) for x in genome_ids])
    return status

def main():
    parser = argparse.ArgumentParser(description='run breseq and a5')
    parser.add_argument('data_path', help='directory containing data', action='store')
    args = parser.parse_args()
    args.data_path = normpath(args.data_path)
    base = basename(args.data_path)

    ## default: do breseq, don't do a5.
    do_breseq = 1
    do_a5 = 0

    ## open the directory, and make an index of all paired-end reads.
    fastq_files = [x for x in listdir(args.data_path) if x.endswith(".fastq")]

    genome_ids = f7([x[:-18] for x in fastq_files])

    ## Keep track of whether the assembly programs have finished or not.
    ## This is a four-bit number. The first bit is 1 if breseq finished.
    ## The second bit is 1 if breseq had an error. The third bit is one
    ## if a5 finished. The fourth bit is one if a5 had an error.
    status_path = join(getcwd(), base + "-status.csv")
    status = get_status(status_path, genome_ids)

    ##The list contains triples of (genome_id, "breseq" or "a5", Popen obj.)
    running_processes = []
    try:
        for i in status.keys():
            paired_end_files = [x for x in fastq_files if x.startswith(i)]
            paired_end_files.sort()
            reads = [ '/'.join( [abspath(args.data_path), x] ) for x in paired_end_files]
            if do_breseq and not (int(status[i]) & 8): # breseq has not finished successfully
                running_processes, status = spawn_child(base, reads, i, "breseq", running_processes, status)
            if do_a5 and not (int(status[i]) & 2): # a5 has not finished successfully
                running_processes, status = spawn_child(base, reads, i, "a5", running_processes, status)
        while len(running_processes):
            running_processes, status = check_processes(base, running_processes, status)
    finally:
        ## Write status to a human-readable text file.
        with open(status_path, 'w') as statusfile:
            the_writer = csv.writer(statusfile)
            the_writer.writerow(["Genome","Status"])
            for k,v in status.iteritems():
                the_writer.writerow([k,v])
        ## upon startup, this file is read to prevent duplicate work in case of errors.

main()
