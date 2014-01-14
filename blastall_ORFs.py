#!/usr/bin/python

##blastall_ORFS.py by Rohan Maddamsetti.

## This script blasts all ORFs in REL606 against a given microbiome dataset.

##Usage: python blastall_ORFs.py > ../data/gut_community_hits.csv

## Questions: 

## 1) Do genes with high thetaS in REL606 tend to have more homologs
## in bacterial sequences in a gut community?

## 2) Do genes with high thetaS in REL606 tend to have more homologs in
## viral genomes in a gut community?

import sys
from subprocess import call
from os import listdir, environ, makedirs, chdir
from os.path import exists
from Bio.Blast import NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import ExactPosition

##put the path from bash into the path in this script.
mypath = environ['PATH']
sys.path.append(mypath)

def get_reference_ORFs(rel606path, orf_dir):
	''' output every ORF in REL606 into separate FASTA files stored in the
	directory named in orf_dif.'''
	## Open the directory in which the ORFs should be written to.
	## If the directory doesn't exist, create it.
	if not exists(orf_dir):
		makedirs(orf_dir)
	ref_genome = SeqIO.parse(rel606path, "genbank").next()
	for feature in ref_genome.features:
		if feature.type != "CDS":
			continue
		else:
			try:
				locus_tag = feature.qualifiers["locus_tag"][0]
			except:
				continue
 			start = feature.location.start.position
			end = feature.location.end.position
			sequence = ref_genome.seq[start:end]
			outpath = orf_dir + locus_tag + ".fasta"
			outfile = open(outpath, "w")
			outfile.write(">" + locus_tag + "\n" + str(sequence) + "\n\n")
			outfile.close()

def runBLAST(working_dir, orf_dir, out_dir, dbname, threshold="0.0000000001"):
    if not exists(out_dir):
        makedirs(out_dir)
    chdir(working_dir)
    queries = [x for x in listdir(orf_dir) if x.endswith(".fasta")]
    for cur_query in queries:
        prefix = cur_query.split('.fasta')[0]
        blast_cmd = ["blastn", "-task", "megablast", "-db", dbname,
                     "-outfmt", "5", "-evalue", threshold,
                     "-query", orf_dir + cur_query, "-out",
                     out_dir+prefix+".xml"]
        blast_cmd_string = ' '.join(blast_cmd)
        call(blast_cmd_string, shell=True)

def printBLASTresults(hit_dir):
    '''print a csv file of locus_tag, hits'''

    results = [x for x in listdir(hit_dir) if x.endswith(".xml")]
    print "locus_tag, hits"
    for cur_blast_output in results:
        tag = cur_blast_output.split('.xml')[0]
        result_handle = open(hit_dir + cur_blast_output, "r")
        blast_record = NCBIXML.parse(result_handle).next()
        alignment_count = 0
        for alignment in blast_record.alignments:
            alignment_count = alignment_count + 1
        print tag + "," + str(alignment_count)

def main():
    
    ## The following lines are commented out because they only need
    ## to be run once for the analysis.
    #get_reference_ORFs("/Users/Rohan/Desktop/REL606.gbk",
    #"/Users/Rohan/Desktop/Projects/gene_population_sizes/data/REL606_ORFs/")
    #runBLAST("/Users/Rohan/Desktop/Projects/gene_population_sizes/data/",
    #"/Users/Rohan/Desktop/Projects/gene_population_sizes/data/REL606_ORFs/",
    #"/Users/Rohan/Desktop/Projects/gene_population_sizes/data/hits/",
    #"gut_db")
    printBLASTresults("/Users/Rohan/Desktop/Projects/gene_population_sizes/data/hits/")

main()
