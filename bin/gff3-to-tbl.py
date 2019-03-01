#!/usr/bin/env python2

# usage: gff3-to-tbl.py genes.gff genome.fa alignment.blastp first_locus_number

import sys
import re
import copy

import argparse

from fastatools.fastatools import fasta_iter as f
from gfftools.gfftools import *

############
# parse arguments

desc = "".join(["Convert a MAKER2 gff to NCBI's tbl format.",
    " Also make adjustments to conform to NCBI's submission requirements,",
    " namely miniumum intron width of 10 bp. This is performed by sliding the",
    " position of the 3' intron boundary rightwards in increments of 3 bp.",
    " so as not to introduce frameshift errors."])

parser = argparse.ArgumentParser(description=desc)

parser.add_argument('gff', action='store', help='MAKER2 gff file')
parser.add_argument('fasta', action='store',
    help='FASTA file of scaffolds annoated by MAKER2')
#parser.add_argument('aln', action='store', help='Tabular blastp alignment of MAKER2 proteins vs. SwissProt')

parser.add_argument('--locus_tag', '-t', action='store', help='NCBI-supplied locus tag prefix [AB205]', default="AB205")
parser.add_argument('--locus_tag_start', '-s', action='store', help='Initial value to use for locus_tag generation [50]', default=50)
parser.add_argument('--tag_width', '-w', action='store', help='Width for numeric portion of locus_tag [7]', default=7)
parser.add_argument('--tag_jump', '-j', action='store', help='Interval between locus_tag numbers [10]', default=10)
parser.add_argument('--min_ident', '-i', action='store', help='Minimum percent identity to accept alignment for annotation [25]', default=25)
parser.add_argument('--min_cov', '-c', action='store', help='Minimum percent coverage to accept alignment for annotation [50]', default=50)
parser.add_argument('--ref_break', '-b', action='store', help='String to split reference protein names (expects SwissProt-style FASTA deflines [OS=]', default='OS=')
parser.add_argument('--prod_break', '-r', action='store', help='Field separator for product line [space]', default=' ')
parser.add_argument('--prefix', '-p', action='store', help='Prefix to use for product name [similar to]', default='similar to')
parser.add_argument('--unknown', '-u', action='store', help='Label to give proteins without an acceptable annotation [hypothetical protein]', default='hypothetical protein')

args = parser.parse_args()

gff = args.gff
fasta = args.fasta
#aln = args.aln

LOCUS = args.locus_tag
LOCS = int(args.locus_tag_start)
LOCW = int(args.tag_width)
LOCJ = int(args.tag_jump)
MIN_IDENT = int(args.min_ident)
MIN_COV = int(args.min_cov)
BRK = args.ref_break
SEPO = args.prod_break
PRE = args.prefix
UNK = args.unknown

############
# hard-coded declarations

# alphabet to use for isoform labeling
ISOALPHA=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S',
            'T','U','V','W','X','Y','Z']

#gffout=open(re.sub(".gff","-prepared.gff",gff),"w")
tblout=open(re.sub(".gff","-prepared.tbl",gff),"w")
#lokey=open(re.sub(".gff","-locus_tag-conversion-key.tsv",gff),"w")

###################
# MAIN

genes = {}
# read in the entries from the gff file
with open(gff,"r") as infile:
    print "Reading the annotations"
    for line in infile:
        if line[0] == "#" or line[0] == "-":
            continue

        feature = GFF(line)

        if feature.type == "gene":
            genes[feature.id] = Gene(feature)
            genes[feature.id].locus_tag = feature.locus_tag

        elif feature.type == "mRNA" or feature.type == "ncRNA":
            # transcripts can only have one parent, so access it explicitly
            genes[feature.parent[0]].add_transcript(feature)
            genes[feature.parent[0]].transcript[feature.id].locus_tag = feature.locus_tag
            genes[feature.parent[0]].transcript[feature.id].product = feature.product
            if feature.prot_desc:
                genes[feature.parent[0]].transcript[feature.id].prot_desc = feature.prot_desc

        elif feature.type == "exon":
            if feature.name:
                genes[feature.name].transcript[feature.parent[0]].add_exon(feature)
                continue
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_exon(feature)

        elif feature.type == "CDS":
            if feature.name:
                genes[feature.name].transcript[feature.parent[0]].add_exon(feature)
                continue
            for parent in feature.parent:
                for gene in genes:
                    if re.findall(genes[gene].id,parent):
                        genes[gene].transcript[parent].add_cds(feature)


###################
# Collect annotations for each transcript
#aseen = set()
#annots = {}
#with open(aln,"r") as blast:
#    print "Loading the annotations"
#    for rec in blast:
#        thisaln = rec.split("\t")
#        if thisaln[0] not in aseen:
#            cov = thisaln[-1]
#            if thisaln[2] >= MIN_IDENT and cov >= MIN_COV:
#                thisannot = thisaln[-2].split(" ")
#                refnam = ''
#                for i in range(0,len(thisannot)):
#                    if i == 0 and len(re.findall(BRK,thisannot[i])) < 1:
#                        refnam += thisannot[i]
#                    elif i > 0 and len(re.findall(BRK,thisannot[i])) < 1:
#                        refnam += SEPO+thisannot[i]
#                    else:
#                        break
#                annots[thisaln[0]] = PRE + SEPO + refnam
#            aseen.add(thisaln[0])

####################
# Add locus_tag to each gene and annotation to each transcript
scaf_order = {}
print "Ordering the predictions along the genomic scaffolds"
for entry in genes:
    this_gene = genes[entry]
    if this_gene.gene.seqid not in scaf_order:
        if this_gene.gene.strand == "-":
            scaf_order[this_gene.gene.seqid] = [[this_gene.gene.id,this_gene.gene.end]]
        else:
            scaf_order[this_gene.gene.seqid] = [[this_gene.gene.id,this_gene.gene.start]]
    else:
        if this_gene.gene.strand == "-":
            scaf_order[this_gene.gene.seqid].append([this_gene.gene.id,this_gene.gene.end])
        else:
            scaf_order[this_gene.gene.seqid].append([this_gene.gene.id,this_gene.gene.start])

#locflag = 0 # flag to indicate first instance of locus ID
#lokey_dict = {}
#print "Adding locus tags and applying functional annotations"
#for scaf in scaf_order:
#    for entry in sorted(scaf_order[scaf],key = lambda x: int(x[1])):
#        gene_name = entry[0]
#        if locflag == 0:
#            nam = LOCUS + "_" + str(LOCS).zfill(LOCW)
#            locs = LOCS + LOCJ
#            locflag = 1
#        else:
#            nam = LOCUS + "_" + str(locs).zfill(LOCW)
#            locs += LOCJ
#
#        genes[gene_name].gene.locus_tag = nam
#        lokey.write(gene_name + "\t" + nam + "\n")
#        lokey_dict[gene_name] = nam
#        # Add locus_tags to each gene's transcripts, with addition of isoform
#        #  letter if needed
#        if len(genes[gene_name].transcript) == 1:
#            tr_id = genes[gene_name].transcript.keys()[0]
#            genes[gene_name].transcript[tr_id].transcript.locus_tag = nam
#            lokey.write(tr_id + "\t" + nam + "\n")
#
#            if genes[gene_name].transcript[tr_id].transcript.type == "mRNA":
#                if tr_id in annots:
#                    genes[gene_name].transcript[tr_id].transcript.product = annots[tr_id]
#                else:
#                    genes[gene_name].transcript[tr_id].transcript.product = UNK
#        else:
#            loc_iter = 0
#            for prod in sorted(genes[gene_name].transcript):
#                genes[gene_name].transcript[prod].transcript.locus_tag = nam + ISOALPHA[loc_iter]
#                lokey.write(prod + "\t" + nam + ISOALPHA[loc_iter] + "\n")
#
#                if genes[gene_name].transcript[prod].transcript.type == "mRNA":
#                    if prod in annots:
#                        genes[gene_name].transcript[prod].transcript.product = annots[prod]
#                    else:
#                        genes[gene_name].transcript[prod].transcript.product = UNK
#                loc_iter += 1
#lokey.close()

###################
#load genomic scaffolds
genome = {}
with open(fasta,"r") as file_object:
    print "Loading genomic scaffolds"
    for line in f(file_object):
        genome[line[0]] = line[1]

###################
# output tbl lines
#scafs = set()
prev_scaf = ''
# order genes by ID; hopefully this corresponds to their position on the scaffold
print "Writing final .tbl file"
for scaf in scaf_order:
    for entry in sorted(scaf_order[scaf], key=lambda x: int(x[1])):
        rec = entry[0]
        if genes[rec].gene.seqid != prev_scaf:
            tblout.write("".join([">Feature ",genes[rec].gene.seqid,"\n"]))
            tblout.write("\t".join(["1",str(len(genome[genes[rec].gene.seqid])),
                                            "REFERENCE\n"]))
            tblout.write("\t".join(["\t\t\tPBARC","12345\n"]))

        tblout.write(str(genes[rec].print_gene()))
        prev_scaf = genes[rec].gene.seqid

        # sort each gene's transcripts by name
        for prod in sorted(genes[rec].transcript.keys()):
#            try:
                # if there is a note about a split feature, fix it
#                if genes[rec].transcript[prod].transcript.note:
#                    this_note = genes[rec].transcript[prod].transcript.note.split(" ")
#                    fixed_note = []
#                    for i in this_note:
#                        if i in lokey_dict:
#                            fixed_note.append(lokey_dict[i])
#                        else:
#                            fixed_note.append(re.sub("end_","end;",i))
#                    genes[rec].transcript[prod].transcript.note = " ".join(fixed_note)
            tblout.write(str(genes[rec].transcript[prod].print_transcript( ### indent one right
                        product_type = genes[rec].transcript[prod].transcript.type,
                        outform = 'tbl',sequence=genome[genes[rec].gene.seqid])))
            #except:
            #    print "Failed to write out " + str(genes[rec].transcript[prod].transcript.id)

tblout.close

### EOF ###
