#!/usr/bin/env python2

# quick and dirty gff3 (maker-style post- sort-gff.py) to bed file converter

import sys
import re

import gfftools

# n.b. s.argv[0] == script name
gff=sys.argv[1]

# dicts to hold the exon and cds records
edic={}
cdic={}

with open(gff,"r") as infile:
    for line in infile:
        if line[0] != "#":

            rec = gfftools.GFF(line)

            if rec.type == "mRNA" or rec.type == "ncRNA":
                mid=rec.id
                fstart=int(rec.start)-1 
                fend=int(rec.end)-1
                edic[mid]=[rec.seqid,fstart,fend,mid,0,rec.strand,fstart,fend,[255,0,0]]
#                if rec.type != "ncRNA":
#                    cdic[mid]=[rec.seqid,fstart,fend,mid,0,rec.strand,fstart,fend,[0,0,255]]
            elif rec.type == "exon":
                mid=re.sub(":exon","",rec.id)
                fstart=int(rec.start)-1-edic[mid][1]
                fend=int(rec.end)-1-edic[mid][1]
                flen=fend-int(fstart)+1
                if len(edic[mid]) == 9:
                    edic[mid].append(1)
                    edic[mid].append([flen])
                    edic[mid].append([fstart])
                else:
                    edic[mid][9]+=1
                    edic[mid][10].append(flen)
                    edic[mid][11].append(fstart)
            elif rec.type == "CDS":
                mid=re.sub(":cds","",rec.id)
                fstart=int(rec.start)-1-cdic[mid][1]
                fend=int(rec.end)-1-cdic[mid][1]
                flen=fend-int(fstart)+1
                if len(cdic[mid]) == 9:
                    cdic[mid][6]=int(rec.start)-1
                    cdic[mid][7]=int(rec.start)-1
                    cdic[mid].append(1)
                    cdic[mid].append([flen])
                    cdic[mid].append([fstart])
                else:
                    cdic[mid][9]+=1
                    cdic[mid][10].append(flen)
                    cdic[mid][11].append(fstart)


for exon in edic:
    # derive bed file name
    bed = re.sub(".gff","_"+str(edic[exon][3])+"_exon.bed",gff)

    with open(bed,"w") as outfile:
        outfile.write('track name="'+str(edic[exon][3])+' exon"'+' description="'+str(edic[exon][3])+'"'+' itemRgb="On"'+'\n')

        for rec in edic[exon][0:8]:
            outfile.write(str(rec)+"\t")

        for i in edic[exon][8]:
            outfile.write(str(i)+",")
        else:
            outfile.write("\t")

        outfile.write(str(edic[exon][9])+"\t")
        
        for i in edic[exon][10]:
            outfile.write(str(i)+",")
        else:
            outfile.write("\t")

        for i in edic[exon][11]:
            outfile.write(str(i)+",")
        else:
            outfile.write("\n")
if len(cdic) > 0:
    for cds in cdic:
        # derive bed file name
        bed = re.sub(".gff","_"+str(cdic[cds][3])+"_cds.bed",gff)

        with open(bed,"w") as outfile:
            outfile.write('track name="'+str(cdic[cds][3])+' CDS"'+' description="'+str(cdic[cds][3])+'"'+' itemRgb="On"'+'\n')

            for rec in cdic[cds][0:8]:
                outfile.write(str(rec)+"\t")

            for i in cdic[cds][8]:
                outfile.write(str(i)+",")
            else:
                outfile.write("\t")

            outfile.write(str(cdic[cds][9])+"\t")
            
            for i in cdic[cds][10]:
                outfile.write(str(i)+",")
            else:
                outfile.write("\t")

            for i in cdic[cds][11]:
                outfile.write(str(i)+",")
            else:
                outfile.write("\n")
### EOF ###
