#!/usr/bin/env python2

import sys
import re

from gfftools import GFF


def selector(gff_record, hit_mrna, hit_gene):
    # hit_mrna is set of desired transcript IDs
    # hit_gene is set of desired gene IDs

    if gff_record.type == 'mRNA':
        if gff_record.id in hit_mrna:
            outbuff = gff_record.raw
        else:
            outbuff = ''
    elif gff_record.type == 'gene':
        if gff_record.id in hit_gene:
            outbuff = gff_record.raw
        else:
            outbuff = ''
    elif gff_record.type == 'CDS' or gff_record.type == 'exon':
        for par in gff_record.parent:
            if par in hit_mrna:
                # drop absent parents
                outbuff = drop_par(gff_record,hit_mrna)
                break
            else:
                outbuff = ''
    else:
        outbuff = ''

    return outbuff


def drop_par(gff_record, hit_mrna):
    raw = gff_record.raw
    for par in gff_record.parent:
        if par not in hit_mrna:
            raw = re.sub(par,'',raw)
    # drop any repeated commas due to deletions
    raw = re.sub(',+',',',raw)
    raw = re.sub('Parent=,','Parent=',raw)
    raw = re.sub(',$','',raw)

    return raw


def load_hits(records):
    # records is file with one mRNA ID per line

    hit_mrna = set()
    hit_gene = set()

    for rec in records:
        mrna = rec.strip()
        gene = re.sub('-mRNA-[0-9+]','',mrna)
        
        hit_mrna.add(mrna)
        hit_gene.add(gene)

    return (hit_mrna, hit_gene)


def main():
    gff_file = sys.argv[1]
    trans = sys.argv[2]

    with open(trans, 'r') as infile:
        hit_mrna, hit_gene = load_hits(infile)

    with open(gff_file, 'r') as gff:
        print '##gff-version 3'
        for line in gff:
            if line[0] == '#':
                continue
            gff_line = GFF(line)
            res = selector(gff_line, hit_mrna, hit_gene)
            # res is an empty string if there was no match
            if res:
                sys.stdout.write(res)


if __name__ == '__main__':
    main()

### EOF ###
