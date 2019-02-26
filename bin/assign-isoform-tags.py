#!/usr/bin/env python2

import sys

from fastatools import fasta_iter


class Gene(object):
    def __init__(self, label):
        self.mrna = []
        self.name = label

    def add_tag(self, locus):
        self.tag = locus

    def add_mrna(self, sid):
        self.mrna.append(sid) 

    def get_tags(self):
        # requires a tag to have been assigned to the gene
        # and all mRNA to have been added too

        # alphabet to use for isoform labeling
        ISOALPHA = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S',
                    'T','U','V','W','X','Y','Z']

        if len(self.mrna) > 1:
            self.mrna.sort()
            letters = ISOALPHA[:len(self.mrna)]
            new_tags = [self.tag + x for x in letters]
            labels = zip(self.mrna, new_tags)
            for label in labels:
                yield label
        elif len(self.mrna) == 1:
            yield (self.mrna[0], self.tag)
        else:
            print 'no mRNA associated with %s locus_tag %s!' % (self.name, self.tag)


def main():
    # read in gene locus_tags
    tags = 'version3-gene-locus-tags.tsv'
    genes = {}
    with open(tags, 'r') as infile:
        for line in infile:
            gid, loc = line.strip().split('\t')
            genes[gid] = Gene(gid)
            genes[gid].add_tag(loc)
    # read in sequence IDs
    transcripts = '/projects/bullfrog_assembly/genome/ARCS/annotation/fresh-Oct2018/shared-with-uvic/compare-proteins/v3-proteins.fa'
    with open(transcripts, 'r') as infile:
        for sid, sqn in fasta_iter(infile):
            gid = sid.split('-mRNA')[0]
            try:
                genes[gid].add_mrna(sid)
            except KeyError:
                print '%s not in gene table!' % sid

    # write out all tags
    outname = 'version3-transcript-locus-tags.tsv'
    with open(outname, 'w') as outfile:
        for entry in genes:
            this_gene = genes[entry]
            for tag in this_gene.get_tags():
                this_mrna, this_tag = tag
                outbuff = ''.join([this_mrna, '\t', this_tag, '\n'])
                outfile.write(outbuff)


if __name__ == '__main__':
    main()
