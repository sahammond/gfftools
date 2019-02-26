#!/usr/bin/env python2

import sys

import gfftools

# read in locus tag table
# read in gff
# add locus tags and protein/pfam annotations

class transcript(object):
    def __init__(self):
        self.blast_hits = {}
        self.pfam_hits = {}

    def add_blast_hit(self, line):
        this_hit = blast_hit(line)
        self.blast_hits[this_hit.sid] = this_hit

    def add_pfam_hit(self, line):
        this_hit = pfam_hit(line)
        self.pfam_hits[this_hit.sname] = this_hit

    def get_best_blast(self):
        best_hit = ''
        def get_score(hit):
            return hit.score

        if len(self.blast_hits.keys()) > 0:
            best_hit = sorted(self.blast_hits.values(), key=lambda x: get_score(x),
                              reverse=True)[0]

        return best_hit

    def get_best_pfam(self):
        best_hit = ''
        def get_score(hit):
            return hit.evalue
        if len(self.pfam_hits.keys()) > 0:
            best_hit = sorted(self.pfam_hits.values(), key=lambda x: get_score(x),
                              reverse=False)[0]

        return best_hit
    

def load_tags(tag_file):
    table = {}
    with open(tag_file, 'r') as infile:
        for line in infile:
            mrna, tag = line.strip().split('\t')
            table[mrna] = tag

    return table


def load_blast(prot_table, POS, QCOV):
    prot_aln = {}
    with open(prot_table, 'r') as infile:
        for line in infile:
            this_line = blast_hit(line)
            if this_line.qid not in prot_aln:
                prot_aln[this_line.qid] = transcript()
            if this_line.ppos > POS and this_line.qcov > QCOV:
                prot_aln[this_line.qid].add_blast_hit(line)

    return prot_aln


class blast_hit(object):
    def __init__(self, line):
        rec = line.strip().split('\t')
        self.qid = rec[0]
        self.sid = rec[1]
        self.pid = float(rec[2])
        self.ppos = float(rec[12])
        self.qcov = float(rec[13])
        self.sname = rec[-1]
        self.score = self.ppos * self.qcov


class pfam_hit(object):
    def __init__(self, line):
        self.qid, self.evalue, self.sname = line.strip().split('\t')        
        self.evalue = float(self.evalue)


def load_pfam(pfam_table, EVALUE):
    pfam_aln = {}
    with open(pfam_table, 'r') as infile:
        for line in infile:
            this_line = pfam_hit(line)
            if this_line.qid not in pfam_aln:
                pfam_aln[this_line.qid] = transcript()
            if this_line.evalue <= EVALUE:
                pfam_aln[this_line.qid].add_pfam_hit(line)

    return pfam_aln


def get_tag(gff_record, mrna_table, gene_table):
    if gff_record.type == 'gene':
        tag = gene_table[gff_record.id]
    elif gff_record.type == 'mRNA':
        tag = mrna_table[gff_record.id]
    else:
        tag = ''
    
    return tag


def load_gv(gv_file):
    gv_good = set()
    with open(gv_file, 'r') as infile:
        for line in infile:
            gv_good.add(line.strip())

    return gv_good


def main():
    gene_loci = 'version3-gene-locus-tags.tsv'
    mrna_loci = 'version3-transcript-locus-tags.tsv'
    gff = '/projects/bullfrog_assembly/genome/ARCS/annotation/fresh-Oct2018/maker-second-round/analysis/third-round/gag_output/genome.gff'
    blastp = '/projects/bullfrog_assembly/genome/ARCS/annotation/fresh-Oct2018/submission_prep/blastp-annotated2.tsv'
    pfam = 'pfam-hits.tsv'
    genevalidator = 'genevalidator-90.txt'
    POS = 25.0
    QCOV = 50.0
    EVALUE = 0.00001
    
    mrna_tags = load_tags(mrna_loci)
    gene_tags = load_tags(gene_loci)
    blast_aln = load_blast(blastp, POS, QCOV)
    pfam_aln = load_pfam(pfam, EVALUE)
    gv_good = load_gv(genevalidator)

    with open(gff, 'r') as infile:
        for line in infile:
            if line[0] == '#':
                sys.stdout.write(line)
            else:
                rec = gfftools.GFF(line)
                tag = get_tag(rec, mrna_tags, gene_tags)
                if rec.type == 'mRNA':
                    if rec.id in blast_aln:
                        hit = blast_aln[rec.id]
                        prot = hit.get_best_blast()
                    else:
                        prot = ''

                    if rec.id in pfam_aln:
                        hit = pfam_aln[rec.id]
                        dom = hit.get_best_pfam()
                    else:
                        dom = ''

                    if prot:
                        hit_name = prot.sname
                        if prot.qid in gv_good:
                            rec.add_prot_desc(hit_name, flag=True)
                        else:
                            rec.add_prot_desc(hit_name)
                    elif dom:
                        pfam_name = dom.sname
                        rec.add_pfam(pfam_name)
                    else:
                        rec.add_prot_desc('hypothetical protein', flag=True)

                    if tag:
                        rec.add_locus_tag(tag)

                    # process attributes and print out
                    rec.extend_attr()
                    print rec.print_record()

                elif rec.type == 'gene':
                    # genes only get locus_tags
                    if tag:
                        rec.add_locus_tag(tag)
                    rec.extend_attr()
                    print rec.print_record()

                else:
                    # non-gene or non-mRNA record
                    rec.extend_attr()
                    print rec.print_record()

if __name__ == '__main__':
    main()
