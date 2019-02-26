#!/usr/bin/env python2

import sys

def gen_tag(curr_count, width, pre):
    newnum = str(int(curr_count) + 10)
    newnum_filled = newnum.zfill(width)
    tag = '_'.join([pre, newnum_filled])

    return (newnum, tag)


def main():
    # presets
    PRE = 'AB205'
    WIDTH = 7
    LOC = 222240 # next available locus_tag suffix

    updates = 'mapped-predictions-genes-uniq2-dedup2-better-maps-locus-tagged-fixed.txt'
    #updates = 'mapped-predictions-genes-uniq2-locustags.txt'
    #updates = 'v2-locus-tags-mapped-to-v3.txt'
    all_proteins = 'version3-gene-IDs.txt'
    outname = 'version3-gene-locus-tags.tsv'

    # read in pre-mapped sequences
    update_dict = dict()
    with open(updates, 'r') as infile:
        for line in infile:
            sid, loc = line.strip().split(' ')
            update_dict[sid] = loc

    # read all sequences and either return mapped locus_tag or make a new one
    curr_tag = LOC
    with open(outname, 'w') as outfile:
        with open(all_proteins,'r') as infile:
            for line in infile:
                sid = line.strip()
                if sid not in update_dict:
                    curr_tag, tag = gen_tag(curr_tag, WIDTH, PRE)
                    outbuff = ''.join([sid,'\t',tag,'\n'])
                else:
                    outbuff = ''.join([sid, '\t', update_dict[sid], '\n'])

                outfile.write(outbuff)

if __name__ == '__main__':
    main()
