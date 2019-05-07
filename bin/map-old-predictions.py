#!/usr/bin/env python2

import sys

#class hit(object):
 #   def __init__(self, line):
#        self.qid, self.status, self.sid, self.pid, self.pcov = line.strip().split("\t")
#        self.score = float(self.pid) * float(self.pcov)
class hit(object):
    def __init__(self, line):
        rec = line.strip().split('\t')
        self.qid = rec[0]
        self.sid = rec[1]
        self.pid = float(rec[2])
        self.qcov = float(rec[12])
        self.score = self.pid * self.qcov


class transcript(object):
    def __init__(self):
        self.hits = {} # key=status, value=(id, pid, pcov)

    def add_hit(self, line):
        this_hit = hit(line)
        self.hits[this_hit.sid] = this_hit

    def get_best(self):
        #best_hit = ''
        def get_score(hit):
            return hit.score

        if len(self.hits.keys()) > 0:
            best_hit = sorted(self.hits.values(), key=lambda x: get_score(x), reverse=True)[0]

        return best_hit
#    def get_best(self):
#        best_hit = ''
#        for rec in self.hits:
#            if best_hit == '':
#                best_hit = self.hits[rec]
#            else:
#                if hits[rec].score > best_hit.score:
#                    best_hit = self.hits[rec]
#
#        return best_hit


def main():
#    results = sys.argv[1] # e.g.expanded-gnavigator-results.tsv
    results = '/projects/bullfrog_assembly/genome/ARCS/annotation/fresh-Oct2018/shared-with-uvic/compare-proteins/v2-v2-v3.blastp.tsv' 
    transcripts = {}
    with open(results, 'r') as infile:
        for line in infile:
            this_line = hit(line)
            if this_line.qid not in transcripts:
                transcripts[this_line.qid] = transcript()

            transcripts[this_line.qid].add_hit(line)
    for trans in transcripts:
        result = transcripts[trans].get_best()
        if result.pid == 100.0:
            # only want cases where the v3 prediction has a perfect overlap with a v2 one
            print(result.qid, result.sid, result.score)


if __name__ == '__main__':
    main()
