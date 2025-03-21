import re
from typing import List, Dict, Tuple, Optional, Any, Union
from intervaltree import Interval

class Read:
    def __init__(self):
        self.qseqid = None
        self.qlen = None
        self.sseqid = None
        self.qstart = None
        self.qend = None
        self.sstart = None
        self.send = None
        self.btop = None

        self.btopl = None
        self.binread = None

        self.donors = []
        self.acceptors = []
        self.weights = []

    def to_interval(self):
        return Interval(min(self.sstart,self.send),max(self.sstart,self.send),self)

    def from_line(self,line):
        self.qseqid,self.qlen,self.sseqid,self.qstart,self.qend,self.sstart,self.send,self.btop = line.strip().split("\t")
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)

        self.parse_btop()
        self.btop_to_list()

    def is_reversed(self):
        return self.sstart > self.send
    
    def reverse(self):
        self.sstart,self.send = self.send,self.sstart
        self.qstart,self.qend = self.qlen-self.qend+1,self.qlen-self.qstart+1
        self.btopl = self.btopl[::-1]
        self.btop_to_list()
        self.donors = [[self.qlen-x[0],x[1],x[2]] for x in self.donors]
        self.acceptors = [[self.qlen-x[0],x[1],x[2]] for x in self.acceptors]
        self.weights = [[self.qlen-x[0],x[1],x[2]] for x in self.weights]

    def parse_btop(self):
        self.btopl = []
        pattern = re.compile(r'(\d+)|([A-Za-z-]{2})')
        matches = re.finditer(pattern, self.btop)

        for match in matches:
            if match.group(1):
                self.btopl.append(int(match.group(1)))
            else:
                self.btopl.append(match.group(2))

    def btop_to_list(self):
        self.binread = [0] * self.qlen
        index = self.qstart-1

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(index,index+b,1):
                    self.binread[i] += 1
                index += b
            elif isinstance(b, str):
                if b[0]=="-": # insertion
                    # if insertion - decrement the score of the next element
                    # this is equivalent to saying, by that position the score on the reference would have decreased by this much due to insertion penalty
                    self.binread[index] -= 1
                elif b[1]=="-": # deletion
                    self.binread[index] -= 1
                    index += 1
                else: # mismatch
                    self.binread[index] = 0
                    index += 1

    def get_sites(self,sites): # start position is on the template this time
        res = []

        # handle the case when the site is upstream of the read start position
        for i in range(self.qstart-2,-1,-1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = self.sstart - (self.qstart-1 - i)
            else:
                genome_pos = self.sstart + (self.qstart-1 - i)
            if genome_pos in sites:
                res.append([i,sites[genome_pos]])
        
        index_read = self.qstart-1
        index_genome = self.sstart
        inc = 1 if self.sstart < self.send else -1

        if index_genome in sites:
            res.append([index_read,sites[index_genome]])

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0,b,1):
                    index_genome += inc
                    index_read += 1
                    if index_genome in sites:
                        res.append([index_read,sites[index_genome]])
            elif isinstance(b, str):
                if b[0]=="-": # insertion in read
                    index_genome += inc
                elif b[1]=="-": # deletion from read
                    index_read += 1
                else: # mismatch - treat as a regular match here
                    index_genome += inc
                    index_read += 1
                    if index_genome in sites:
                        res.append([index_read,sites[index_genome]])

        # handle the case when the site is downstream of the read end position
        for i in range(index_read,self.qlen,1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = index_genome + (i - index_read)
            else:
                genome_pos = index_genome - (i - index_read)
            if genome_pos in sites:
                res.append([i,sites[genome_pos]])
            if genome_pos in sites:
                res.append([i,sites[genome_pos]])

        return res
    
    def load_donors(self,donors):
        if self.sseqid in donors:
            self.donors = self.get_sites(donors[self.sseqid])
    def load_acceptors(self,acceptors):
        if self.sseqid in acceptors:
            self.acceptors = self.get_sites(acceptors[self.sseqid])

    def load_weights(self,weights):
        # iterates over the positions on the read and extracts weights that belong to the read

        # handle the case when the insertion is upstream of the read start position
        for i in range(self.qstart-1,-1,-1):
            genome_pos = self.sstart - (self.qstart - i)
            weight = weights.get((self.sseqid,genome_pos),None)
            if weight is not None:
                self.weights.append((i,weight))

        index_read = self.qstart-1
        index_genome = self.sstart

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0,b,1):
                    weight = weights.get((self.sseqid,index_genome+i),None)
                    if weight is not None:
                        self.weights.append((index_read+i,weight))
                index_genome += b
                index_read += b

            elif isinstance(b, str):
                if b[0]=="-": # insertion in read
                    weight = weights.get((self.sseqid,index_genome),None)
                    if weight is not None:
                        self.weights.append((index_read,weight))
                    index_genome += 1

                elif b[1]=="-": # deletion from read
                    index_read += 1
                else: # mismatch - treat as a regular match here
                    weight = weights.get((self.sseqid,index_genome),None)
                    if weight is not None:
                        self.weights.append((index_read,weight))
                    index_genome += 1
                    index_read += 1

        # handle the case when the insertion is downstream of the read end position
        for i in range(index_read,self.qlen,1):
            genome_pos = index_genome + (i - index_read)
            weight = weights.get((self.sseqid,genome_pos),None)
            if weight is not None:
                self.weights.append((i,weight))

        return
    
    def read2genome(self,pos):
        for i in range(self.qstart-2,-1,-1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = self.sstart - (self.qstart-1 - i)
            else:
                genome_pos = self.sstart + (self.qstart-1 - i)
            if i == pos:
                return genome_pos
        
        index_read = self.qstart-1
        index_genome = self.sstart
        inc = 1 if self.sstart < self.send else -1

        if index_read == pos:
            return index_genome

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0,b,1):
                    index_genome += inc
                    index_read += 1
                    if index_read == pos:
                        return index_genome
            elif isinstance(b, str):
                if b[0]=="-": # insertion in read
                    index_genome += inc
                elif b[1]=="-": # deletion from read
                    index_read += 1
                else: # mismatch - treat as a regular match here
                    index_genome += inc
                    index_read += 1
                    if index_read == pos:
                        return index_genome

        # handle the case when the site is downstream of the read end position
        for i in range(index_read,self.qlen,1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = index_genome + (i - index_read)
            else:
                genome_pos = index_genome - (i - index_read)
            if i == pos:
                return genome_pos