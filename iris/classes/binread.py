import copy
from typing import List, Dict, Tuple, Optional, Any, Set, Union
from ..classes.read import Read

class Binread:
    def __init__(self):
        self.read1 = Read()
        self.read2 = Read()

        self.sj1 = None
        self.sj2 = None

        self.binread = []
        self.breakpoint = None
        self.score = None
        self.genes = None

        self.orientation = None # value 1 means that in the RIS read1 is on the left, read2 is on the right. Value 2 means the opposite.
        self.reversed = False

    def __str__(self):
        genome1_read_breakpoint = self.breakpoint-1 if self.orientation==1 else self.breakpoint
        genome2_read_breakpoint = self.breakpoint if self.orientation==1 else self.breakpoint-1

        genome1_pos = self.read1.read2genome(genome1_read_breakpoint)
        genome2_pos = self.read2.read2genome(genome2_read_breakpoint)
        gene1 = self.sj1[0] if self.sj1 is not None else "-"
        gene2 = self.sj2[0] if self.sj2 is not None else "-"

        # compute start and end positions for each genome
        genome1_read_start =  self.read1.qstart if self.orientation == 1 else genome1_read_breakpoint
        genome1_read_end = genome1_read_breakpoint if self.orientation == 1 else self.read2.qstart
        genome2_read_start = genome2_read_breakpoint if self.orientation == 1 else self.read1.qend
        genome2_read_end = self.read2.qend if self.orientation == 1 else genome1_read_breakpoint

        genome1_start_pos = self.read1.read2genome(genome1_read_start) if self.orientation == 1 else self.read2.read2genome(genome1_read_start)
        genome1_end_pos = self.read1.read2genome(genome1_read_end) if self.orientation == 1 else self.read2.read2genome(genome1_read_end)
        genome2_start_pos = self.read2.read2genome(genome2_read_start) if self.orientation == 1 else self.read1.read2genome(genome2_read_start)
        genome2_end_pos = self.read2.read2genome(genome2_read_end) if self.orientation == 1 else self.read1.read2genome(genome2_read_end)

        # parse genes. for every gene if both exonic and intronic labels exist - select exonic
        genes1 = dict()
        for gene in sorted(list(self.genes[0])):
            if gene == "-":
                continue
            genes1.setdefault(gene[0],"")
            if gene[1] == "exon":
                genes1[gene[0]] = "exon"
            elif gene[1] == "transcript" and genes1[gene[0]] != "exon":
                genes1[gene[0]] = "intron"
            else:
                continue
        genes1 = [x[0]+":"+x[1] for x in genes1.items()]
        genes1 = ";".join(genes1) if len(genes1)>0 else "-"

        genes2 = dict()
        for gene in sorted(list(self.genes[1])):
            if gene == "-":
                continue
            genes2.setdefault(gene[0],"")
            if gene[1] == "exon":
                genes2[gene[0]] = "exon"
            elif gene[1] == "transcript" and genes2[gene[0]] != "exon":
                genes2[gene[0]] = "intron"
            else:
                continue
        genes2 = [x[0]+":"+x[1] for x in genes2.items()]
        genes2 = ";".join(genes2) if len(genes2)>0 else "-"

        
        return "\t".join([self.read1.qseqid,
                          str(genome1_read_breakpoint),
                          str(genome2_read_breakpoint),
                          str(genome1_read_start),
                          str(genome1_read_end),
                          str(genome2_read_start),
                          str(genome2_read_end),
                          self.read1.sseqid,
                          self.read2.sseqid,
                          str(genome1_pos),
                          str(genome2_pos),
                          str(genome1_start_pos),
                          str(genome1_end_pos),
                          str(genome2_start_pos),
                          str(genome2_end_pos),
                          str(self.orientation),
                          str(int(self.reversed)),
                          str(self.score),
                          gene1,
                          gene2,
                          genes1,
                          genes2,
                          "".join([str(x) for x in self.binread])])

    def add_read1(self,line):
        self.read1.from_line(line)
        assert self.read2.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read2.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."
    def add_read2(self,line):
        self.read2.from_line(line)
        assert self.read1.qlen is None or self.read1.qlen == self.read2.qlen, "Read lengths do not match."
        assert self.read1.qseqid is None or self.read1.qseqid == self.read2.qseqid, "Read names do not match."

    def is_reversed(self):
        return self.read1.is_reversed() == self.read2.is_reversed() and self.read1.is_reversed()
    
    def reverse(self):
        self.read1.reverse()
        self.read2.reverse()

    @staticmethod
    def _find_breakpoint(list1, list2, orientation):
        max_sum = float('-inf')
        max_indices = []
        for i in range(len(list1)):
            current_sum = sum(list1[:i]) + sum(list2[i:])
            if current_sum > max_sum:
                max_sum = current_sum
                max_indices = [i]
            elif current_sum == max_sum:
                max_indices.append(i)

        breakpoint = sum(max_indices) // len(max_indices) if max_indices else -1
        binread = list1[:breakpoint]+list2[breakpoint:] if orientation==1 else list2[:breakpoint]+list1[breakpoint:]
        score = sum(binread)/len(binread)
        return [breakpoint,binread,score]
    
    def break_at(self,pos,orientation):
        # returns the result if read is broken at specific position
        # can be used to score breaks at junction sites by explicitely splitting at a given position
        # and incrementing the returned score by the wight of the junction
        binread = None
        if orientation == 1:
            binread = self.read1.binread[:pos+1]+self.read2.binread[pos+1:] # +1 because we want to include the breakpoint
        else:
            binread = self.read2.binread[:pos+1]+self.read1.binread[pos+1:]
        score = sum(binread)/len(binread)
        return [pos,binread,score]

    def find_breakpoint(self,weight_pair,orientation):
        results = []
        
        wp1 = weight_pair[0]
        wp2 = weight_pair[1]

        self.sj1,self.sj2 = None,None

        if wp1[0] is None and wp2[0] is None:
            pos,binread,score = self._find_breakpoint(self.read1.binread, self.read2.binread,orientation)
            results.append([pos,binread,score,(None,None)])
        
        if wp1[0] is not None and wp1[0]<len(self.read1.binread):
            if orientation == 1:
                self.read1.binread[wp1[0]] = wp1[2]
            else:
                self.read2.binread[wp1[0]] = wp1[2]
        if wp2[0] is not None and wp2[0]<len(self.read2.binread):
            if orientation == 1:
                self.read2.binread[wp2[0]] = wp2[2]
            else:
                self.read1.binread[wp2[0]] = wp2[2]

        if wp1[0] is not None:
            bp = self.break_at(wp1[0],orientation)
            bp.append((wp1[1],wp2[1]))
            results.append(bp)
            self.sj1,self.sj2 = [bp[3][0],bp[3][1]]
        if wp2[0] is not None:
            bp = self.break_at(wp2[0],orientation)
            bp.append((wp1[1],wp2[1]))
            results.append(bp)
            self.sj1,self.sj2 = [bp[3][0],bp[3][1]]

        if len(results)==0:
            self.breakpoint,self.binread,self.score,self.genes = None,None,None,None
        else:
            self.breakpoint,self.binread,self.score,self.genes = max(results, key=lambda x: x[2])
    