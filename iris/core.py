import sys
import copy
import tempfile
import argparse

from .utils.common import *

from .classes.txgroup import Transcriptome
from .classes.binread import Binread
from .classes.read import Read

from typing import List, Dict, Tuple, Optional, Any, Set, Union

from intervaltree import Interval, IntervalTree

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

class IRIS:
    def __init__(self, args):         
        # INPUT FILES
        self.input1 = args.input1
        self.input2 = args.input2
        self.annotation1 = args.annotation1
        self.annotation2 = args.annotation2
        self.genome1 = args.genome1
        self.genome2 = args.genome2
        self.chim_genome = args.chim_genome

        self.output = args.output
        self.outfname_pass2 = self.output + ".corrected"
        self.outfname_grouped = self.output + ".grouped"

        self.two_pass = args.two_pass
        self.group = args.group
        self.max_dist = args.max_dist
        self.max_weight = args.max_weight
        self.full_weight = args.full_weight
        self.half_weight = args.half_weight
        self.overhang = args.overhang
        
        self.tome1 = None if self.annotation1 is None else Transcriptome()
        if self.annotation1 is not None:
            self.tome1.build_from_file(self.annotation1)
        self.tome2 = None if self.annotation2 is None else Transcriptome()
        if self.annotation2 is not None:
            self.tome2.build_from_file(self.annotation2)

        self.gene_trees1 = self.extract_genes(self.tome1)
        self.gene_trees2 = self.extract_genes(self.tome2)
        self.donors1, self.acceptors1 = self.extract_donor_acceptor(self.tome1)
        self.donors2, self.acceptors2 = self.extract_donor_acceptor(self.tome2)

        # sort inputs into temporary files
        self.sorted_input1 = None
        self.sorted_input2 = None
        with tempfile.NamedTemporaryFile(delete=False) as tmp_file1, tempfile.NamedTemporaryFile(delete=False) as tmp_file2:
            self.sorted_input1 = tmp_file1.name
            self.sorted_input2 = tmp_file2.name

        sort_by_n_column(self.input1, self.sorted_input1, 0)
        sort_by_n_column(self.input2, self.sorted_input2, 0)

    def cleanup(self):
        try:
            if self.sorted_input1 is not None:
                os.remove(self.sorted_input1)
            if self.sorted_input2 is not None:
                os.remove(self.sorted_input2)
        except Exception as e:
            print(f"Error cleaning up temporary files: {e}")

    def _run(self,outfname,prev_bps=None):
        outFP = open(outfname,"w+")
        outFP.write("read_name\t" +
                    "genome1_read_breakpoint\t" +
                    "genome2_read_breakpoint\t" +
                    "genome1_read_start\t" +
                    "genome1_read_end\t" +
                    "genome2_read_start\t" +
                    "genome2_read_end\t" +
                    "genome1_seqid\t" +
                    "genome2_seqid\t" +
                    "genome1_breakpoint\t" +
                    "genome2_breakpoint\t" +
                    "genome1_start\t" +
                    "genome1_end\t" +
                    "genome2_start\t" +
                    "genome2_end\t" +
                    "orientation\t" +
                    "reversed\t" +
                    "score\t" +
                    "junction1\t" +
                    "junction2\t" +
                    "gene1\t" +
                    "gene2\t" +
                    "binread\n")

        # iterate over read groups
        for read_name, lines1, lines2 in self.next_read_group():
            if len(lines1)==0 or len(lines2)==0:
                continue

            # create all unique combinations of two lists
            breakpoints = []
            for line1,line2 in [(x,y) for x in lines1 for y in lines2]:
                res = self.process(line1,line2,prev_bps)
                if res is not None:
                    it1 = res.read1.to_interval()
                    it2 = res.read2.to_interval()
                    found_genes1 = set()
                    found_genes2 = set()
                    for seqid,genes in self.gene_trees1.items():
                        inter = genes.overlap(it1)
                        if len(inter)>0 and seqid[0]==res.read1.sseqid:
                            found_genes1.update([x[2] for x in inter])
                    for seqid,genes in self.gene_trees2.items():
                        inter = genes.overlap(it2)
                        if len(inter)>0 and seqid[0]==res.read2.sseqid:
                            found_genes2.update([x[2] for x in inter])
                    
                    if len(found_genes1)==0:
                        found_genes1.add("-")
                    if len(found_genes2)==0:
                        found_genes2.add("-")
                    
                    res.genes = (found_genes1,found_genes2)
                    breakpoints.append(res)

            if len(breakpoints)==0:
                continue

            breakpoint = max(breakpoints, key=lambda x: x.score)
            outFP.write(str(breakpoint)+"\n")

        outFP.close()
        return

    def _run_pass2(self):
        # read in the output of the first pass and collect all breakpoints
        # extra weight is assigned if the breakpoint matches donor/acceptor or both
        pass1_bps = {} # (seqid,bp):weighted count
        with open(self.output,"r") as inFP:
            next(inFP) # skip header
            for line in inFP:
                lcs = line.strip().split("\t")
                seqid1 = lcs[7]
                seqid2 = lcs[8]
                bp1 = int(lcs[9])
                bp2 = int(lcs[10])
                sj1 = lcs[18]!="-"
                sj2 = lcs[19]!="-"
                weight1 = 1 + (1 if sj1 else 0) + (1 if sj2 else 0)
                weight2 = 1 + (1 if sj2 else 0) + (1 if sj1 else 0)
                pass1_bps.setdefault((seqid1,bp1),0)
                pass1_bps.setdefault((seqid2,bp2),0)
                pass1_bps[(seqid1,bp1)] += weight1
                pass1_bps[(seqid2,bp2)] += weight2

        # rerun ris detection, this time adding weights according to the previous round
        self._run(self.outfname_pass2,pass1_bps)

    def run(self):
        group_in_fname = None
        try:
            self._run(self.output)

            group_in_fname = self.output

            if self.two_pass:
                self._run_pass2()
                group_in_fname = self.output + ".corrected"

            groups = None
            if self.group:
                with open(self.output + ".grouped", "w+") as group_outFP:
                    group_outFP.write("genome1_seqid\t" +
                                    "genome2_seqid\t" +
                                    "genome1_breakpoint\t" +
                                    "genome2_breakpoint\t" +
                                    "orientation\t" +
                                    "genome1_start\t" +
                                    "genome1_end\t" +
                                    "genome2_start\t" +
                                    "genome2_end\t" +
                                    "count\t" +
                                    "junction1\t" +
                                    "junction2\t" +
                                    "gene1\t" +
                                    "gene2\n")
                    groups = self.group_breakpoints(group_in_fname, group_outFP)

            if self.chim_genome:
                with open(self.output + ".chimera.fasta", "w+") as chim_fastaFP, open(self.output + ".chimera.gtf", "w+") as chim_gtfFP:
                    self.generate_chimeras(chim_fastaFP, chim_gtfFP, groups)

        except Exception as e:
            print(e)
            sys.exit(1)

    def extract_region(self,fasta_path, seqid, start, end, reverse_complement=False):
        """Extract the sequence region with overhang, handling strand orientation."""
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id == seqid:
                seq = record.seq[start : end]
                if reverse_complement:
                    seq = seq.reverse_complement()
                return str(seq)
        raise ValueError(f"Sequence ID {seqid} not found in FASTA file {fasta_path}")

    def truncate_transcripts(self,interval_tree, seqid, strand, start, end, offset, reversed=False):
        """
        Truncate transcripts in the interval tree to the specified region
        and adjust their coordinates based on the offset in the chimeric genome.

        Parameters:
        - interval_tree: Interval tree for the genome.
        - seqid (str): Sequence ID of the region in the genome.
        - strand (str): Strand of the region (`+` or `-`).
        - start (int): Start position in the original genome.
        - end (int): End position in the original genome.
        - offset (int): Offset in the chimeric genome for this region.
        - reversed (bool): Whether the sequence has been reverse complemented.

        Returns:
        - list: List of truncated and offset GTF lines.
        """
        truncated = []
        if not (seqid, strand) in interval_tree:
            return truncated
        for interval in interval_tree[(seqid, strand)].overlap(start, end):
            gid, feature_type = interval.data
            original_start = max(start, interval.begin)
            original_end = min(end, interval.end)
            
            # Calculate new coordinates
            new_strand = strand
            if reversed:
                new_start = offset + (end - original_end) + 1
                new_end = offset + (end - original_start) + 1
                # Flip strand if reversed
                new_strand = "+" if strand == "-" else "-"
            else:
                new_start = original_start - start + offset + 1
                new_end = original_end - start + offset + 1
            
            # Ensure correct order of start and end
            if new_start > new_end:
                new_start, new_end = new_end, new_start

            truncated.append(
                f"chimeric\t{feature_type}\t{new_start}\t{new_end}\t.\t{new_strand}\t.\ttranscript_id \"{gid}\""
            )
        return truncated

    def create_chimeric_genome(self,
        seqid1, start1, end1, seqid2, start2, end2, orientation
    ):
        """
        Create a chimeric genome from two genome FASTA files and merge GTF annotations.

        Parameters:
        - seqid1 (str): Sequence ID of the region in the first genome.
        - start1 (int): Start position of the first genome region (1-based inclusive).
        - end1 (int): End position of the first genome region (1-based inclusive).
        - seqid2 (str): Sequence ID of the region in the second genome.
        - start2 (int): Start position of the second genome region (1-based inclusive).
        - end2 (int): End position of the second genome region (1-based inclusive).
        - orientation (int): 1 if genome1 is first and genome2 is second; 2 otherwise.

        Returns:
        - tuple: (chimeric sequence string, merged GTF records).
        """

        # Adjust coordinates for overhang
        adj_start1 = max(0, start1 - self.overhang) if start1 < end1 else max(0, end1 - self.overhang)
        adj_end1 = end1 + self.overhang if start1 < end1 else start1 + self.overhang
        adj_start2 = max(0, start2 - self.overhang) if start2 < end2 else max(0, end2 - self.overhang)
        adj_end2 = end2 + self.overhang if start2 < end2 else start2 + self.overhang
        strand1 = "+" if start1 < end1 else "-"
        strand2 = "+" if start2 < end2 else "-"

        # Extract sequences with overhang
        seq1 = self.extract_region(self.genome1, seqid1, adj_start1, adj_end1, strand1 == "-")
        seq2 = self.extract_region(self.genome2, seqid2, adj_start2, adj_end2, strand2 == "-")

        # Combine sequences based on orientation
        if orientation == 1:
            chimeric_seq = seq1 + seq2
            offset1 = 0
            offset2 = len(seq1)
        elif orientation == 2:
            chimeric_seq = seq2 + seq1
            offset1 = len(seq2)
            offset2 = 0
        else:
            raise ValueError("Invalid orientation: must be 1 or 2.")

        # Create transcript records for full genome regions
        merged_gtf = []
        merged_gtf.append(f"chimeric\ttranscript\t{offset1 + 1}\t{offset1 + len(seq1)}\t.\t{strand1}\t.\ttranscript_id \"genome1\"")
        merged_gtf.append(f"chimeric\texon\t{offset1 + 1}\t{offset1 + len(seq1)}\t.\t{strand1}\t.\ttranscript_id \"genome1\"")
        merged_gtf.append(f"chimeric\ttranscript\t{offset2 + 1}\t{offset2 + len(seq2)}\t.\t{strand2}\t.\ttranscript_id \"genome2\"")
        merged_gtf.append(f"chimeric\texon\t{offset2 + 1}\t{offset2 + len(seq2)}\t.\t{strand2}\t.\ttranscript_id \"genome2\"")

        # create transcript records defined by the breakpoint without the overhang
        merged_gtf.append(f"chimeric\ttranscript\t{offset1 + 1 + self.overhang}\t{offset1 + len(seq1) - self.overhang}\t.\t{strand1}\t.\ttranscript_id \"genome1_breakpoint\"")
        merged_gtf.append(f"chimeric\texon\t{offset1 + 1 + self.overhang}\t{offset1 + len(seq1) - self.overhang}\t.\t{strand1}\t.\ttranscript_id \"genome1_breakpoint\"")
        merged_gtf.append(f"chimeric\ttranscript\t{offset2 + 1 + self.overhang}\t{offset2 + len(seq2) - self.overhang}\t.\t{strand2}\t.\ttranscript_id \"genome2_breakpoint\"")
        merged_gtf.append(f"chimeric\texon\t{offset2 + 1 + self.overhang}\t{offset2 + len(seq2) - self.overhang}\t.\t{strand2}\t.\ttranscript_id \"genome2_breakpoint\"")

        # Add truncated transcripts from both genomes
        merged_gtf.extend(self.truncate_transcripts(self.gene_trees1, seqid1, strand1, adj_start1, adj_end1, offset1, strand1 == "-"))
        merged_gtf.extend(self.truncate_transcripts(self.gene_trees2, seqid2, strand2, adj_start2, adj_end2, offset2, strand2 == "-"))

        return chimeric_seq, merged_gtf

    def generate_chimeras(self,chim_fastaFP,chim_gtfFP,groups):
        for k,v in groups.items():
            chim_seq, merged_gtf = self.create_chimeric_genome(k[0],v["g1_start"],v["g1_end"],k[1],v["g2_start"],v["g2_end"],k[4])
            chim_seqid = f"{k[0]}_{v['g1_start']}_{v['g1_end']}_{k[1]}_{v['g2_start']}_{v['g2_end']}_{k[4]}"
            chim_fastaFP.write(f">{chim_seqid}\n{chim_seq}\n")
            for line in merged_gtf:
                chim_gtfFP.write(chim_seqid+"\t"+line+"\n")

    def group_breakpoints(self,in_fname,outFP):
        # read in the input, count the number of times each breakpoint is observed, output the result

        groups = {} # (seqid1,seqid2,bp1,bp2,orientation,sj1,sj2,gene1,gene2):{"count":count,"g1_start":[],"g1_end":[],"g2_start":[],"g2_end":[]}

        with open(in_fname,"r") as inFP:
            next(inFP)
            for line in inFP:
                lcs = line.strip().split("\t")
                seqid1 = lcs[7]
                seqid2 = lcs[8]
                bp1 = int(lcs[9])
                bp2 = int(lcs[10])
                g1_start = int(lcs[11])
                g1_end = int(lcs[12])
                g2_start = int(lcs[13])
                g2_end = int(lcs[14])
                orientation = int(lcs[15])
                sj1 = lcs[18]
                sj2 = lcs[19]
                gene1 = lcs[20]
                gene2 = lcs[21]

                k = (seqid1,seqid2,bp1,bp2,orientation,sj1,sj2,gene1,gene2)
                groups.setdefault(k,{"count":0,
                                    "g1_start":[g1_start],
                                    "g1_end":[g1_end],
                                    "g2_start":[g2_start],
                                    "g2_end":[g2_end]})
                                                                                        
                groups[k]["count"] +=1
                groups[k]["g1_start"].append(g1_start)
                groups[k]["g1_end"].append(g1_end)
                groups[k]["g2_start"].append(g2_start)
                groups[k]["g2_end"].append(g2_end)


        for k,v in groups.items():
            g1_start = min(v["g1_start"]) if v["g1_start"] < v["g1_end"] else max(v["g1_start"])
            g1_end = max(v["g1_end"]) if v["g1_start"] < v["g1_end"] else min(v["g1_end"])
            g2_start = min(v["g2_start"]) if v["g2_start"] < v["g2_end"] else max(v["g2_start"])
            g2_end = max(v["g2_end"]) if v["g2_start"] < v["g2_end"] else min(v["g2_end"])

            # overwrite the start and end positions
            v["g1_start"] = g1_start
            v["g1_end"] = g1_end
            v["g2_start"] = g2_start
            v["g2_end"] = g2_end

            outFP.write(k[0]+"\t"+
                        k[1]+"\t"+
                        str(k[2])+"\t"+
                        str(k[3])+"\t"+
                        str(k[4])+"\t"+
                        str(g1_start)+"\t"+
                        str(g1_end)+"\t"+
                        str(g2_start)+"\t"+
                        str(g2_end)+"\t"+
                        str(v["count"])+"\t"+
                        k[5]+"\t"+
                        k[6]+"\t"+
                        k[7]+"\t"+
                        k[8]+"\n")

        return groups

    def match_donor_acceptor(self,donors,acceptors):
        # finds pairs of donors, acceptors that are adjacent to each other
        res = []
        if len(donors)==0:
            for y in acceptors:
                res.append([[None,None,None],y])
                res[-1][1].append(self.half_weight)
        if len(acceptors)==0:
            for x in donors:
                res.append([x,[None,None,None]])
                res[-1][0].append(self.half_weight)
        for x in donors:
            for y in acceptors:
                if y[0] - x[0] == 1:
                    res.append([x,y])
                    res[-1][0].append(self.full_weight)
                    res[-1][1].append(self.full_weight)
                else:
                    res.append([x,[None,None,None]])
                    res.append([[None,None,None],y])
                    res[-2][0].append(self.half_weight)
                    res[-1][1].append(self.half_weight)
        return res

    def merge_weights(self, l1, l2):
        res = []
        used_y_idxs = set()
        
        for x in l1:
            r = copy.deepcopy(x)
            for y in l2:
                if y[0][0] == x[0][0]:
                    r[1] += y[1]
                    used_y_idxs.add(l2.index(y))
                if y[1] == x[0]:
                    r[1] += y[0]
                    used_y_idxs.add(l2.index(y))
                if y[0] == x[1]:
                    r[1] += y[1]
                    used_y_idxs.add(l2.index(y))
        
        return res

    def _process(self,binread:Binread,forward,orientation):
        copy_binread = copy.deepcopy(binread) # create a master copy of the binread
        copy_binread.orientation = orientation
        if not forward:
            copy_binread.reversed = True
            copy_binread.reverse()

        # add donor/acceptor sites
        copy_binread.read1.load_donors(self.donors1)
        copy_binread.read1.load_acceptors(self.acceptors1)
        copy_binread.read2.load_donors(self.donors2)
        copy_binread.read2.load_acceptors(self.acceptors2)

        da = None
        if orientation==1:
            da = self.match_donor_acceptor(copy_binread.read1.donors,copy_binread.read2.acceptors)
        else:
            da = self.match_donor_acceptor(copy_binread.read2.donors,copy_binread.read1.acceptors)
        weight_pairs = [[[None,None,None],[None,None,None]]] # weights are stored as: [[read1_pos,read1_label,read1_weight],[read2_pos,read2_label,read2_weight]]
        for pair in da:
            weight_pairs.append([pair[0],pair[1]])
        
        results = []
        for weight_pair in weight_pairs:
            base_binread = copy.deepcopy(copy_binread) # create a master copy of the binread
            base_binread.find_breakpoint(weight_pair,orientation)
            results.append(base_binread)

        if len(results)==0:
            return None
        else:
            return max(results, key=lambda x: x.score)

    def process(self,m1,m2,pass1_bps=None):
        # take two mappings of the same read and find the breakpoint between them

        binread = Binread()
        binread.add_read1(m1)
        binread.add_read2(m2)

        # orientation 1 means that in the RIS read1 is on the left, read2 is on the right
        # orientation 2 means the opposite
        f1 = self._process(binread,True,1)
        f2 = self._process(binread,True,2)
        r1 = self._process(binread,False,1)
        r2 = self._process(binread,False,2)

        return max([f1,f2,r1,r2], key=lambda x: x.score)

    def extract_genes(self, tome):
        gene_trees = {}

        for tx in tome.transcript_it():
            gid = tx.get_gid()
            seqid = tx.get_seqid()
            strand = tx.get_strand()
            start = tx.get_start()
            end = tx.get_end()
            gene_trees.setdefault((seqid,strand),IntervalTree()).addi(start,end,(gid,"transcript"))
            for exon in tx.get_exons():
                start = exon[2].get_start()
                end = exon[2].get_end()
                gene_trees.setdefault((seqid,strand),IntervalTree()).addi(start,end,(gid,"exon"))

        return gene_trees
    
    def extract_donor_acceptor(self, tome):
        donors = {}
        acceptors = {}

        for tx in tome.transcript_it():
            gid = tx.get_gid()
            seqid = tx.get_seqid()
            strand = tx.get_strand()
            
            for it in tx.introns_it():
                donor = it[0]-1 if strand=="+" else it[1]
                acceptor = it[1] if strand=="+" else it[0]-1

                donors.setdefault(seqid,{})[donor] = (gid,strand)
                acceptors.setdefault(seqid,{})[acceptor] = (gid,strand)
        
        return donors, acceptors
    
    def next_read_group(self):
        # iterate over lines of two sorted files
        # group lines by read name (1st column)
        # yield two lists of lines, one from each file, that have the same read name
        with open(self.input1, 'r') as inFP1, open(self.input2, 'r') as inFP2:
            iter1, iter2 = iter(inFP1), iter(inFP2)
            line1, line2 = next(iter1, None), next(iter2, None)
            while line1 == "\n":
                line1 = next(iter1, None)
            while line2 == "\n":
                line2 = next(iter2, None)

            while line1 is not None or line2 is not None:
                current_read_name = None
                lines1, lines2 = [], []

                if line1 is not None and (line2 is None or line1.split("\t")[0] <= line2.split("\t")[0]):
                    current_read_name = line1.split("\t")[0]
                    while line1 is not None and line1.split("\t")[0] == current_read_name:
                        lines1.append(line1.strip())
                        line1 = next(iter1, None)

                if line2 is not None and (line1 is None or line2.split("\t")[0] <= line1.split("\t")[0]):
                    current_read_name = line2.split("\t")[0] if current_read_name is None else current_read_name
                    while line2 is not None and line2.split("\t")[0] == current_read_name:
                        lines2.append(line2.strip())
                        line2 = next(iter2, None)

                yield current_read_name, lines1, lines2
    
def main():
    parser = argparse.ArgumentParser(description="Detect Precise Chimeric Breakpoints.")

    parser.add_argument('-i1',
                        '--input1',
                        required=True,
                        type=str,
                        help="File containing mapping of reads to genome #1.")
    parser.add_argument('-i2',
                        '--input2',
                        required=True,
                        type=str,
                        help="File containing mapping of reads to genome #2.")
    parser.add_argument('-a1',
                        '--annotation1',
                        required=True,
                        type=str,
                        help="GTF file containing gene annotations for genome #1.")
    parser.add_argument('-a2',
                        '--annotation2',
                        required=True,
                        type=str,
                        help="GTF file containing gene annotations for genome #2.")
    parser.add_argument('--two-pass',
                        required=False,
                        action='store_true',
                        help="Run two pass approach. First pass will find all possible breakpoints. Second pass will try to match breakpoints that are close to each other.")
    parser.add_argument('-g',
                        '--group',
                        required=False,
                        action='store_true',
                        help="If enabled, will output a file with breakpoints groupped by position.")
    parser.add_argument('--chim-genome',
                        required=False,
                        action='store_true',
                        help="(Requires -group). If enabled, will generate a fasta file with chimeric genome sequences, stitching together the two genomes at the breakpoints.")
    parser.add_argument('-g1',
                        "--genome1",
                        required=False,
                        type=str,
                        help="Path to the first genome.")
    parser.add_argument('-g2',
                        "--genome2",
                        required=False,
                        type=str,
                        help="Path to the second genome.")
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help="Output file.")
    parser.add_argument('--max-dist',
                        required=False,
                        type=int,
                        default=5,
                        help="Maximum distance between breakpoints of the two segments. Default: 5.")
    parser.add_argument('--max-weight',
                        required=False,
                        type=int,
                        default=5,
                        help="Maximum weight of a breakpoint when biasing the 2nd pass. Default: 5.")
    parser.add_argument('--full-weight',
                        required=False,
                        type=int,
                        default=5,
                        help="Weight of a breakpoint that matches donor and acceptor. Default: 5.")
    parser.add_argument('--half-weight',
                        required=False,
                        type=int,
                        default=3,
                        help="Weight of a breakpoint that matches either donor or acceptor. Default: 3.")
    parser.add_argument('--overhang',
                        required=False,
                        type=int,
                        default=1000,
                        help="Number of bases to include in the chimeric genome overhang. Default: 1000.")

    args = parser.parse_args()

    iris = IRIS(args)
    iris.run()

if __name__ == "__main__":
    main()