"""
Module containing the Read class for representing and analyzing BLAST alignments.
"""

import re
from typing import List, Dict, Tuple, Optional, Any, Union


class Read:
    """
    Represents a BLAST alignment for a single read.
    
    Attributes:
        qseqid (str): Query sequence ID
        qlen (int): Query sequence length
        sseqid (str): Subject sequence ID
        qstart (int): Start position in query sequence
        qend (int): End position in query sequence
        sstart (int): Start position in subject sequence
        send (int): End position in subject sequence
        btop (str): BLAST traceback operations
        btopl (list): Parsed BTOP string
        binread (list): Binary representation of the alignment
        donors (list): Donor sites
        acceptors (list): Acceptor sites
        weights (list): Weights associated with positions
    """

    def __init__(self):
        """Initialize a new Read object with default values."""
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

    def to_interval(self) -> Tuple[int, int, 'Read']:
        """
        Convert the read alignment to an interval.
        
        Returns:
            Tuple[int, int, Read]: A tuple containing the start, end positions and the Read object
        """
        return (min(self.sstart, self.send), max(self.sstart, self.send), self)

    def from_line(self, line: str) -> None:
        """
        Parse a BLAST output line to populate the Read object.
        
        Args:
            line (str): A line from the BLAST output
        """
        self.qseqid, self.qlen, self.sseqid, self.qstart, self.qend, self.sstart, self.send, self.btop = line.strip().split("\t")
        self.qlen = int(self.qlen)
        self.qstart = int(self.qstart)
        self.qend = int(self.qend)
        self.sstart = int(self.sstart)
        self.send = int(self.send)

        self.parse_btop()
        self.btop_to_list()

    def is_reversed(self) -> bool:
        """
        Check if the alignment is reversed.
        
        Returns:
            bool: True if the alignment is reversed, False otherwise
        """
        return self.sstart > self.send
    
    def reverse(self) -> None:
        """Reverse the alignment coordinates and associated data."""
        self.sstart, self.send = self.send, self.sstart
        self.qstart, self.qend = self.qlen - self.qend + 1, self.qlen - self.qstart + 1
        self.btopl = self.btopl[::-1]
        self.btop_to_list()
        self.donors = [[self.qlen - x[0], x[1], x[2]] for x in self.donors]
        self.acceptors = [[self.qlen - x[0], x[1], x[2]] for x in self.acceptors]
        self.weights = [[self.qlen - x[0], x[1], x[2]] for x in self.weights]

    def parse_btop(self) -> None:
        """Parse the BTOP string into a list of operations."""
        self.btopl = []
        pattern = re.compile(r'(\d+)|([A-Za-z-]{2})')
        matches = re.finditer(pattern, self.btop)

        for match in matches:
            if match.group(1):
                self.btopl.append(int(match.group(1)))
            else:
                self.btopl.append(match.group(2))

    def btop_to_list(self) -> None:
        """Convert BTOP operations to a binary representation of the alignment."""
        self.binread = [0] * self.qlen
        index = self.qstart - 1

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(index, index + b, 1):
                    self.binread[i] += 1
                index += b
            elif isinstance(b, str):
                if b[0] == "-":  # insertion
                    # if insertion - decrement the score of the next element
                    # this is equivalent to saying, by that position the score on the reference would have decreased by this much due to insertion penalty
                    self.binread[index] -= 1
                elif b[1] == "-":  # deletion
                    self.binread[index] -= 1
                    index += 1
                else:  # mismatch
                    self.binread[index] = 0
                    index += 1

    def get_sites(self, sites: Dict[int, Any]) -> List[List[Any]]:
        """
        Get positions of specific sites in the read.
        
        Args:
            sites (Dict[int, Any]): Dictionary mapping genome positions to annotations
            
        Returns:
            List[List[Any]]: List of sites found in the read
        """
        res = []

        # handle the case when the site is upstream of the read start position
        for i in range(self.qstart - 2, -1, -1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = self.sstart - (self.qstart - 1 - i)
            else:
                genome_pos = self.sstart + (self.qstart - 1 - i)
            if genome_pos in sites:
                res.append([i, sites[genome_pos]])
        
        index_read = self.qstart - 1
        index_genome = self.sstart
        inc = 1 if self.sstart < self.send else -1

        if index_genome in sites:
            res.append([index_read, sites[index_genome]])

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0, b, 1):
                    index_genome += inc
                    index_read += 1
                    if index_genome in sites:
                        res.append([index_read, sites[index_genome]])
            elif isinstance(b, str):
                if b[0] == "-":  # insertion in read
                    index_genome += inc
                elif b[1] == "-":  # deletion from read
                    index_read += 1
                else:  # mismatch - treat as a regular match here
                    index_genome += inc
                    index_read += 1
                    if index_genome in sites:
                        res.append([index_read, sites[index_genome]])

        # handle the case when the site is downstream of the read end position
        for i in range(index_read, self.qlen, 1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = index_genome + (i - index_read)
            else:
                genome_pos = index_genome - (i - index_read)
            if genome_pos in sites:
                res.append([i, sites[genome_pos]])

        return res
    
    def load_donors(self, donors: Dict[str, Dict[int, Any]]) -> None:
        """
        Load donor sites from a dictionary.
        
        Args:
            donors (Dict[str, Dict[int, Any]]): Dictionary mapping sequence IDs to donors
        """
        if self.sseqid in donors:
            self.donors = self.get_sites(donors[self.sseqid])

    def load_acceptors(self, acceptors: Dict[str, Dict[int, Any]]) -> None:
        """
        Load acceptor sites from a dictionary.
        
        Args:
            acceptors (Dict[str, Dict[int, Any]]): Dictionary mapping sequence IDs to acceptors
        """
        if self.sseqid in acceptors:
            self.acceptors = self.get_sites(acceptors[self.sseqid])

    def read2genome(self, pos: int) -> Optional[int]:
        """
        Convert a position in the read to a position in the genome.
        
        Args:
            pos (int): Position in the read
            
        Returns:
            Optional[int]: Corresponding position in the genome, or None if not found
        """
        for i in range(self.qstart - 2, -1, -1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = self.sstart - (self.qstart - 1 - i)
            else:
                genome_pos = self.sstart + (self.qstart - 1 - i)
            if i == pos:
                return genome_pos
        
        index_read = self.qstart - 1
        index_genome = self.sstart
        inc = 1 if self.sstart < self.send else -1

        if index_read == pos:
            return index_genome

        for b in self.btopl:
            if isinstance(b, int):
                for i in range(0, b, 1):
                    index_genome += inc
                    index_read += 1
                    if index_read == pos:
                        return index_genome
            elif isinstance(b, str):
                if b[0] == "-":  # insertion in read
                    index_genome += inc
                elif b[1] == "-":  # deletion from read
                    index_read += 1
                else:  # mismatch - treat as a regular match here
                    index_genome += inc
                    index_read += 1
                    if index_read == pos:
                        return index_genome

        # handle the case when the site is downstream of the read end position
        for i in range(index_read, self.qlen, 1):
            genome_pos = None
            if self.sstart < self.send:
                genome_pos = index_genome + (i - index_read)
            else:
                genome_pos = index_genome - (i - index_read)
            if i == pos:
                return genome_pos
                
        return None