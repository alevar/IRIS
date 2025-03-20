import argparse

from .utils.common import *

from .classes.txgroup import Transcriptome
from .classes.binread import Binread
from .classes.read import Read

class IRIS:
    def __init__(self, args):         
        # INPUT FILES
        self.input1 = args.i1
        self.input2 = args.i2
        self.annotation1 = args.annotation1
        self.annotation2 = args.annotation2
        self.output = args.output
        self.two_pass = args.two_pass
        self.group = args.group
        self.max_dist = args.max_dist
        self.max_weight = args.max_weight
        self.full_weight = args.full_weight
        self.half_weight = args.half_weight
        
        self.tome1 = None if self.annotation1 is None else Transcriptome()
        if self.annotation1 is not None:
            self.tome1.build_from_file(self.annotation1)
        self.tome2 = None if self.annotation2 is None else Transcriptome()
        if self.annotation2 is not None:
            self.tome2.build_from_file(self.annotation2)

    def run(self):
        return
        
    def match_donor_acceptor(
        donors: List[List[Any]], 
        acceptors: List[List[Any]], 
        orientation: int,
        full_weight: int, 
        half_weight: int
    ) -> List[List[List[Any]]]:
        """
        Find pairs of donors and acceptors that are adjacent to each other.
        
        Args:
            donors (List[List[Any]]): List of donor sites
            acceptors (List[List[Any]]): List of acceptor sites
            orientation (int): Orientation of the alignment
            full_weight (int): Weight for full matches
            half_weight (int): Weight for partial matches
            
        Returns:
            List[List[List[Any]]]: Matched donor-acceptor pairs
        """
        res = []
        if len(donors) == 0:
            for y in acceptors:
                res.append([[None, None, None], y])
                res[-1][1].append(half_weight)
        if len(acceptors) == 0:
            for x in donors:
                res.append([x, [None, None, None]])
                res[-1][0].append(half_weight)
        for x in donors:
            for y in acceptors:
                if y[0] - x[0] == 1:
                    res.append([x, y])
                    res[-1][0].append(full_weight)
                    res[-1][1].append(full_weight)
                else:
                    res.append([x, [None, None, None]])
                    res.append([[None, None, None], y])
                    res[-2][0].append(half_weight)
                    res[-1][1].append(half_weight)
        return res


    def _process(
        binread: Binread,
        forward: bool,
        orientation: int,
        donors1: Dict[str, Dict[int, Any]],
        acceptors1: Dict[str, Dict[int, Any]],
        donors2: Dict[str, Dict[int, Any]],
        acceptors2: Dict[str, Dict[int, Any]],
        args: Any
    ) -> Optional[Binread]:
        """
        Process a binread object to find integration sites.
        
        Args:
            binread (Binread): Binread object to process
            forward (bool): Whether to process in forward direction
            orientation (int): Orientation of the alignment
            donors1 (Dict[str, Dict[int, Any]]): Donors for the first genome
            acceptors1 (Dict[str, Dict[int, Any]]): Acceptors for the first genome
            donors2 (Dict[str, Dict[int, Any]]): Donors for the second genome
            acceptors2 (Dict[str, Dict[int, Any]]): Acceptors for the second genome
            args (Any): Arguments with configuration options
            
        Returns:
            Optional[Binread]: Processed Binread object, or None if no results found
        """
        copy_binread = copy.deepcopy(binread)  # create a master copy of the binread
        if not forward:
            copy_binread.reverse()

        # add donor/acceptor sites
        copy_binread.read1.load_donors(donors1)
        copy_binread.read1.load_acceptors(acceptors1)
        copy_binread.read2.load_donors(donors2)
        copy_binread.read2.load_acceptors(acceptors2)

        da = None
        if orientation == 1:
            da = match_donor_acceptor(
                copy_binread.read1.donors,
                copy_binread.read2.acceptors,
                orientation,
                args.full_weight,
                args.half_weight
            )
        else:
            da = match_donor_acceptor(
                copy_binread.read2.donors,
                copy_binread.read1.acceptors,
                orientation,
                args.full_weight,
                args.half_weight
            )
        
        weight_pairs = [[[None, None, None], [None, None, None]]]  # weights are stored as: [[read1_pos,read1_label,read1_weight],[read2_pos,read2_label,read2_weight]]
        for pair in da:
            weight_pairs.append([pair[0], pair[1]])
        
        results = []
        for weight_pair in weight_pairs:
            base_binread = copy.deepcopy(copy_binread)  # create a copy of the binread
            base_binread.orientation = orientation
            base_binread.find_breakpoint(weight_pair, orientation)
            results.append(base_binread)

        if len(results) == 0:
            return None
        else:
            return max(results, key=lambda x: x.score)


    def process(
        m1: str,
        m2: str,
        donors1: Dict[str, Dict[int, Any]],
        acceptors1: Dict[str, Dict[int, Any]],
        donors2: Dict[str, Dict[int, Any]],
        acceptors2: Dict[str, Dict[int, Any]],
        args: Any,
        pass1_bps: Optional[Dict[Tuple[str, int], int]] = None
    ) -> Binread:
        """
        Process two mappings of the same read to find integration sites.
        
        Args:
            m1 (str): First mapping line
            m2 (str): Second mapping line
            donors1 (Dict[str, Dict[int, Any]]): Donors for the first genome
            acceptors1 (Dict[str, Dict[int, Any]]): Acceptors for the first genome
            donors2 (Dict[str, Dict[int, Any]]): Donors for the second genome
            acceptors2 (Dict[str, Dict[int, Any]]): Acceptors for the second genome
            args (Any): Arguments with configuration options
            pass1_bps (Optional[Dict[Tuple[str, int], int]]): First pass breakpoints
            
        Returns:
            Binread: Processed Binread object with the highest score
        """
        binread = Binread()
        binread.add_read1(m1)
        binread.add_read2(m2)

        f1 = _process(binread, True, 1, donors1, acceptors1, donors2, acceptors2, args)
        f2 = _process(binread, True, 2, donors1, acceptors1, donors2, acceptors2, args)
        r1 = _process(binread, False, 1, donors1, acceptors1, donors2, acceptors2, args)
        r2 = _process(binread, False, 2, donors1, acceptors1, donors2, acceptors2, args)

        results = [r for r in [f1, f2, r1, r2] if r is not None]
        if not results:
            return None
            
        return max(results, key=lambda x: x.score)


    def next_read_group(fname1: str, fname2: str):
        """
        Iterate over lines of two sorted files, grouping by read name.
        
        Args:
            fname1 (str): First file path
            fname2 (str): Second file path
            
        Yields:
            Tuple[str, List[str], List[str]]: Read name, lines from file 1, lines from file 2
        """
        with open(fname1, 'r') as inFP1, open(fname2, 'r') as inFP2:
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
    parser.add_argument('--two_pass',
                        required=False,
                        action='store_true',
                        help="Run two pass approach. First pass will find all possible breakpoints. Second pass will try to match breakpoints that are close to each other.")
    parser.add_argument('--group',
                        required=False,
                        action='store_true',
                        help="If enabled, will output a file with breakpoints groupped by position.")
    parser.add_argument('-o',
                        '--output',
                        required=True,
                        type=str,
                        help="Output file.")
    parser.add_argument('-max_dist',
                        required=False,
                        type=int,
                        default=5,
                        help="Maximum distance between breakpoints of the two segments. Default: 5.")
    parser.add_argument('-max_weight',
                        required=False,
                        type=int,
                        default=5,
                        help="Maximum weight of a breakpoint when biasing the 2nd pass. Default: 5.")
    parser.add_argument('-full_weight',
                        required=False,
                        type=int,
                        default=5,
                        help="Weight of a breakpoint that matches donor and acceptor. Default: 5.")
    parser.add_argument('-half_weight',
                        required=False,
                        type=int,
                        default=3,
                        help="Weight of a breakpoint that matches either donor or acceptor. Default: 3.")

    args = parser.parse_args()

    iris = IRIS(args)
    iris.run()

if __name__ == "__main__":
    main()