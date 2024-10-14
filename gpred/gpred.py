import argparse
import sys
import os
import csv
import re
from textwrap import fill
from re import Pattern
from pathlib import Path
from typing import List, Union, Optional


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='genome_file', type=isfile, required=True,
                        help="Complete genome file in fasta format")
    parser.add_argument('-g', dest='min_gene_len', type=int,
                        default=50, help="Minimum gene length to consider (default 50).")
    parser.add_argument('-s', dest='max_shine_dalgarno_distance', type=int,
                        default=16, help="Maximum distance from start codon "
                        "where to look for a Shine-Dalgarno motif (default 16).")
    parser.add_argument('-d', dest='min_gap', type=int, default=40,
                        help="Minimum gap between two genes - shine box not included (default 40).")
    parser.add_argument('-p', dest='predicted_genes_file', type=Path,
                        default=Path("predict_genes.csv"),
                        help="Tabular file giving position of predicted genes")
    parser.add_argument('-o', dest='fasta_file', type=Path,
                        default=Path("genes.fna"),
                        help="Fasta file giving sequence of predicted genes")
    return parser.parse_args()


def read_fasta(fasta_file: Path) -> str:
    """Extract genome sequence from fasta files.

    :param fasta_file: (Path) Path to the fasta file.
    :return: (str) Sequence from the genome. 
    """
    sequence = ""
    with open(fasta_file, "r") as f:
        for line in f:
            if not line.startswith(">"):
                sequence += line.strip().upper()

    return sequence


def find_start(start_regex: Pattern, sequence: str, start: int, stop: int) -> Union[int, None]:
    """Find next start codon before a end position.

    :param start_regexp: A regex object that identifies a start codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :param stop: (int) Stop position of the research
    :return: (int) If exist, position of the start codon. Otherwise None. 
    """
    starts = start_regex.search(sequence, start, stop)
    if not starts:
        return None
    else:
        return starts.start(0)


def find_stop(stop_regex: Pattern, sequence: str, start: int) -> Union[int, None]:
    """Find next stop codon that should be in the same reading phase as the start.

    :param stop_regexp: A regex object that identifies a stop codon.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Start position of the research
    :return: (int) If exist, position of the stop codon. Otherwise None. 
    """
    stops = stop_regex.finditer(sequence, start, len(sequence))

    # If iterator is empty
    if stops is None:
        return None
    # Check if the stop codon is in frame with the start codon
    is_in_frame = False
    for stop in stops:
        if (stop.start(0)-start) % 3 == 0:
            is_in_frame = True
            return stop.start(0)
    if not is_in_frame:
        return None


def has_shine_dalgarno(
    shine_regex: Pattern,
    sequence: str, start: int,
    max_shine_dalgarno_distance: int) -> bool:
    """Find a shine dalgarno motif before the start codon

    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param sequence: (str) Sequence from the genome
    :param start: (int) Position of the start in the genome
    :param max_shine_dalgarno_distance: (int) Max dist of the shine dalgarno to the start position
    :return: (boolean) true -> has a shine dalgarno upstream to the gene, false -> no
    """
    search_start = start - max_shine_dalgarno_distance

    if search_start < 0:
        return False

    shine_dalgarno = shine_regex.search(sequence, search_start, start - 6)
    if shine_dalgarno is not None:
        return True
    else:
        return False


def predict_genes(
    sequence: str,
    start_regex: Pattern,
    stop_regex: Pattern,
    shine_regex: Pattern,
    min_gene_len: int,
    max_shine_dalgarno_distance: int,
    min_gap: int) -> List:
    """Predict most probable genes

    :param sequence: (str) Sequence from the genome.
    :param start_regexp: A regex object that identifies a start codon.
    :param stop_regexp: A regex object that identifies a stop codon.
    :param shine_regexp: A regex object that identifies a shine-dalgarno motif.
    :param min_gene_len: (int) Minimum gene length.
    :param max_shine_dalgarno_distance: (int) Max dist of the shine dalgarno to the start position.
    :param min_gap: (int) Minimum distance between two genes.
    :return: (list) List of [start, stop] position of each predicted genes.
    """
    current_pos = 0
    gene_list = []

    while (len(sequence) - current_pos) >= min_gap:
        current_pos = find_start(start_regex,
                                 sequence,
                                 current_pos,
                                 len(sequence))
        if current_pos is not None:
            stop = find_stop(stop_regex,
                             sequence,
                             current_pos)
            if stop is not None:
                gene_len = stop - current_pos + 1
                if gene_len >= min_gene_len:
                    if has_shine_dalgarno(shine_regex,
                                          sequence,
                                          current_pos,
                                          max_shine_dalgarno_distance):
                        # Probable gene identified, save it:
                        gene_list.append([current_pos+1, stop+3])
                        current_pos = stop + 3 + min_gap
                    else:
                        current_pos = current_pos + 1
                else:
                    current_pos = current_pos + 1
            else:
                current_pos = current_pos + 1
    return gene_list


def write_genes_pos(predicted_genes_file: Path, probable_genes: List[List[int]]) -> None:
    """Write list of gene positions.

    :param predicted_genes_file: (Path) Output file of gene positions.
    :param probable_genes: List of [start, stop] position of each predicted genes.
    """
    try:
        with predicted_genes_file.open("wt") as predict_genes:
            predict_genes_writer = csv.writer(predict_genes, delimiter=",")
            predict_genes_writer.writerow(["Start", "Stop"])
            predict_genes_writer.writerows(probable_genes)
    except IOError:
        sys.exit("Error cannot open {}".format(predicted_genes_file))


def write_genes(
    fasta_file: Path,
    sequence: str,
    probable_genes: List[List[int]],
    sequence_rc: str,
    probable_genes_comp: List[List[int]]):
    """Write gene sequence in fasta format

    :param fasta_file: (Path) Output fasta file.
    :param sequence: (str) Sequence of genome file in 5'->3'.
    :param probable_genes: (list) List of [start, stop] position of each predicted genes in 5'->3'.
    :param sequence_rc: (str) Sequence of genome file in 3' -> 5'.
    :param probable_genes_comp: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    """
    try:
        with open(fasta_file, "wt") as fasta:
            for i,gene_pos in enumerate(probable_genes):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                    i+1, os.linesep,
                    fill(sequence[gene_pos[0]-1:gene_pos[1]])))
            i = i+1
            for j,gene_pos in enumerate(probable_genes_comp):
                fasta.write(">gene_{0}{1}{2}{1}".format(
                            i+1+j, os.linesep,
                            fill(sequence_rc[gene_pos[0]-1:gene_pos[1]])))
    except IOError:
        sys.exit("Error cannot open {}".format(fasta_file))


def reverse_complement(sequence: str) -> str:
    """Get the reverse complement

    :param sequence: (str) DNA Sequence.
    :return: (str) Reverse complemented sequence.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in sequence[::-1]])


def reorder_gene_list(gene_list : list, sequence) -> list:
    """Reorder the gene list from a 3' -> 5' to 5' -> 3'

    :param gene_list: (list)List of [start, stop] position of each predicted genes in 3' -> 5'.
    :param sequence: (str) Sequence from the genome.
    :return: (list) List of [start, stop] position of each predicted genes in the 5' -> 3'.
    """
    new_gene_list = []
    for gene in gene_list:
        start = len(sequence) - gene[1] + 1
        stop = len(sequence) - gene[0] + 1
        new_gene_list.append([start,stop])

    return new_gene_list
    # # switch back to a list and sort it
    # genes = genes.tolist()
    # genes.sort()

    # # return the sorted list
    # return genes



#==============================================================
# Main program
#==============================================================
def main() -> None: # pragma: no cover
    """
    Main program function
    """
    # Gene detection over genome involves to consider a thymine instead of
    # an uracile that we would find on the expressed RNA
    #start_codons = ['TTG', 'CTG', 'ATT', 'ATG', 'GTG']
    #stop_codons = ['TAA', 'TAG', 'TGA']
    start_regex = re.compile('AT[TG]|[ATCG]TG')
    stop_regex = re.compile('TA[GA]|TGA')
    # Shine AGGAGGUAA
    #AGGA ou GGAGG 
    shine_regex = re.compile('A?G?GAGG|GGAG|GG.{1}GG')

    # Get the arguments
    args = get_arguments()

    # Read fasta file
    sequence = read_fasta(args.genome_file)

    # Forward gene search
    search_5_3 = predict_genes(sequence,
                              start_regex,
                              stop_regex,
                              shine_regex,
                              args.min_gene_len,
                              args.max_shine_dalgarno_distance,
                              args.min_gap)
    # Reverse gene search
    rv_sequence = reverse_complement(sequence)
    search_3_5 = predict_genes(rv_sequence,
                              start_regex,
                              stop_regex,
                              shine_regex,
                              args.min_gene_len,
                              args.max_shine_dalgarno_distance,
                              args.min_gap)

    search_3_5_rv = reorder_gene_list(search_3_5, sequence)

    probable_genes = search_5_3 + search_3_5_rv
    probable_genes.sort()

    write_genes_pos(args.predicted_genes_file, probable_genes)
    write_genes(args.fasta_file, sequence, search_5_3, rv_sequence, search_3_5)


if __name__ == '__main__':
    main()
