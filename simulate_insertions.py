#!/usr/bin/env python

import argparse
import sys

from pathlib import Path

from src.simulation import (get_sequences, 
                            insert_into_nucleus)

IDENTITIES = [(0.41, 0.5), (0.51, 0.60), 
              (0.61, 0.7), (0.71, 0.8),
              (0.81, 0.9), (0.91, 1)]


#Generating program options
def parse_arguments():
    desc = "Insert sequences with random mutations into a genome"
    parser = argparse.ArgumentParser(description=desc)
    
    
    help_insertion_source = "(Required) FASTA insertion source."
    parser.add_argument("--input_sequence", "-is", type=str,
                        help=help_insertion_source)

    help_insertion_positions = "(Required) FASTA insertion positions."
    parser.add_argument("--input_positions", "-ip", type=str,
                        help=help_insertion_positions)
    
    help_destiny_fasta = "(Required) FASTA destiny."
    parser.add_argument("--destiny", "-d", type=str,
                        help=help_destiny_fasta)

    help_insertion_length = "(Required) Insertion lengths."
    parser.add_argument("--insertion_length", "-l", type=int,
                        help=help_insertion_length)
    
    help_out_dir = "(Required) Out dir"
    parser.add_argument("--out_dir", "-o", type=str,
                        help=help_out_dir)
    if len(sys.argv)==1:
        parser.print_help()
        exit()
    return parser.parse_args()


#Parse and return values given to options when running this program
def get_arguments():
    parser = parse_arguments()
    input_sequence = parser.input_sequence
    input_positions = parser.input_positions
    destiny = parser.destiny
    insertion_length = parser.insertion_length
    out_dir = Path(parser.out_dir)
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    return {"input_sequence": input_sequence, "input_positions": input_positions,
            "destiny": destiny, "insertion_length": insertion_length, 
            "out_dir": out_dir}


def main():
    arguments = get_arguments()
    for indentity_range in IDENTITIES:
        freq1 = indentity_range[0]
        freq2 = indentity_range[1]
        out_dirname = "{}_{}".format(str(freq1), str(freq2))
        out_fpath = arguments["out_dir"] / out_dirname
        if not out_fpath.exists():
            out_fpath.mkdir(parents=True)
        with open(arguments["input_sequence"]) as input_fhand:
            with open(arguments["input_positions"]) as starting_fhand:
                sequences = get_sequences(input_fhand, starting_fhand, freq1, freq2, 
                                          arguments["insertion_length"], out_fpath)
        with open(arguments["destiny"]) as destiny_fhand:
            insertions = insert_into_nucleus(destiny_fhand, sequences, out_fpath)
        with open(out_fpath / "summary.txt", "w") as out_fhand:
            for line in insertions:
                write_ = "{}\t{}\t{}\n".format(str(line[0]), str(line[1].replace(">", "")), str(line[2]))
                out_fhand.write(write_)

if __name__ == "__main__":
    main()